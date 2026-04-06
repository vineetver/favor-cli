//! MetaSTAAR cross-study analysis and summary statistics export.

use std::collections::HashMap;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, AsArray, BooleanArray, BooleanBuilder, Float64Array, Float64Builder,
    Int32Array, Int32Builder, Int64Array, ListArray, ListBuilder, StringArray, StringBuilder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use faer::Mat;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use serde::{Deserialize, Serialize};

use super::carrier::sparse_score;
use super::carrier::AnalysisVectors;
use super::sparse_g::SparseG;
use super::masks::MaskGroup;
use super::score;
use super::GeneResult;
use crate::column::{Col, STAAR_WEIGHTS};
use crate::data::parquet_column_names;
use crate::engine::DfEngine;
use crate::error::FavorError;
use crate::ingest::{ColumnContract, ColumnRequirement};
use crate::output::Output;
use crate::types::{
    AnnotatedVariant, AnnotationWeights, Chromosome, Consequence, FunctionalAnnotation,
    MetaVariant, RegionType, RegulatoryFlags,
};

const META_SUMSTATS_COLUMNS: &[ColumnRequirement] = &[
    ColumnRequirement { name: "segment_id", source: "summary statistics", used_by: "segment identification" },
    ColumnRequirement { name: "positions", source: "summary statistics", used_by: "variant positions" },
    ColumnRequirement { name: "refs", source: "summary statistics", used_by: "reference alleles" },
    ColumnRequirement { name: "alts", source: "summary statistics", used_by: "alternate alleles" },
    ColumnRequirement { name: "cov_lower", source: "summary statistics", used_by: "covariance matrix" },
];

/// Load ALL segments for a (study, chromosome) by reading the parquet file directly.
pub fn load_all_segments(
    study: &StudyHandle,
    chrom: &str,
) -> Result<HashMap<i32, SegmentCov>, FavorError> {
    let seg_path = study.path.join(format!("chromosome={chrom}/segments.parquet"));
    if !seg_path.exists() { return Ok(HashMap::new()); }

    // Upfront schema validation
    let cols = parquet_column_names(&seg_path)?;
    let contract = ColumnContract { command: "meta-staar", required: META_SUMSTATS_COLUMNS };
    let missing = contract.check(&cols);
    if !missing.is_empty() {
        return Err(FavorError::Analysis(format!(
            "Invalid summary statistics in {}:\n{}\n\
             Re-export with: `favor staar --emit-sumstats`",
            seg_path.display(),
            ColumnContract::format_missing(&missing),
        )));
    }

    let file = File::open(&seg_path)
        .map_err(|e| FavorError::Analysis(format!("Cannot open {}: {e}", seg_path.display())))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| FavorError::Analysis(format!("Bad parquet {}: {e}", seg_path.display())))?
        .build()
        .map_err(|e| FavorError::Analysis(format!("Parquet reader error {}: {e}", seg_path.display())))?;

    let mut result = HashMap::new();
    for batch in reader {
        let batch = batch.map_err(|e| FavorError::Analysis(format!("Parquet read error: {e}")))?;
        let n = batch.num_rows();

        let seg_ids = batch.column_by_name("segment_id")
            .ok_or_else(|| FavorError::Analysis("Missing segment_id column".into()))?
            .as_any().downcast_ref::<Int32Array>()
            .ok_or_else(|| FavorError::Analysis("segment_id is not Int32".into()))?;
        let positions_col = batch.column_by_name("positions")
            .ok_or_else(|| FavorError::Analysis("Missing positions column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("positions is not List".into()))?;
        let refs_col = batch.column_by_name("refs")
            .ok_or_else(|| FavorError::Analysis("Missing refs column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("refs is not List".into()))?;
        let alts_col = batch.column_by_name("alts")
            .ok_or_else(|| FavorError::Analysis("Missing alts column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("alts is not List".into()))?;
        let cov_col = batch.column_by_name("cov_lower")
            .ok_or_else(|| FavorError::Analysis("Missing cov_lower column".into()))?
            .as_any().downcast_ref::<ListArray>()
            .ok_or_else(|| FavorError::Analysis("cov_lower is not List".into()))?;

        for row in 0..n {
            let seg_id = seg_ids.value(row);

            let pos_arr = positions_col.value(row);
            let pos_vals = pos_arr.as_any().downcast_ref::<Int32Array>()
                .ok_or_else(|| FavorError::Analysis("positions inner is not Int32".into()))?;
            let positions: Vec<u32> = (0..pos_vals.len()).map(|i| pos_vals.value(i) as u32).collect();

            let refs_arr = refs_col.value(row);
            let refs_str = refs_arr.as_string::<i32>();
            let refs: Vec<String> = (0..refs_str.len()).map(|i| refs_str.value(i).to_string()).collect();

            let alts_arr = alts_col.value(row);
            let alts_str = alts_arr.as_string::<i32>();
            let alts: Vec<String> = (0..alts_str.len()).map(|i| alts_str.value(i).to_string()).collect();

            let cov_arr = cov_col.value(row);
            let cov_vals = cov_arr.as_any().downcast_ref::<Float64Array>()
                .ok_or_else(|| FavorError::Analysis("cov_lower inner is not Float64".into()))?;
            let cov_lower: Vec<f64> = (0..cov_vals.len()).map(|i| cov_vals.value(i)).collect();

            let mut pos_index: HashMap<u32, Vec<usize>> = HashMap::new();
            for (i, &p) in positions.iter().enumerate() {
                pos_index.entry(p).or_default().push(i);
            }

            result.insert(seg_id, SegmentCov { refs, alts, cov_lower, pos_index });
        }
    }
    Ok(result)
}

pub struct StudyHandle {
    pub path: std::path::PathBuf,
    pub meta: StudyMeta,
}

/// Segment covariance read from one study.
pub struct SegmentCov {
    refs: Vec<String>,
    alts: Vec<String>,
    cov_lower: Vec<f64>,
    pos_index: HashMap<u32, Vec<usize>>,
}

impl SegmentCov {
    /// Extract sub-matrix for a subset of variants identified by (position, ref, alt).
    fn extract_submatrix(&self, keys: &[(u32, &str, &str)]) -> Mat<f64> {
        let m = keys.len();
        let mut mat = Mat::zeros(m, m);

        let mut key_to_local: Vec<Option<usize>> = Vec::with_capacity(m);
        for &(pos, ref_a, alt_a) in keys {
            let local = self.pos_index.get(&pos).and_then(|candidates| {
                candidates.iter().find(|&&i| self.refs[i] == ref_a && self.alts[i] == alt_a).copied()
            });
            key_to_local.push(local);
        }

        for i in 0..m {
            let Some(li) = key_to_local[i] else { continue };
            for j in 0..=i {
                let Some(lj) = key_to_local[j] else { continue };
                let (row, col) = if li >= lj { (li, lj) } else { (lj, li) };
                let idx = row * (row + 1) / 2 + col;
                if idx < self.cov_lower.len() {
                    mat[(i, j)] = self.cov_lower[idx];
                    mat[(j, i)] = self.cov_lower[idx];
                }
            }
        }
        mat
    }
}

/// Load and validate study directories.
pub fn load_studies(paths: &[std::path::PathBuf]) -> Result<Vec<StudyHandle>, FavorError> {
    let mut studies = Vec::with_capacity(paths.len());
    for path in paths {
        let meta_path = path.join("meta_staar.json");
        if !meta_path.exists() {
            return Err(FavorError::Input(format!(
                "Not a MetaSTAAR study directory: {}. Missing meta_staar.json. \
                 Run `favor staar --emit-sumstats` first.", path.display()
            )));
        }
        let content = std::fs::read_to_string(&meta_path)?;
        let meta: StudyMeta = serde_json::from_str(&content)
            .map_err(|e| FavorError::Input(format!("Invalid meta_staar.json in {}: {e}", path.display())))?;
        if meta.favor_meta_version != 1 {
            return Err(FavorError::Input(format!(
                "Unsupported meta version {} in {}. Expected 1.", meta.favor_meta_version, path.display()
            )));
        }
        studies.push(StudyHandle { path: path.clone(), meta });
    }

    if studies.len() < 2 {
        return Err(FavorError::Input("MetaSTAAR requires at least 2 studies.".into()));
    }

    let first_type = &studies[0].meta.trait_type;
    let first_seg = studies[0].meta.segment_size;
    for s in &studies[1..] {
        if s.meta.trait_type != *first_type {
            return Err(FavorError::Input(format!(
                "Trait type mismatch: {} vs {}. All studies must have the same trait type.",
                first_type, s.meta.trait_type
            )));
        }
        if s.meta.segment_size != first_seg {
            return Err(FavorError::Input(format!(
                "Segment size mismatch: {} vs {}. All studies must use the same segment size.",
                first_seg, s.meta.segment_size
            )));
        }
    }

    Ok(studies)
}

/// Build union variant catalog for one chromosome across all studies.
pub fn merge_chromosome(
    engine: &DfEngine,
    studies: &[StudyHandle],
    chrom: &str,
    maf_cutoff: f64,
) -> Result<Vec<MetaVariant>, FavorError> {
    let mut union_parts = Vec::new();
    for (idx, study) in studies.iter().enumerate() {
        let var_path = study.path.join(format!("chromosome={chrom}/variants.parquet"));
        if !var_path.exists() { continue; }
        engine.register_parquet_file(&format!("_study_{idx}"), &var_path)?;
        union_parts.push(format!("SELECT {idx} AS study_idx, * FROM _study_{idx}"));
    }
    if union_parts.is_empty() { return Ok(Vec::new()); }

    engine.execute(&format!(
        "CREATE OR REPLACE TABLE _study_variants AS {}",
        union_parts.join(" UNION ALL ")
    ))?;

    // Weight columns: first_value(w_col) AS w_col for each weight
    let weight_aggs: String = STAAR_WEIGHTS.iter()
        .map(|c| format!("first_value({col}) AS {col}", col = c.as_str()))
        .collect::<Vec<_>>()
        .join(", ");
    let weight_select: String = STAAR_WEIGHTS.iter()
        .map(|c| c.as_str())
        .collect::<Vec<_>>()
        .join(", ");

    engine.execute(&format!(
        "CREATE OR REPLACE TABLE _meta_variants AS \
         SELECT \
             {pos}, {ref_a}, {alt_a}, \
             SUM(u_stat) AS u_meta, \
             SUM(mac) AS mac_total, \
             SUM(n_obs) AS n_total, \
             first_value({gene}) FILTER (WHERE {gene} != '') AS {gene}, \
             first_value({region}) FILTER (WHERE {region} != '') AS {region}, \
             first_value({csq}) FILTER (WHERE {csq} != '') AS {csq}, \
             first_value({cadd}) AS {cadd}, \
             first_value({revel}) AS {revel}, \
             bool_or({cage_p}) AS {cage_p}, \
             bool_or({cage_e}) AS {cage_e}, \
             bool_or({ccre_p}) AS {ccre_p}, \
             bool_or({ccre_e}) AS {ccre_e}, \
             {weight_aggs}, \
             CAST(array_agg(named_struct('s', study_idx, 'seg', segment_id)) AS VARCHAR) AS study_segs \
         FROM _study_variants \
         WHERE {maf} < {maf_cutoff} \
         GROUP BY {pos}, {ref_a}, {alt_a} \
         ORDER BY {pos}",
        pos = Col::Position, ref_a = Col::RefAllele, alt_a = Col::AltAllele,
        maf = Col::Maf,
        gene = Col::GeneName, region = Col::RegionType, csq = Col::Consequence,
        cadd = Col::CaddPhred, revel = Col::Revel,
        cage_p = Col::IsCagePromoter, cage_e = Col::IsCageEnhancer,
        ccre_p = Col::IsCcrePromoter, ccre_e = Col::IsCcreEnhancer,
    ))?;

    let batches = engine.collect(
        &format!("SELECT {pos}, {ref_a}, {alt_a}, u_meta, \
         mac_total, n_total, {gene}, {region}, {csq}, \
         {cadd}, {revel}, \
         {cage_p}, {cage_e}, {ccre_p}, {ccre_e}, \
         {weight_select}, \
         study_segs \
         FROM _meta_variants ORDER BY {pos}",
        pos = Col::Position, ref_a = Col::RefAllele, alt_a = Col::AltAllele,
        gene = Col::GeneName, region = Col::RegionType, csq = Col::Consequence,
        cadd = Col::CaddPhred, revel = Col::Revel,
        cage_p = Col::IsCagePromoter, cage_e = Col::IsCageEnhancer,
        ccre_p = Col::IsCcrePromoter, ccre_e = Col::IsCcreEnhancer,
    ))?;

    let mut result = Vec::new();
    let chrom_parsed: Chromosome = chrom.parse().unwrap_or(Chromosome::Autosome(1));

    for batch in &batches {
        let n = batch.num_rows();
        let col = |i: usize| batch.column(i);
        let i32_col = |i: usize| col(i).as_any().downcast_ref::<Int32Array>()
            .ok_or_else(|| FavorError::Analysis(format!("_meta_variants col {i}: expected Int32")));
        let i64_col = |i: usize| col(i).as_any().downcast_ref::<Int64Array>()
            .ok_or_else(|| FavorError::Analysis(format!("_meta_variants col {i}: expected Int64")));
        let f64_col = |i: usize| col(i).as_any().downcast_ref::<Float64Array>()
            .ok_or_else(|| FavorError::Analysis(format!("_meta_variants col {i}: expected Float64")));
        let str_col = |i: usize| col(i).as_any().downcast_ref::<StringArray>()
            .ok_or_else(|| FavorError::Analysis(format!("_meta_variants col {i}: expected Utf8")));
        let bool_col = |i: usize| col(i).as_any().downcast_ref::<BooleanArray>()
            .ok_or_else(|| FavorError::Analysis(format!("_meta_variants col {i}: expected Boolean")));

        let pos_arr = i32_col(0)?;
        let ref_arr = str_col(1)?;
        let alt_arr = str_col(2)?;
        let u_meta_arr = f64_col(3)?;
        let mac_arr = i64_col(4)?;
        let n_obs_arr = i64_col(5)?;
        let gene_arr = str_col(6)?;
        let rt_arr = str_col(7)?;
        let csq_arr = str_col(8)?;
        let cadd_arr = f64_col(9)?;
        let revel_arr = f64_col(10)?;
        let cp_arr = bool_col(11)?;
        let ce_arr = bool_col(12)?;
        let crp_arr = bool_col(13)?;
        let cre_arr = bool_col(14)?;

        let mut w_arrs: Vec<&Float64Array> = Vec::with_capacity(11);
        for i in 0..11 {
            w_arrs.push(f64_col(15 + i)?);
        }
        let segs_arr = str_col(26)?;

        for i in 0..n {
            let mut weights = [0.0f64; 11];
            for (ch, wa) in w_arrs.iter().enumerate() {
                weights[ch] = if wa.is_null(i) { 0.0 } else { wa.value(i) };
            }

            let segs_str = if segs_arr.is_null(i) { "" } else { segs_arr.value(i) };
            let study_segments = parse_study_segments(segs_str);

            let mac_total = if mac_arr.is_null(i) { 0 } else { mac_arr.value(i) };
            let n_total = if n_obs_arr.is_null(i) { 0 } else { n_obs_arr.value(i) };
            // MAF = MAC / (2*N) — diploid genomes have 2 alleles per sample
            let maf = if n_total > 0 { mac_total as f64 / (2.0 * n_total as f64) } else { 0.0 };

            result.push(MetaVariant {
                variant: AnnotatedVariant {
                    chromosome: chrom_parsed,
                    position: pos_arr.value(i) as u32,
                    ref_allele: if ref_arr.is_null(i) { "".into() } else { ref_arr.value(i).into() },
                    alt_allele: if alt_arr.is_null(i) { "".into() } else { alt_arr.value(i).into() },
                    maf,
                    gene_name: if gene_arr.is_null(i) { "".into() } else { gene_arr.value(i).into() },
                    annotation: FunctionalAnnotation {
                        region_type: RegionType::from_str_lossy(if rt_arr.is_null(i) { "" } else { rt_arr.value(i) }),
                        consequence: Consequence::from_str_lossy(if csq_arr.is_null(i) { "" } else { csq_arr.value(i) }),
                        cadd_phred: if cadd_arr.is_null(i) { 0.0 } else { cadd_arr.value(i) },
                        revel: if revel_arr.is_null(i) { 0.0 } else { revel_arr.value(i) },
                        regulatory: RegulatoryFlags {
                            cage_promoter: if cp_arr.is_null(i) { false } else { cp_arr.value(i) },
                            cage_enhancer: if ce_arr.is_null(i) { false } else { ce_arr.value(i) },
                            ccre_promoter: if crp_arr.is_null(i) { false } else { crp_arr.value(i) },
                            ccre_enhancer: if cre_arr.is_null(i) { false } else { cre_arr.value(i) },
                        },
                        weights: AnnotationWeights(weights),
                    },
                },
                u_meta: if u_meta_arr.is_null(i) { 0.0 } else { u_meta_arr.value(i) },
                mac_total,
                n_total,
                study_segments,
            });
        }
    }

    engine.execute("DROP TABLE IF EXISTS _study_variants")?;
    engine.execute("DROP TABLE IF EXISTS _meta_variants")?;

    Ok(result)
}

/// Run meta-analysis score tests for one gene/mask group.
pub fn meta_score_gene(
    group: &MaskGroup,
    meta_variants: &[MetaVariant],
    studies: &[StudyHandle],
    _chrom: &str,
    segment_cache: &HashMap<(usize, i32), SegmentCov>,
) -> Option<GeneResult> {
    let indices: Vec<usize> = group.variant_indices.iter()
        .filter(|&&i| i < meta_variants.len())
        .copied()
        .collect();
    if indices.len() < 2 { return None; }

    let m = indices.len();

    let mut u = Mat::zeros(m, 1);
    for (local, &gi) in indices.iter().enumerate() {
        u[(local, 0)] = meta_variants[gi].u_meta;
    }

    let keys: Vec<(u32, &str, &str)> = indices.iter()
        .map(|&gi| (meta_variants[gi].variant.position, &*meta_variants[gi].variant.ref_allele, &*meta_variants[gi].variant.alt_allele))
        .collect();

    let mut cov = Mat::zeros(m, m);
    for study_idx in 0..studies.len() {
        let mut needed_segments: std::collections::HashSet<i32> = std::collections::HashSet::new();
        for &gi in &indices {
            for &(sidx, seg_id) in &meta_variants[gi].study_segments {
                if sidx == study_idx {
                    needed_segments.insert(seg_id);
                }
            }
        }

        for seg_id in needed_segments {
            let cache_key = (study_idx, seg_id);
            if let Some(seg) = segment_cache.get(&cache_key) {
                let sub = seg.extract_submatrix(&keys);
                for i in 0..m {
                    for j in 0..m {
                        cov[(i, j)] += sub[(i, j)];
                    }
                }
            }
        }
    }

    let mafs: Vec<f64> = indices.iter().map(|&gi| {
        let mv = &meta_variants[gi];
        if mv.n_total > 0 { mv.mac_total as f64 / (2.0 * mv.n_total as f64) } else { 0.0 }
    }).collect();

    // Use max N across variants in this gene for MAC-based ACAT-V grouping
    let n_total: usize = indices.iter()
        .map(|&gi| meta_variants[gi].n_total as usize)
        .max()
        .unwrap_or(0);

    let ann_matrix: Vec<Vec<f64>> = (0..11).map(|ch| {
        indices.iter().map(|&gi| {
            meta_variants[gi].variant.annotation.weights.0[ch]
        }).collect()
    }).collect();

    let sr = score::run_staar_from_sumstats(&u, &cov, &ann_matrix, &mafs, n_total);

    let cmac: i64 = indices.iter().map(|&gi| meta_variants[gi].mac_total).sum();

    Some(GeneResult {
        ensembl_id: group.name.clone(),
        gene_symbol: group.name.clone(),
        chromosome: group.chromosome,
        start: group.start,
        end: group.end,
        n_variants: m as u32,
        cumulative_mac: cmac as u32,
        staar: sr,
    })
}

fn parse_study_segments(s: &str) -> Vec<(usize, i32)> {
    let mut result = Vec::new();
    for part in s.split('{') {
        let part = part.trim_matches(|c: char| c == '[' || c == ']' || c == ',' || c == ' ' || c == '}');
        if part.is_empty() { continue; }
        let mut study = None;
        let mut seg = None;
        for kv in part.split(',') {
            let kv = kv.trim();
            if let Some(val) = kv.strip_prefix("s:").or_else(|| kv.strip_prefix("s :")) {
                study = val.trim().parse().ok();
            } else if let Some(val) = kv.strip_prefix("seg:").or_else(|| kv.strip_prefix("seg :")) {
                seg = val.trim().parse().ok();
            }
        }
        if let (Some(s), Some(g)) = (study, seg) {
            result.push((s, g));
        }
    }
    result
}

// ──────────────────────────────────────────────────────────────────────────────
// Summary statistics export (formerly sumstats.rs)
// ──────────────────────────────────────────────────────────────────────────────

const SEGMENT_BP: u32 = 500_000;
const MAX_SEGMENT_VARIANTS: usize = 2000;

#[derive(Serialize, Deserialize)]
pub struct StudyMeta {
    pub favor_meta_version: u32,
    pub trait_type: String,
    pub trait_name: String,
    pub n_samples: usize,
    pub sigma2: f64,
    pub maf_cutoff: f64,
    pub covariates: Vec<String>,
    pub segment_size: u32,
}

/// Export MetaSTAAR summary statistics using carrier-indexed sparse scoring.
///
/// Computes U and K per segment directly from carrier lists in O(total_MAC)
/// instead of building dense genotype matrices.
pub fn emit_sumstats(
    store_dir: &Path,
    analysis: &AnalysisVectors,
    variants: &[AnnotatedVariant],
    output_dir: &Path,
    meta: &StudyMeta,
    out: &dyn Output,
) -> Result<(), FavorError> {
    out.status("MetaSTAAR: computing summary statistics (carrier-indexed)...");

    let chrom_set: Vec<String> = variants
        .iter()
        .map(|v| v.chromosome.label())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();

    for chrom in &chrom_set {
        let indices: Vec<usize> = variants
            .iter()
            .enumerate()
            .filter(|(_, v)| v.chromosome.label() == *chrom)
            .map(|(i, _)| i)
            .collect();
        if indices.is_empty() {
            continue;
        }

        let dir = output_dir.join(format!("chromosome={chrom}"));
        std::fs::create_dir_all(&dir)
            .map_err(|e| FavorError::Resource(format!("Cannot create '{}': {e}", dir.display())))?;

        let chrom_dir = store_dir.join(format!("chromosome={chrom}"));
        let sparse_g = SparseG::open(&chrom_dir)?;
        let variant_index = super::carrier::VariantIndex::load(&chrom_dir)?;
        emit_chromosome_sparse(&sparse_g, &variant_index, analysis, variants, &indices, chrom, &dir, out)?;
    }

    let meta_json = serde_json::to_string_pretty(meta)
        .map_err(|e| FavorError::Resource(format!("Failed to serialize meta_staar.json: {e}")))?;
    let meta_path = output_dir.join("meta_staar.json");
    std::fs::write(&meta_path, meta_json)
        .map_err(|e| FavorError::Resource(format!("Cannot write '{}': {e}", meta_path.display())))?;

    out.success(&format!("Summary statistics -> {}", output_dir.display()));
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn emit_chromosome_sparse(
    sparse_g: &SparseG,
    variant_index: &super::carrier::VariantIndex,
    analysis: &AnalysisVectors,
    variants: &[AnnotatedVariant],
    chrom_indices: &[usize],
    chrom: &str,
    dir: &Path,
    out: &dyn Output,
) -> Result<(), FavorError> {
    let n = analysis.n_pheno;

    // Build position-based segments (same binning as v1 for output compat)
    let mut coarse: std::collections::BTreeMap<i32, Vec<usize>> =
        std::collections::BTreeMap::new();
    for &gi in chrom_indices {
        coarse
            .entry((variants[gi].position / SEGMENT_BP) as i32)
            .or_default()
            .push(gi);
    }
    let mut segments: Vec<(i32, Vec<usize>)> = Vec::new();
    let mut next_id = 0i32;
    for (_, indices) in coarse {
        for chunk in indices.chunks(MAX_SEGMENT_VARIANTS) {
            segments.push((next_id, chunk.to_vec()));
            next_id += 1;
        }
    }
    out.status(&format!(
        "  chr{chrom}: {} variants, {} segments (sparse)",
        chrom_indices.len(),
        segments.len()
    ));

    let total_variants = chrom_indices.len();
    let mut b_position = Int32Builder::with_capacity(total_variants);
    let mut b_ref = StringBuilder::with_capacity(total_variants, total_variants * 4);
    let mut b_alt = StringBuilder::with_capacity(total_variants, total_variants * 4);
    let mut b_maf = Float64Builder::with_capacity(total_variants);
    let mut b_mac = Int32Builder::with_capacity(total_variants);
    let mut b_n_obs = Int32Builder::with_capacity(total_variants);
    let mut b_u_stat = Float64Builder::with_capacity(total_variants);
    let mut b_v_stat = Float64Builder::with_capacity(total_variants);
    let mut b_segment_id = Int32Builder::with_capacity(total_variants);
    let mut b_gene = StringBuilder::with_capacity(total_variants, total_variants * 8);
    let mut b_region = StringBuilder::with_capacity(total_variants, total_variants * 8);
    let mut b_consequence = StringBuilder::with_capacity(total_variants, total_variants * 8);
    let mut b_cadd = Float64Builder::with_capacity(total_variants);
    let mut b_revel = Float64Builder::with_capacity(total_variants);
    let mut b_cage_prom = BooleanBuilder::with_capacity(total_variants);
    let mut b_cage_enh = BooleanBuilder::with_capacity(total_variants);
    let mut b_ccre_prom = BooleanBuilder::with_capacity(total_variants);
    let mut b_ccre_enh = BooleanBuilder::with_capacity(total_variants);
    let mut b_weights: [Float64Builder; 11] =
        std::array::from_fn(|_| Float64Builder::with_capacity(total_variants));

    let n_segments = segments.len();
    let mut s_segment_id = Int32Builder::with_capacity(n_segments);
    let mut s_n_variants = Int32Builder::with_capacity(n_segments);
    let mut s_positions = ListBuilder::new(Int32Builder::with_capacity(total_variants));
    let mut s_refs = ListBuilder::new(StringBuilder::with_capacity(
        total_variants,
        total_variants * 4,
    ));
    let mut s_alts = ListBuilder::new(StringBuilder::with_capacity(
        total_variants,
        total_variants * 4,
    ));
    let mut s_cov_lower = ListBuilder::new(Float64Builder::with_capacity(
        total_variants * (total_variants + 1) / 2,
    ));

    for (seg_id, seg_indices) in &segments {
        let seg_id = *seg_id;
        let m = seg_indices.len();

        // Resolve each segment variant to variant_vcf via VariantIndex, then
        // batch-load carrier lists from SparseG.
        let variant_vcfs: Vec<u32> = seg_indices
            .iter()
            .filter_map(|&gi| {
                let v = &variants[gi];
                let vid = crate::types::format_vid(chrom, v.position, &v.ref_allele, &v.alt_allele);
                variant_index.resolve_vid(&vid)
            })
            .collect();
        let seg_carriers = sparse_g.load_variants(&variant_vcfs);

        let (u, k) = sparse_score::score_gene_sparse(&seg_carriers, analysis);

        for (j, &gi) in seg_indices.iter().enumerate() {
            let v = &variants[gi];
            b_position.append_value(v.position as i32);
            b_ref.append_value(&v.ref_allele);
            b_alt.append_value(&v.alt_allele);
            b_maf.append_value(v.maf);
            b_mac.append_value((2.0 * v.maf * n as f64).round() as i32);
            b_n_obs.append_value(n as i32);
            b_u_stat.append_value(u[(j, 0)]);
            b_v_stat.append_value(k[(j, j)]);
            b_segment_id.append_value(seg_id);
            b_gene.append_value(&v.gene_name);
            b_region.append_value(v.annotation.region_type.as_str());
            b_consequence.append_value(v.annotation.consequence.as_str());
            b_cadd.append_value(v.annotation.cadd_phred);
            b_revel.append_value(v.annotation.revel);
            b_cage_prom.append_value(v.annotation.regulatory.cage_promoter);
            b_cage_enh.append_value(v.annotation.regulatory.cage_enhancer);
            b_ccre_prom.append_value(v.annotation.regulatory.ccre_promoter);
            b_ccre_enh.append_value(v.annotation.regulatory.ccre_enhancer);
            for (i, builder) in b_weights.iter_mut().enumerate() {
                builder.append_value(v.annotation.weights.0[i]);
            }
        }

        s_segment_id.append_value(seg_id);
        s_n_variants.append_value(m as i32);

        let pos_builder = s_positions.values();
        for &gi in seg_indices {
            pos_builder.append_value(variants[gi].position as i32);
        }
        s_positions.append(true);

        let ref_builder = s_refs.values();
        for &gi in seg_indices {
            ref_builder.append_value(&variants[gi].ref_allele);
        }
        s_refs.append(true);

        let alt_builder = s_alts.values();
        for &gi in seg_indices {
            alt_builder.append_value(&variants[gi].alt_allele);
        }
        s_alts.append(true);

        let cov_builder = s_cov_lower.values();
        for i in 0..m {
            for j in 0..=i {
                cov_builder.append_value(k[(i, j)]);
            }
        }
        s_cov_lower.append(true);
    }

    write_variants_parquet(
        dir,
        &mut b_position,
        &mut b_ref,
        &mut b_alt,
        &mut b_maf,
        &mut b_mac,
        &mut b_n_obs,
        &mut b_u_stat,
        &mut b_v_stat,
        &mut b_segment_id,
        &mut b_gene,
        &mut b_region,
        &mut b_consequence,
        &mut b_cadd,
        &mut b_revel,
        &mut b_cage_prom,
        &mut b_cage_enh,
        &mut b_ccre_prom,
        &mut b_ccre_enh,
        &mut b_weights,
    )?;

    write_segments_parquet(
        dir,
        &mut s_segment_id,
        &mut s_n_variants,
        &mut s_positions,
        &mut s_refs,
        &mut s_alts,
        &mut s_cov_lower,
    )?;

    out.status(&format!("    chr{chrom} done"));
    Ok(())
}

fn variant_schema() -> Schema {
    let mut fields = vec![
        Field::new(Col::Position.as_str(), DataType::Int32, false),
        Field::new(Col::RefAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::AltAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::Maf.as_str(), DataType::Float64, false),
        Field::new("mac", DataType::Int32, false),
        Field::new("n_obs", DataType::Int32, false),
        Field::new("u_stat", DataType::Float64, false),
        Field::new("v_stat", DataType::Float64, false),
        Field::new("segment_id", DataType::Int32, false),
        Field::new(Col::GeneName.as_str(), DataType::Utf8, false),
        Field::new(Col::RegionType.as_str(), DataType::Utf8, false),
        Field::new(Col::Consequence.as_str(), DataType::Utf8, false),
        Field::new(Col::CaddPhred.as_str(), DataType::Float64, false),
        Field::new(Col::Revel.as_str(), DataType::Float64, false),
        Field::new(Col::IsCagePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCageEnhancer.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcrePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcreEnhancer.as_str(), DataType::Boolean, false),
    ];
    for col in &STAAR_WEIGHTS {
        fields.push(Field::new(col.as_str(), DataType::Float64, false));
    }
    Schema::new(fields)
}

fn segment_schema() -> Schema {
    Schema::new(vec![
        Field::new("segment_id", DataType::Int32, false),
        Field::new("n_variants", DataType::Int32, false),
        Field::new(
            "positions",
            DataType::List(Arc::new(Field::new("item", DataType::Int32, true))),
            false,
        ),
        Field::new(
            "refs",
            DataType::List(Arc::new(Field::new("item", DataType::Utf8, true))),
            false,
        ),
        Field::new(
            "alts",
            DataType::List(Arc::new(Field::new("item", DataType::Utf8, true))),
            false,
        ),
        Field::new(
            "cov_lower",
            DataType::List(Arc::new(Field::new("item", DataType::Float64, true))),
            false,
        ),
    ])
}

fn parquet_props() -> WriterProperties {
    WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_row_count(Some(5000))
        .build()
}

#[allow(clippy::too_many_arguments)]
fn write_variants_parquet(
    dir: &Path,
    b_position: &mut Int32Builder,
    b_ref: &mut StringBuilder,
    b_alt: &mut StringBuilder,
    b_maf: &mut Float64Builder,
    b_mac: &mut Int32Builder,
    b_n_obs: &mut Int32Builder,
    b_u_stat: &mut Float64Builder,
    b_v_stat: &mut Float64Builder,
    b_segment_id: &mut Int32Builder,
    b_gene: &mut StringBuilder,
    b_region: &mut StringBuilder,
    b_consequence: &mut StringBuilder,
    b_cadd: &mut Float64Builder,
    b_revel: &mut Float64Builder,
    b_cage_prom: &mut BooleanBuilder,
    b_cage_enh: &mut BooleanBuilder,
    b_ccre_prom: &mut BooleanBuilder,
    b_ccre_enh: &mut BooleanBuilder,
    b_weights: &mut [Float64Builder; 11],
) -> Result<(), FavorError> {
    let schema = Arc::new(variant_schema());
    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(b_position.finish()),
        Arc::new(b_ref.finish()),
        Arc::new(b_alt.finish()),
        Arc::new(b_maf.finish()),
        Arc::new(b_mac.finish()),
        Arc::new(b_n_obs.finish()),
        Arc::new(b_u_stat.finish()),
        Arc::new(b_v_stat.finish()),
        Arc::new(b_segment_id.finish()),
        Arc::new(b_gene.finish()),
        Arc::new(b_region.finish()),
        Arc::new(b_consequence.finish()),
        Arc::new(b_cadd.finish()),
        Arc::new(b_revel.finish()),
        Arc::new(b_cage_prom.finish()),
        Arc::new(b_cage_enh.finish()),
        Arc::new(b_ccre_prom.finish()),
        Arc::new(b_ccre_enh.finish()),
    ];
    for b in b_weights.iter_mut() {
        columns.push(Arc::new(b.finish()));
    }

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;

    let file = File::create(dir.join("variants.parquet"))
        .map_err(|e| FavorError::Resource(format!("Create variants.parquet: {e}")))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(parquet_props()))
        .map_err(|e| FavorError::Resource(format!("Parquet writer init: {e}")))?;
    writer
        .write(&batch)
        .map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
    writer
        .close()
        .map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

    Ok(())
}

fn write_segments_parquet(
    dir: &Path,
    s_segment_id: &mut Int32Builder,
    s_n_variants: &mut Int32Builder,
    s_positions: &mut ListBuilder<Int32Builder>,
    s_refs: &mut ListBuilder<StringBuilder>,
    s_alts: &mut ListBuilder<StringBuilder>,
    s_cov_lower: &mut ListBuilder<Float64Builder>,
) -> Result<(), FavorError> {
    let schema = Arc::new(segment_schema());
    let columns: Vec<ArrayRef> = vec![
        Arc::new(s_segment_id.finish()),
        Arc::new(s_n_variants.finish()),
        Arc::new(s_positions.finish()),
        Arc::new(s_refs.finish()),
        Arc::new(s_alts.finish()),
        Arc::new(s_cov_lower.finish()),
    ];

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;

    let file = File::create(dir.join("segments.parquet"))
        .map_err(|e| FavorError::Resource(format!("Create segments.parquet: {e}")))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(parquet_props()))
        .map_err(|e| FavorError::Resource(format!("Parquet writer init: {e}")))?;
    writer
        .write(&batch)
        .map_err(|e| FavorError::Resource(format!("Parquet write: {e}")))?;
    writer
        .close()
        .map_err(|e| FavorError::Resource(format!("Parquet close: {e}")))?;

    Ok(())
}
