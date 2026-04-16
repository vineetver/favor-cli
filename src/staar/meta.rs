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
use super::masks::MaskGroup;
use super::score;
use super::GeneResult;

/// Meta-analysis per-gene result. Wraps a `GeneResult` with the
/// unweighted burden effect size — only meaningful for meta runs, so it
/// lives here instead of bloating `GeneResult` with NaN-80%-of-the-time
/// fields on the single-study path.
#[derive(Debug, Clone)]
pub struct MetaGeneResult {
    pub gene: GeneResult,
    /// Unweighted burden coefficient β̂ = (1ᵀU)/(1ᵀK1).
    pub burden_beta: f64,
    /// SE of β̂: sqrt(1 / 1ᵀK1).
    pub burden_se: f64,
}
use crate::column::{Col, STAAR_PHRED_CHANNELS};
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::ingest::{ColumnContract, ColumnRequirement};
use crate::output::Output;
use crate::runtime::Engine;
use crate::store::cohort::{ChromosomeView, CohortId};
use crate::store::list::parquet_column_names;
use crate::types::{
    AnnotatedVariant, AnnotationWeights, Chromosome, Consequence, FunctionalAnnotation,
    MetaVariant, RegionType, RegulatoryFlags,
};

const META_SUMSTATS_COLUMNS: &[ColumnRequirement] = &[
    ColumnRequirement {
        name: "segment_id",
        source: "summary statistics",
        used_by: "segment identification",
    },
    ColumnRequirement {
        name: "positions",
        source: "summary statistics",
        used_by: "variant positions",
    },
    ColumnRequirement {
        name: "refs",
        source: "summary statistics",
        used_by: "reference alleles",
    },
    ColumnRequirement {
        name: "alts",
        source: "summary statistics",
        used_by: "alternate alleles",
    },
    ColumnRequirement {
        name: "cov_lower",
        source: "summary statistics",
        used_by: "covariance matrix",
    },
];

/// Load ALL segments for a (study, chromosome) by reading the parquet file directly.
pub fn load_all_segments(
    study: &StudyHandle,
    chrom: &str,
) -> Result<HashMap<i32, SegmentCov>, CohortError> {
    let seg_path = study
        .path
        .join(format!("chromosome={chrom}/segments.parquet"));
    if !seg_path.exists() {
        return Ok(HashMap::new());
    }

    // Upfront schema validation
    let cols = parquet_column_names(&seg_path)?;
    let contract = ColumnContract {
        command: "meta-staar",
        required: META_SUMSTATS_COLUMNS,
    };
    let missing = contract.check(&cols);
    if !missing.is_empty() {
        return Err(CohortError::Analysis(format!(
            "Invalid summary statistics in {}:\n{}\n\
             Re-export with: `favorstaar --emit-sumstats`",
            seg_path.display(),
            ColumnContract::format_missing(&missing),
        )));
    }

    let file = File::open(&seg_path)
        .map_err(|e| CohortError::Analysis(format!("Cannot open {}: {e}", seg_path.display())))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| CohortError::Analysis(format!("Bad parquet {}: {e}", seg_path.display())))?
        .build()
        .map_err(|e| {
            CohortError::Analysis(format!("Parquet reader error {}: {e}", seg_path.display()))
        })?;

    let mut result = HashMap::new();
    for batch in reader {
        let batch = batch.map_err(|e| CohortError::Analysis(format!("Parquet read error: {e}")))?;
        let n = batch.num_rows();

        let seg_ids = batch
            .column_by_name("segment_id")
            .ok_or_else(|| CohortError::Analysis("Missing segment_id column".into()))?
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or_else(|| CohortError::Analysis("segment_id is not Int32".into()))?;
        let positions_col = batch
            .column_by_name("positions")
            .ok_or_else(|| CohortError::Analysis("Missing positions column".into()))?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| CohortError::Analysis("positions is not List".into()))?;
        let refs_col = batch
            .column_by_name("refs")
            .ok_or_else(|| CohortError::Analysis("Missing refs column".into()))?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| CohortError::Analysis("refs is not List".into()))?;
        let alts_col = batch
            .column_by_name("alts")
            .ok_or_else(|| CohortError::Analysis("Missing alts column".into()))?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| CohortError::Analysis("alts is not List".into()))?;
        let cov_col = batch
            .column_by_name("cov_lower")
            .ok_or_else(|| CohortError::Analysis("Missing cov_lower column".into()))?
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| CohortError::Analysis("cov_lower is not List".into()))?;

        for row in 0..n {
            let seg_id = seg_ids.value(row);

            let pos_arr = positions_col.value(row);
            let pos_vals = pos_arr
                .as_any()
                .downcast_ref::<Int32Array>()
                .ok_or_else(|| CohortError::Analysis("positions inner is not Int32".into()))?;
            let positions: Vec<u32> = (0..pos_vals.len())
                .map(|i| pos_vals.value(i) as u32)
                .collect();

            let refs_arr = refs_col.value(row);
            let refs_str = refs_arr.as_string::<i32>();
            let refs: Vec<String> = (0..refs_str.len())
                .map(|i| refs_str.value(i).to_string())
                .collect();

            let alts_arr = alts_col.value(row);
            let alts_str = alts_arr.as_string::<i32>();
            let alts: Vec<String> = (0..alts_str.len())
                .map(|i| alts_str.value(i).to_string())
                .collect();

            let cov_arr = cov_col.value(row);
            let cov_vals = cov_arr
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| CohortError::Analysis("cov_lower inner is not Float64".into()))?;
            let cov_lower: Vec<f64> = (0..cov_vals.len()).map(|i| cov_vals.value(i)).collect();

            let mut pos_index: HashMap<u32, Vec<usize>> = HashMap::new();
            for (i, &p) in positions.iter().enumerate() {
                pos_index.entry(p).or_default().push(i);
            }

            result.insert(
                seg_id,
                SegmentCov {
                    refs,
                    alts,
                    cov_lower,
                    pos_index,
                },
            );
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

/// Lex-min orientation of a biallelic site. Two studies that disagree on
/// which allele is REF vs ALT must hash to the same `(pos, ref, alt)` group
/// so their summary statistics can be combined.
fn canonical_alleles<'a>(ref_b: &'a str, alt_b: &'a str) -> (&'a str, &'a str) {
    if ref_b <= alt_b {
        (ref_b, alt_b)
    } else {
        (alt_b, ref_b)
    }
}

impl SegmentCov {
    /// Extract sub-matrix for a subset of variants identified by `(pos, ref,
    /// alt)` in the **caller's canonical orientation**.
    ///
    /// Lookup is orientation-blind so a study that stored a variant flipped
    /// still resolves. The caller's K is signed against canonical-orientation
    /// dosage; the stored covariance is signed against local-orientation
    /// dosage. With `G_local = ±G_canonical`, the j-th diagonal of K is
    /// invariant (`G²` = `(−G)²`) but each off-diagonal `K[i,j] = G_i·G_j`
    /// flips sign whenever exactly one of `(i, j)` is locally flipped. We
    /// track a per-row sign and multiply through, so the returned sub-matrix
    /// is consistently in canonical orientation.
    fn extract_submatrix(&self, keys: &[(u32, &str, &str)]) -> Mat<f64> {
        let m = keys.len();
        let mut mat = Mat::zeros(m, m);

        let mut key_to_local: Vec<Option<(usize, f64)>> = Vec::with_capacity(m);
        for &(pos, ref_a, alt_a) in keys {
            let canon = canonical_alleles(ref_a, alt_a);
            let resolved = self.pos_index.get(&pos).and_then(|candidates| {
                candidates.iter().find_map(|&i| {
                    let local = (self.refs[i].as_str(), self.alts[i].as_str());
                    if canonical_alleles(local.0, local.1) != canon {
                        return None;
                    }
                    // Local matches canonical when its REF is the lex-min,
                    // i.e. local.0 == canon.0. Otherwise dosage is negated.
                    let sign = if local.0 == canon.0 { 1.0 } else { -1.0 };
                    Some((i, sign))
                })
            });
            key_to_local.push(resolved);
        }

        for i in 0..m {
            let Some((li, si)) = key_to_local[i] else {
                continue;
            };
            for j in 0..=i {
                let Some((lj, sj)) = key_to_local[j] else {
                    continue;
                };
                let (row, col) = if li >= lj { (li, lj) } else { (lj, li) };
                let idx = row * (row + 1) / 2 + col;
                if idx < self.cov_lower.len() {
                    let v = self.cov_lower[idx] * si * sj;
                    mat[(i, j)] = v;
                    mat[(j, i)] = v;
                }
            }
        }
        mat
    }
}

/// Load and validate study directories.
pub fn load_studies(paths: &[std::path::PathBuf]) -> Result<Vec<StudyHandle>, CohortError> {
    let mut studies = Vec::with_capacity(paths.len());
    for path in paths {
        let meta_path = path.join("meta_staar.json");
        if !meta_path.exists() {
            return Err(CohortError::Input(format!(
                "Not a MetaSTAAR study directory: {}. Missing meta_staar.json. \
                 Run `favorstaar --emit-sumstats` first.",
                path.display()
            )));
        }
        let content = std::fs::read_to_string(&meta_path)?;
        let meta: StudyMeta = serde_json::from_str(&content).map_err(|e| {
            CohortError::Input(format!(
                "Invalid meta_staar.json in {}: {e}",
                path.display()
            ))
        })?;
        if meta.cohort_meta_version != 1 {
            return Err(CohortError::Input(format!(
                "Unsupported meta version {} in {}. Expected 1.",
                meta.cohort_meta_version,
                path.display()
            )));
        }
        studies.push(StudyHandle {
            path: path.clone(),
            meta,
        });
    }

    if studies.is_empty() {
        return Err(CohortError::Input(
            "MetaSTAAR requires at least 1 study directory.".into(),
        ));
    }

    let first_type = &studies[0].meta.trait_type;
    let first_seg = studies[0].meta.segment_size;
    for s in &studies[1..] {
        if s.meta.trait_type != *first_type {
            return Err(CohortError::Input(format!(
                "Trait type mismatch: {} vs {}. All studies must have the same trait type.",
                first_type, s.meta.trait_type
            )));
        }
        if s.meta.segment_size != first_seg {
            return Err(CohortError::Input(format!(
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
) -> Result<Vec<MetaVariant>, CohortError> {
    let mut union_parts = Vec::new();
    for (idx, study) in studies.iter().enumerate() {
        let var_path = study
            .path
            .join(format!("chromosome={chrom}/variants.parquet"));
        if !var_path.exists() {
            continue;
        }
        engine.register_parquet_file(&format!("_study_{idx}"), &var_path)?;
        union_parts.push(format!("SELECT {idx} AS study_idx, * FROM _study_{idx}"));
    }
    if union_parts.is_empty() {
        return Ok(Vec::new());
    }

    engine.execute(&format!(
        "CREATE OR REPLACE TABLE _study_variants AS {}",
        union_parts.join(" UNION ALL ")
    ))?;

    // Weight columns: first_value(w_col) AS w_col for each weight
    let weight_aggs: String = STAAR_PHRED_CHANNELS
        .iter()
        .map(|c| format!("first_value({col}) AS {col}", col = c.as_str()))
        .collect::<Vec<_>>()
        .join(", ");
    let weight_select: String = STAAR_PHRED_CHANNELS
        .iter()
        .map(|c| c.as_str())
        .collect::<Vec<_>>()
        .join(", ");

    // Group by lex-min(ref,alt) so two studies that disagree on REF/ALT
    // orientation collapse into one row. The score statistic is signed by
    // alt-allele dosage, so flipping ref↔alt flips the sign of `u_stat`.
    // The covariance K is invariant under that flip and stays put — the
    // segment-cache lookup canonicalises on its own side.
    // `cadd_phred_raw` carries the numeric PHRED used by mask predicates;
    // the transformed [0,1] `cadd_phred` weight arrives via {weight_aggs}
    // as channel 0. Keeping them under distinct names avoids a duplicate
    // column error from DuckDB.
    engine.execute(&format!(
        "CREATE OR REPLACE TABLE _meta_variants AS \
         SELECT \
             {pos}, \
             CASE WHEN {ref_a} <= {alt_a} THEN {ref_a} ELSE {alt_a} END AS {ref_a}, \
             CASE WHEN {ref_a} <= {alt_a} THEN {alt_a} ELSE {ref_a} END AS {alt_a}, \
             SUM(CASE WHEN {ref_a} <= {alt_a} THEN u_stat ELSE -u_stat END) AS u_meta, \
             SUM(mac) AS mac_total, \
             SUM(n_obs) AS n_total, \
             first_value({gene}) FILTER (WHERE {gene} != '') AS {gene}, \
             first_value({region}) FILTER (WHERE {region} != '') AS {region}, \
             first_value({csq}) FILTER (WHERE {csq} != '') AS {csq}, \
             first_value(cadd_phred_raw) AS cadd_phred_raw, \
             first_value({revel}) AS {revel}, \
             bool_or({cage_p}) AS {cage_p}, \
             bool_or({cage_e}) AS {cage_e}, \
             bool_or({ccre_p}) AS {ccre_p}, \
             bool_or({ccre_e}) AS {ccre_e}, \
             {weight_aggs}, \
             CAST(array_agg(named_struct('s', study_idx, 'seg', segment_id)) AS VARCHAR) AS study_segs \
         FROM _study_variants \
         WHERE {maf} < {maf_cutoff} \
         GROUP BY {pos}, \
             CASE WHEN {ref_a} <= {alt_a} THEN {ref_a} ELSE {alt_a} END, \
             CASE WHEN {ref_a} <= {alt_a} THEN {alt_a} ELSE {ref_a} END \
         ORDER BY {pos}",
        pos = Col::Position, ref_a = Col::RefAllele, alt_a = Col::AltAllele,
        maf = Col::Maf,
        gene = Col::GeneName, region = Col::RegionType, csq = Col::Consequence,
        revel = Col::Revel,
        cage_p = Col::IsCagePromoter, cage_e = Col::IsCageEnhancer,
        ccre_p = Col::IsCcrePromoter, ccre_e = Col::IsCcreEnhancer,
    ))?;

    let batches = engine.collect(&format!(
        "SELECT {pos}, {ref_a}, {alt_a}, u_meta, \
         mac_total, n_total, {gene}, {region}, {csq}, \
         cadd_phred_raw, {revel}, \
         {cage_p}, {cage_e}, {ccre_p}, {ccre_e}, \
         {weight_select}, \
         study_segs \
         FROM _meta_variants ORDER BY {pos}",
        pos = Col::Position,
        ref_a = Col::RefAllele,
        alt_a = Col::AltAllele,
        gene = Col::GeneName,
        region = Col::RegionType,
        csq = Col::Consequence,
        revel = Col::Revel,
        cage_p = Col::IsCagePromoter,
        cage_e = Col::IsCageEnhancer,
        ccre_p = Col::IsCcrePromoter,
        ccre_e = Col::IsCcreEnhancer,
    ))?;

    let mut result = Vec::new();
    let chrom_parsed: Chromosome = chrom.parse().unwrap_or(Chromosome::Autosome(1));

    for batch in &batches {
        let n = batch.num_rows();
        let col = |i: usize| batch.column(i);
        let i32_col = |i: usize| {
            col(i).as_any().downcast_ref::<Int32Array>().ok_or_else(|| {
                CohortError::Analysis(format!("_meta_variants col {i}: expected Int32"))
            })
        };
        let i64_col = |i: usize| {
            col(i).as_any().downcast_ref::<Int64Array>().ok_or_else(|| {
                CohortError::Analysis(format!("_meta_variants col {i}: expected Int64"))
            })
        };
        let f64_col = |i: usize| {
            col(i)
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| {
                    CohortError::Analysis(format!("_meta_variants col {i}: expected Float64"))
                })
        };
        let str_col = |i: usize| {
            col(i)
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| {
                    CohortError::Analysis(format!("_meta_variants col {i}: expected Utf8"))
                })
        };
        let bool_col = |i: usize| {
            col(i)
                .as_any()
                .downcast_ref::<BooleanArray>()
                .ok_or_else(|| {
                    CohortError::Analysis(format!("_meta_variants col {i}: expected Boolean"))
                })
        };

        let pos_arr = i32_col(0)?;
        let ref_arr = str_col(1)?;
        let alt_arr = str_col(2)?;
        let u_meta_arr = f64_col(3)?;
        let mac_arr = i64_col(4)?;
        let n_obs_arr = i64_col(5)?;
        let gene_arr = str_col(6)?;
        let rt_arr = str_col(7)?;
        let csq_arr = str_col(8)?;
        let cadd_raw_arr = f64_col(9)?;
        let revel_arr = f64_col(10)?;
        let cp_arr = bool_col(11)?;
        let ce_arr = bool_col(12)?;
        let crp_arr = bool_col(13)?;
        let cre_arr = bool_col(14)?;

        // Sumstats parquet stores transformed [0,1] weights under the
        // FAVOR-native channel names; pass them through so the
        // PHRED -> probability transform is not applied a second time.
        let mut w_arrs: Vec<&Float64Array> = Vec::with_capacity(11);
        for i in 0..11 {
            w_arrs.push(f64_col(15 + i)?);
        }
        let segs_arr = str_col(26)?;

        for i in 0..n {
            let mut weights = [0.0f64; 11];
            for (ch, wa) in w_arrs.iter().enumerate() {
                weights[ch] = f64_or(wa, i, 0.0);
            }
            let cadd_phred = f64_or(cadd_raw_arr, i, 0.0);

            let study_segments = parse_study_segments(str_or(segs_arr, i, ""));

            let mac_total = i64_or(mac_arr, i, 0);
            let n_total = i64_or(n_obs_arr, i, 0);
            // MAF = MAC / (2*N) — diploid genomes have 2 alleles per sample
            let maf = if n_total > 0 {
                mac_total as f64 / (2.0 * n_total as f64)
            } else {
                0.0
            };

            result.push(MetaVariant {
                variant: AnnotatedVariant {
                    chromosome: chrom_parsed,
                    position: pos_arr.value(i) as u32,
                    ref_allele: str_or(ref_arr, i, "").into(),
                    alt_allele: str_or(alt_arr, i, "").into(),
                    maf,
                    gene_name: str_or(gene_arr, i, "").into(),
                    annotation: FunctionalAnnotation {
                        region_type: RegionType::from_str_lossy(str_or(rt_arr, i, "")),
                        consequence: Consequence::from_str_lossy(str_or(csq_arr, i, "")),
                        cadd_phred,
                        revel: f64_or(revel_arr, i, 0.0),
                        regulatory: RegulatoryFlags {
                            cage_promoter: bool_or(cp_arr, i, false),
                            cage_enhancer: bool_or(ce_arr, i, false),
                            ccre_promoter: bool_or(crp_arr, i, false),
                            ccre_enhancer: bool_or(cre_arr, i, false),
                        },
                        weights: AnnotationWeights(weights),
                    },
                },
                u_meta: f64_or(u_meta_arr, i, 0.0),
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
    segment_cache: &HashMap<(usize, i32), SegmentCov>,
) -> Option<MetaGeneResult> {
    let indices: Vec<usize> = group
        .variant_indices
        .iter()
        .filter(|&&i| i < meta_variants.len())
        .copied()
        .collect();
    if indices.len() < 2 {
        return None;
    }

    let m = indices.len();

    let mut u = Mat::zeros(m, 1);
    for (local, &gi) in indices.iter().enumerate() {
        u[(local, 0)] = meta_variants[gi].u_meta;
    }

    let keys: Vec<(u32, &str, &str)> = indices
        .iter()
        .map(|&gi| {
            (
                meta_variants[gi].variant.position,
                &*meta_variants[gi].variant.ref_allele,
                &*meta_variants[gi].variant.alt_allele,
            )
        })
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

    let mafs: Vec<f64> = indices
        .iter()
        .map(|&gi| {
            let mv = &meta_variants[gi];
            if mv.n_total > 0 {
                mv.mac_total as f64 / (2.0 * mv.n_total as f64)
            } else {
                0.0
            }
        })
        .collect();

    // Use max N across variants in this gene for MAC-based ACAT-V grouping
    let n_total: usize = indices
        .iter()
        .map(|&gi| meta_variants[gi].n_total as usize)
        .max()
        .unwrap_or(0);

    let ann_matrix: Vec<Vec<f64>> = (0..11)
        .map(|ch| {
            indices
                .iter()
                .map(|&gi| meta_variants[gi].variant.annotation.weights.0[ch])
                .collect()
        })
        .collect();

    let sr = score::run_staar_from_sumstats(&u, &cov, &ann_matrix, &mafs, n_total);

    let (burden_beta, burden_se) = unweighted_burden_estimate(&u, &cov);
    let cmac: i64 = indices.iter().map(|&gi| meta_variants[gi].mac_total).sum();

    Some(MetaGeneResult {
        gene: GeneResult {
            ensembl_id: group.name.clone(),
            gene_symbol: group.name.clone(),
            chromosome: group.chromosome,
            start: group.start,
            end: group.end,
            n_variants: m as u32,
            cumulative_mac: cmac as u32,
            staar: sr,
        },
        burden_beta,
        burden_se,
    })
}

/// Unweighted (collapsing) burden coefficient and standard error from
/// a meta-analysis score vector and its covariance:
///
///   β̂  = 1ᵀU / 1ᵀK1
///   SE = √(1 / 1ᵀK1)
///
/// Returns NaN when 1ᵀK1 ≤ 0 (degenerate gene with no usable variance).
fn unweighted_burden_estimate(u: &Mat<f64>, k: &Mat<f64>) -> (f64, f64) {
    let m = u.nrows();
    let mut one_t_u = 0.0;
    for i in 0..m {
        one_t_u += u[(i, 0)];
    }
    let mut one_t_k_one = 0.0;
    for i in 0..m {
        for j in 0..m {
            one_t_k_one += k[(i, j)];
        }
    }
    if one_t_k_one <= 0.0 || !one_t_k_one.is_finite() {
        return (f64::NAN, f64::NAN);
    }
    (one_t_u / one_t_k_one, (1.0 / one_t_k_one).sqrt())
}

/// Small arrow-null unwrap helpers. Replace the
/// `if arr.is_null(i) { default } else { arr.value(i) }` ternary that
/// used to repeat ~15 times across `merge_chromosome`'s row loop.
#[inline]
fn f64_or(arr: &Float64Array, i: usize, default: f64) -> f64 {
    if arr.is_null(i) { default } else { arr.value(i) }
}

#[inline]
fn i64_or(arr: &Int64Array, i: usize, default: i64) -> i64 {
    if arr.is_null(i) { default } else { arr.value(i) }
}

#[inline]
fn str_or<'a>(arr: &'a StringArray, i: usize, default: &'a str) -> &'a str {
    if arr.is_null(i) { default } else { arr.value(i) }
}

#[inline]
fn bool_or(arr: &BooleanArray, i: usize, default: bool) -> bool {
    if arr.is_null(i) { default } else { arr.value(i) }
}

fn parse_study_segments(s: &str) -> Vec<(usize, i32)> {
    let mut result = Vec::new();
    for part in s.split('{') {
        let part =
            part.trim_matches(|c: char| c == '[' || c == ']' || c == ',' || c == ' ' || c == '}');
        if part.is_empty() {
            continue;
        }
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

const SEGMENT_BP: u32 = 500_000;
const MAX_SEGMENT_VARIANTS: usize = 2000;

#[derive(Serialize, Deserialize)]
pub struct StudyMeta {
    pub cohort_meta_version: u32,
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
    engine: &Engine,
    cohort_id: &CohortId,
    analysis: &AnalysisVectors,
    variants: &[AnnotatedVariant],
    output_dir: &Path,
    meta: &StudyMeta,
    out: &dyn Output,
) -> Result<(), CohortError> {
    out.status("MetaSTAAR: computing summary statistics (carrier-indexed)...");

    let cohort = engine.cohort(cohort_id);
    let chrom_set: Vec<String> = variants
        .iter()
        .map(|v| v.chromosome.label())
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();

    for chrom_label in &chrom_set {
        let indices: Vec<usize> = variants
            .iter()
            .enumerate()
            .filter(|(_, v)| v.chromosome.label() == *chrom_label)
            .map(|(i, _)| i)
            .collect();
        if indices.is_empty() {
            continue;
        }

        let dir = output_dir.join(format!("chromosome={chrom_label}"));
        std::fs::create_dir_all(&dir)
            .map_err(|e| CohortError::Resource(format!("Cannot create '{}': {e}", dir.display())))?;

        let chrom: Chromosome = chrom_label
            .parse()
            .map_err(|e: String| CohortError::Input(e))?;
        let view = cohort.chromosome(&chrom)?;
        emit_chromosome_sparse(&view, analysis, variants, &indices, chrom_label, &dir, out)?;
    }

    let meta_json = serde_json::to_string_pretty(meta)
        .map_err(|e| CohortError::Resource(format!("Failed to serialize meta_staar.json: {e}")))?;
    let meta_path = output_dir.join("meta_staar.json");
    std::fs::write(&meta_path, meta_json).map_err(|e| {
        CohortError::Resource(format!("Cannot write '{}': {e}", meta_path.display()))
    })?;

    out.success(&format!("Summary statistics -> {}", output_dir.display()));
    Ok(())
}

fn emit_chromosome_sparse(
    view: &ChromosomeView<'_>,
    analysis: &AnalysisVectors,
    variants: &[AnnotatedVariant],
    chrom_indices: &[usize],
    chrom: &str,
    dir: &Path,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let variant_index = view.index()?;
    let sparse_g = view.sparse_g()?;
    let n = analysis.n_pheno;

    // Build position-based segments (same binning as v1 for output compat)
    let mut coarse: std::collections::BTreeMap<i32, Vec<usize>> = std::collections::BTreeMap::new();
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
    let mut writer = SumstatsWriter::with_capacity(total_variants, segments.len());

    for (seg_id, seg_indices) in &segments {
        let seg_id = *seg_id;
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
            writer.push_variant_row(&variants[gi], n, u[(j, 0)], k[(j, j)], seg_id);
        }
        writer.push_segment(seg_id, seg_indices, variants, &k);
    }

    writer.write_variants(dir)?;
    writer.write_segments(dir)?;

    out.status(&format!("    chr{chrom} done"));
    Ok(())
}

fn variant_schema() -> Schema {
    // `cadd_phred_raw` holds the numeric PHRED used by mask predicates
    // (MissenseDS, SpliceDS). The 11 channel columns — including one
    // named `cadd_phred` — carry the transformed [0,1] STAAR weights.
    // Two columns with distinct semantics, so distinct names.
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
        Field::new("cadd_phred_raw", DataType::Float64, false),
        Field::new(Col::Revel.as_str(), DataType::Float64, false),
        Field::new(Col::IsCagePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCageEnhancer.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcrePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcreEnhancer.as_str(), DataType::Boolean, false),
    ];
    for col in &STAAR_PHRED_CHANNELS {
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

/// Owns every Arrow builder for one chromosome's MetaSTAAR sumstats output.
/// Each column is a field; `push_variant_row` / `push_segment` keep the
/// 19 variant columns and the 6 segment columns in lock-step instead of
/// scattering the parameter list across helper boundaries.
struct SumstatsWriter {
    // Variant-row builders.
    position: Int32Builder,
    ref_allele: StringBuilder,
    alt_allele: StringBuilder,
    maf: Float64Builder,
    mac: Int32Builder,
    n_obs: Int32Builder,
    u_stat: Float64Builder,
    v_stat: Float64Builder,
    segment_id: Int32Builder,
    gene: StringBuilder,
    region: StringBuilder,
    consequence: StringBuilder,
    cadd: Float64Builder,
    revel: Float64Builder,
    cage_prom: BooleanBuilder,
    cage_enh: BooleanBuilder,
    ccre_prom: BooleanBuilder,
    ccre_enh: BooleanBuilder,
    weights: [Float64Builder; 11],

    // Segment-row builders.
    seg_id: Int32Builder,
    seg_n_variants: Int32Builder,
    seg_positions: ListBuilder<Int32Builder>,
    seg_refs: ListBuilder<StringBuilder>,
    seg_alts: ListBuilder<StringBuilder>,
    seg_cov_lower: ListBuilder<Float64Builder>,
}

impl SumstatsWriter {
    fn with_capacity(total_variants: usize, n_segments: usize) -> Self {
        Self {
            position: Int32Builder::with_capacity(total_variants),
            ref_allele: StringBuilder::with_capacity(total_variants, total_variants * 4),
            alt_allele: StringBuilder::with_capacity(total_variants, total_variants * 4),
            maf: Float64Builder::with_capacity(total_variants),
            mac: Int32Builder::with_capacity(total_variants),
            n_obs: Int32Builder::with_capacity(total_variants),
            u_stat: Float64Builder::with_capacity(total_variants),
            v_stat: Float64Builder::with_capacity(total_variants),
            segment_id: Int32Builder::with_capacity(total_variants),
            gene: StringBuilder::with_capacity(total_variants, total_variants * 8),
            region: StringBuilder::with_capacity(total_variants, total_variants * 8),
            consequence: StringBuilder::with_capacity(total_variants, total_variants * 8),
            cadd: Float64Builder::with_capacity(total_variants),
            revel: Float64Builder::with_capacity(total_variants),
            cage_prom: BooleanBuilder::with_capacity(total_variants),
            cage_enh: BooleanBuilder::with_capacity(total_variants),
            ccre_prom: BooleanBuilder::with_capacity(total_variants),
            ccre_enh: BooleanBuilder::with_capacity(total_variants),
            weights: std::array::from_fn(|_| Float64Builder::with_capacity(total_variants)),

            seg_id: Int32Builder::with_capacity(n_segments),
            seg_n_variants: Int32Builder::with_capacity(n_segments),
            seg_positions: ListBuilder::new(Int32Builder::with_capacity(total_variants)),
            seg_refs: ListBuilder::new(StringBuilder::with_capacity(
                total_variants,
                total_variants * 4,
            )),
            seg_alts: ListBuilder::new(StringBuilder::with_capacity(
                total_variants,
                total_variants * 4,
            )),
            seg_cov_lower: ListBuilder::new(Float64Builder::with_capacity(
                total_variants * (total_variants + 1) / 2,
            )),
        }
    }

    fn push_variant_row(
        &mut self,
        v: &AnnotatedVariant,
        n_samples: usize,
        u: f64,
        k_diag: f64,
        segment_id: i32,
    ) {
        self.position.append_value(v.position as i32);
        self.ref_allele.append_value(&v.ref_allele);
        self.alt_allele.append_value(&v.alt_allele);
        self.maf.append_value(v.maf);
        self.mac
            .append_value((2.0 * v.maf * n_samples as f64).round() as i32);
        self.n_obs.append_value(n_samples as i32);
        self.u_stat.append_value(u);
        self.v_stat.append_value(k_diag);
        self.segment_id.append_value(segment_id);
        self.gene.append_value(&v.gene_name);
        self.region.append_value(v.annotation.region_type.as_str());
        self.consequence
            .append_value(v.annotation.consequence.as_str());
        self.cadd.append_value(v.annotation.cadd_phred);
        self.revel.append_value(v.annotation.revel);
        self.cage_prom
            .append_value(v.annotation.regulatory.cage_promoter);
        self.cage_enh
            .append_value(v.annotation.regulatory.cage_enhancer);
        self.ccre_prom
            .append_value(v.annotation.regulatory.ccre_promoter);
        self.ccre_enh
            .append_value(v.annotation.regulatory.ccre_enhancer);
        for (i, builder) in self.weights.iter_mut().enumerate() {
            builder.append_value(v.annotation.weights.0[i]);
        }
    }

    fn push_segment(
        &mut self,
        seg_id: i32,
        seg_indices: &[usize],
        variants: &[AnnotatedVariant],
        k: &Mat<f64>,
    ) {
        let m = seg_indices.len();
        self.seg_id.append_value(seg_id);
        self.seg_n_variants.append_value(m as i32);

        let pos_builder = self.seg_positions.values();
        for &gi in seg_indices {
            pos_builder.append_value(variants[gi].position as i32);
        }
        self.seg_positions.append(true);

        let ref_builder = self.seg_refs.values();
        for &gi in seg_indices {
            ref_builder.append_value(&variants[gi].ref_allele);
        }
        self.seg_refs.append(true);

        let alt_builder = self.seg_alts.values();
        for &gi in seg_indices {
            alt_builder.append_value(&variants[gi].alt_allele);
        }
        self.seg_alts.append(true);

        let cov_builder = self.seg_cov_lower.values();
        for i in 0..m {
            for j in 0..=i {
                cov_builder.append_value(k[(i, j)]);
            }
        }
        self.seg_cov_lower.append(true);
    }

    fn write_variants(&mut self, dir: &Path) -> Result<(), CohortError> {
        let schema = Arc::new(variant_schema());
        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(self.position.finish()),
            Arc::new(self.ref_allele.finish()),
            Arc::new(self.alt_allele.finish()),
            Arc::new(self.maf.finish()),
            Arc::new(self.mac.finish()),
            Arc::new(self.n_obs.finish()),
            Arc::new(self.u_stat.finish()),
            Arc::new(self.v_stat.finish()),
            Arc::new(self.segment_id.finish()),
            Arc::new(self.gene.finish()),
            Arc::new(self.region.finish()),
            Arc::new(self.consequence.finish()),
            Arc::new(self.cadd.finish()),
            Arc::new(self.revel.finish()),
            Arc::new(self.cage_prom.finish()),
            Arc::new(self.cage_enh.finish()),
            Arc::new(self.ccre_prom.finish()),
            Arc::new(self.ccre_enh.finish()),
        ];
        for b in self.weights.iter_mut() {
            columns.push(Arc::new(b.finish()));
        }
        write_record_batch(dir, "variants.parquet", schema, columns)
    }

    fn write_segments(&mut self, dir: &Path) -> Result<(), CohortError> {
        let schema = Arc::new(segment_schema());
        let columns: Vec<ArrayRef> = vec![
            Arc::new(self.seg_id.finish()),
            Arc::new(self.seg_n_variants.finish()),
            Arc::new(self.seg_positions.finish()),
            Arc::new(self.seg_refs.finish()),
            Arc::new(self.seg_alts.finish()),
            Arc::new(self.seg_cov_lower.finish()),
        ];
        write_record_batch(dir, "segments.parquet", schema, columns)
    }
}

fn write_record_batch(
    dir: &Path,
    file_name: &str,
    schema: Arc<Schema>,
    columns: Vec<ArrayRef>,
) -> Result<(), CohortError> {
    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| CohortError::Resource(format!("Arrow batch: {e}")))?;
    let path = dir.join(file_name);
    let file = File::create(&path)
        .map_err(|e| CohortError::Resource(format!("Create {}: {e}", path.display())))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(parquet_props()))
        .map_err(|e| CohortError::Resource(format!("Parquet writer init: {e}")))?;
    writer
        .write(&batch)
        .map_err(|e| CohortError::Resource(format!("Parquet write: {e}")))?;
    writer
        .close()
        .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
    Ok(())
}

/// Parse a known-loci file into (chrom_label, position, ref, alt) tuples.
/// Format: one `chr:pos:ref:alt` per line; `#` comments and empty lines skipped.
/// Lines with only `chr:pos` are accepted (ref/alt set to "*" wildcard).
pub fn parse_known_loci_file(
    path: &Path,
) -> Result<Vec<(String, u32, String, String)>, CohortError> {
    let content = std::fs::read_to_string(path).map_err(|e| {
        CohortError::Resource(format!(
            "Cannot read known loci file '{}': {e}",
            path.display()
        ))
    })?;
    let loci: Vec<(String, u32, String, String)> = content
        .lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter_map(|l| {
            let parts: Vec<&str> = l.split(':').collect();
            if parts.len() < 2 {
                return None;
            }
            let chrom: crate::types::Chromosome = parts[0].parse().ok()?;
            let pos = parts[1].parse::<u32>().ok()?;
            let ref_a = parts.get(2).unwrap_or(&"*").to_string();
            let alt_a = parts.get(3).unwrap_or(&"*").to_string();
            Some((chrom.label().to_string(), pos, ref_a, alt_a))
        })
        .collect();
    if loci.is_empty() {
        return Err(CohortError::Input(format!(
            "Known loci file '{}' is empty or unparseable.",
            path.display()
        )));
    }
    Ok(loci)
}

/// Conditional meta-analysis: condition gene-level U/K on known loci
/// before running STAAR tests.
///
/// Homogeneous model: condition the merged (cross-study) U and K.
/// Heterogeneous model: condition per-study U and K before merging.
///
/// The conditioning step uses Schur complement:
///   U_cond = U_t - K_tc * K_cc^{-1} * U_c
///   K_cond = K_tt - K_tc * K_cc^{-1} * K_ct
///
/// where t = test (gene) variants, c = conditioning (known loci) variants.
/// Conditional meta-scoring uses the homogeneous model: condition the
/// merged (cross-study) U and K on known-loci variants via Schur complement.
/// The heterogeneous model (per-study conditioning) is rejected at config
/// time because --emit-sumstats does not yet persist per-study U vectors.
pub fn meta_score_gene_conditional(
    group: &MaskGroup,
    meta_variants: &[MetaVariant],
    studies: &[StudyHandle],
    segment_cache: &HashMap<(usize, i32), SegmentCov>,
    known_loci_indices: &[usize],
    _heterogeneous: bool,
) -> Option<MetaGeneResult> {
    let gene_indices: Vec<usize> = group
        .variant_indices
        .iter()
        .filter(|&&i| i < meta_variants.len())
        .copied()
        .collect();
    if gene_indices.len() < 2 {
        return None;
    }

    // Filter known loci to those actually present in meta_variants
    // and not overlapping with the gene's own variants.
    let gene_set: std::collections::HashSet<usize> = gene_indices.iter().copied().collect();
    let cond_indices: Vec<usize> = known_loci_indices
        .iter()
        .filter(|&&i| i < meta_variants.len() && !gene_set.contains(&i))
        .copied()
        .collect();

    // No conditioning loci in range: fall back to unconditional.
    if cond_indices.is_empty() {
        return meta_score_gene(group, meta_variants, studies, segment_cache);
    }

    let m_t = gene_indices.len();
    let m_c = cond_indices.len();

    // Build combined keys: [gene variants | conditioning variants].
    let combined: Vec<usize> = gene_indices
        .iter()
        .chain(cond_indices.iter())
        .copied()
        .collect();
    let m_all = combined.len();

    // Homogeneous: condition merged U/K.
    let mut u_all = Mat::zeros(m_all, 1);
    for (local, &gi) in combined.iter().enumerate() {
        u_all[(local, 0)] = meta_variants[gi].u_meta;
    }

    let keys: Vec<(u32, &str, &str)> = combined
        .iter()
        .map(|&gi| {
            (
                meta_variants[gi].variant.position,
                &*meta_variants[gi].variant.ref_allele,
                &*meta_variants[gi].variant.alt_allele,
            )
        })
        .collect();

    let mut cov_all = Mat::zeros(m_all, m_all);
    for study_idx in 0..studies.len() {
        let mut needed_segments: std::collections::HashSet<i32> =
            std::collections::HashSet::new();
        for &gi in &combined {
            for &(sidx, seg_id) in &meta_variants[gi].study_segments {
                if sidx == study_idx {
                    needed_segments.insert(seg_id);
                }
            }
        }
        for seg_id in needed_segments {
            if let Some(seg) = segment_cache.get(&(study_idx, seg_id)) {
                let sub = seg.extract_submatrix(&keys);
                for i in 0..m_all {
                    for j in 0..m_all {
                        cov_all[(i, j)] += sub[(i, j)];
                    }
                }
            }
        }
    }

    let (u_t, k_tt, u_c, k_cc, k_tc) = partition_u_cov(&u_all, &cov_all, m_t, m_c);
    let (u_cond, k_cond) = schur_condition(&u_t, &k_tt, &u_c, &k_cc, &k_tc);

    finish_conditional(&gene_indices, meta_variants, &u_cond, &k_cond, group)
}

/// Partition combined U and K into test (t) and conditioning (c) blocks.
#[allow(clippy::type_complexity)]
fn partition_u_cov(
    u: &Mat<f64>,
    k: &Mat<f64>,
    m_t: usize,
    m_c: usize,
) -> (Mat<f64>, Mat<f64>, Mat<f64>, Mat<f64>, Mat<f64>) {
    let u_t = u.subrows(0, m_t).to_owned();
    let u_c = u.subrows(m_t, m_c).to_owned();
    let k_tt = k.submatrix(0, 0, m_t, m_t).to_owned();
    let k_cc = k.submatrix(m_t, m_t, m_c, m_c).to_owned();
    let k_tc = k.submatrix(0, m_t, m_t, m_c).to_owned();
    (u_t, k_tt, u_c, k_cc, k_tc)
}

/// Schur complement conditioning:
///   U_cond = U_t - K_tc * K_cc^{-1} * U_c
///   K_cond = K_tt - K_tc * K_cc^{-1} * K_ct
fn schur_condition(
    u_t: &Mat<f64>,
    k_tt: &Mat<f64>,
    u_c: &Mat<f64>,
    k_cc: &Mat<f64>,
    k_tc: &Mat<f64>,
) -> (Mat<f64>, Mat<f64>) {
    use faer::linalg::solvers::Solve;

    let m_t = u_t.nrows();
    let m_c = u_c.nrows();

    if m_c == 0 {
        return (u_t.to_owned(), k_tt.to_owned());
    }

    // Regularize K_cc for numerical stability.
    let mut k_cc_reg = k_cc.to_owned();
    for i in 0..m_c {
        k_cc_reg[(i, i)] += 1e-8;
    }

    // K_cc^{-1} via column-pivoted QR (robust for small m_c).
    let k_cc_inv = {
        let eye = Mat::identity(m_c, m_c);
        k_cc_reg.col_piv_qr().solve(&eye)
    };

    // K_tc * K_cc^{-1}
    let k_tc_kinv = k_tc * &k_cc_inv;

    // U_cond = U_t - K_tc * K_cc^{-1} * U_c
    let mut u_cond = u_t.to_owned();
    let correction_u = &k_tc_kinv * u_c;
    for i in 0..m_t {
        u_cond[(i, 0)] -= correction_u[(i, 0)];
    }

    // K_cond = K_tt - K_tc * K_cc^{-1} * K_ct
    let k_ct = k_tc.transpose().to_owned();
    let correction_k = &k_tc_kinv * &k_ct;
    let mut k_cond = k_tt.to_owned();
    for i in 0..m_t {
        for j in 0..m_t {
            k_cond[(i, j)] -= correction_k[(i, j)];
        }
    }

    (u_cond, k_cond)
}

/// Finish conditional scoring: compute MAFs, annotation weights, run STAAR.
fn finish_conditional(
    gene_indices: &[usize],
    meta_variants: &[MetaVariant],
    u_cond: &Mat<f64>,
    k_cond: &Mat<f64>,
    group: &MaskGroup,
) -> Option<MetaGeneResult> {
    let m_t = gene_indices.len();

    let mafs: Vec<f64> = gene_indices
        .iter()
        .map(|&gi| {
            let mv = &meta_variants[gi];
            if mv.n_total > 0 {
                mv.mac_total as f64 / (2.0 * mv.n_total as f64)
            } else {
                0.0
            }
        })
        .collect();

    let n_total: usize = gene_indices
        .iter()
        .map(|&gi| meta_variants[gi].n_total as usize)
        .max()
        .unwrap_or(0);

    let ann_matrix: Vec<Vec<f64>> = (0..11)
        .map(|ch| {
            gene_indices
                .iter()
                .map(|&gi| meta_variants[gi].variant.annotation.weights.0[ch])
                .collect()
        })
        .collect();

    let sr = score::run_staar_from_sumstats(u_cond, k_cond, &ann_matrix, &mafs, n_total);
    let (burden_beta, burden_se) = unweighted_burden_estimate(u_cond, k_cond);
    let cmac: i64 = gene_indices
        .iter()
        .map(|&gi| meta_variants[gi].mac_total)
        .sum();

    Some(MetaGeneResult {
        gene: GeneResult {
            ensembl_id: group.name.clone(),
            gene_symbol: group.name.clone(),
            chromosome: group.chromosome,
            start: group.start,
            end: group.end,
            n_variants: m_t as u32,
            cumulative_mac: cmac as u32,
            staar: sr,
        },
        burden_beta,
        burden_se,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::staar::masks::MaskGroup;
    use crate::staar::score;
    use crate::types::{
        AnnotationWeights, Consequence, FunctionalAnnotation, RegionType, RegulatoryFlags,
    };

    /// Lower-triangular packing must match what `extract_submatrix` reads:
    /// `cov_lower[i*(i+1)/2 + j]` for `i ≥ j`.
    fn pack_lower(k: &Mat<f64>) -> Vec<f64> {
        let m = k.nrows();
        let mut packed = Vec::with_capacity(m * (m + 1) / 2);
        for i in 0..m {
            for j in 0..=i {
                packed.push(k[(i, j)]);
            }
        }
        packed
    }

    fn make_variant(pos: u32, ref_b: &str, alt_b: &str, maf: f64) -> AnnotatedVariant {
        AnnotatedVariant {
            chromosome: Chromosome::Autosome(22),
            position: pos,
            ref_allele: ref_b.into(),
            alt_allele: alt_b.into(),
            maf,
            gene_name: "GENE".into(),
            annotation: FunctionalAnnotation {
                region_type: RegionType::Exonic,
                consequence: Consequence::NonsynonymousSNV,
                cadd_phred: 20.0,
                revel: 0.5,
                regulatory: RegulatoryFlags::default(),
                weights: AnnotationWeights([1.0; 11]),
            },
        }
    }

    fn segment_cov(positions: &[u32], refs: &[&str], alts: &[&str], k: &Mat<f64>) -> SegmentCov {
        let mut pos_index: HashMap<u32, Vec<usize>> = HashMap::new();
        for (i, &p) in positions.iter().enumerate() {
            pos_index.entry(p).or_default().push(i);
        }
        SegmentCov {
            refs: refs.iter().map(|s| s.to_string()).collect(),
            alts: alts.iter().map(|s| s.to_string()).collect(),
            cov_lower: pack_lower(k),
            pos_index,
        }
    }

    fn dummy_study(n: usize) -> StudyHandle {
        StudyHandle {
            path: std::path::PathBuf::new(),
            meta: StudyMeta {
                cohort_meta_version: 1,
                trait_type: "Continuous".into(),
                trait_name: "TRAIT".into(),
                n_samples: n,
                sigma2: 1.0,
                maf_cutoff: 0.01,
                covariates: vec![],
                segment_size: SEGMENT_BP,
            },
        }
    }

    type Fixture = (Mat<f64>, Mat<f64>, Vec<f64>, Vec<Vec<f64>>, Vec<u32>);
    fn fixture(m: usize, n: usize) -> Fixture {
        let u = Mat::<f64>::from_fn(m, 1, |i, _| 0.3 * (i as f64 + 1.0));
        // Diagonal-dominant ⇒ PSD by construction.
        let mut k = Mat::<f64>::zeros(m, m);
        for i in 0..m {
            k[(i, i)] = 1.0 + i as f64 * 0.1;
            for j in 0..i {
                let v = 0.05 * (i + j) as f64;
                k[(i, j)] = v;
                k[(j, i)] = v;
            }
        }
        // 2*maf*n integer ⇒ mac round-trips through (mac/(2n)) without drift.
        let mafs: Vec<f64> = (0..m)
            .map(|i| (i as f64 + 1.0) * 10.0 / (2.0 * n as f64))
            .collect();
        let ann: Vec<Vec<f64>> = (0..11).map(|_| vec![1.0; m]).collect();
        let positions: Vec<u32> = (0..m as u32).map(|i| 1000 * (i + 1)).collect();
        (u, k, mafs, ann, positions)
    }

    /// One study through `meta_score_gene` must reproduce direct
    /// `run_staar_from_sumstats`: SUM-of-one is identity, sub-matrix
    /// extraction is identity.
    #[test]
    fn k1_meta_matches_direct_sumstats() {
        let m = 4;
        let n = 10_000usize;
        let (u, k, mafs, ann, positions) = fixture(m, n);
        let direct = score::run_staar_from_sumstats(&u, &k, &ann, &mafs, n);
        let (expected_beta, expected_se) = unweighted_burden_estimate(&u, &k);

        let refs: Vec<&str> = vec!["A"; m];
        let alts: Vec<&str> = vec!["C"; m];
        let variants: Vec<MetaVariant> = (0..m)
            .map(|i| MetaVariant {
                variant: make_variant(positions[i], refs[i], alts[i], mafs[i]),
                u_meta: u[(i, 0)],
                mac_total: (2.0 * mafs[i] * n as f64).round() as i64,
                n_total: n as i64,
                study_segments: vec![(0, 0)],
            })
            .collect();

        let group = MaskGroup {
            name: "GENE".into(),
            chromosome: Chromosome::Autosome(22),
            start: positions[0],
            end: *positions.last().unwrap(),
            variant_indices: (0..m).collect(),
        };

        let mut cache = HashMap::new();
        cache.insert((0, 0), segment_cov(&positions, &refs, &alts, &k));
        let studies = vec![dummy_study(n)];

        let result = meta_score_gene(&group, &variants, &studies, &cache).expect("score");
        assert!(
            (result.gene.staar.staar_o - direct.staar_o).abs() < 1e-12,
            "K=1 meta {} vs direct {}",
            result.gene.staar.staar_o,
            direct.staar_o
        );
        assert_eq!(result.gene.n_variants, m as u32);
        assert!((result.burden_beta - expected_beta).abs() < 1e-12);
        assert!((result.burden_se - expected_se).abs() < 1e-12);
    }

    #[test]
    fn canonical_orientation_is_lex_min() {
        assert_eq!(canonical_alleles("A", "C"), ("A", "C"));
        assert_eq!(canonical_alleles("C", "A"), ("A", "C"));
        // Indels: lex order extends naturally; "C" < "CT".
        assert_eq!(canonical_alleles("CT", "C"), ("C", "CT"));
        assert_eq!(canonical_alleles("C", "CT"), ("C", "CT"));
        // Equal alleles never occur in real VCF; the function is total anyway.
        assert_eq!(canonical_alleles("A", "A"), ("A", "A"));
    }

    /// A segment that stored a variant flipped (`alt`,`ref`) must still
    /// resolve when the caller asks in canonical (`ref`,`alt`) order.
    #[test]
    fn extract_submatrix_is_orientation_blind() {
        let positions = vec![100u32];
        let k = Mat::<f64>::from_fn(1, 1, |_, _| 4.2);
        let local_flipped = segment_cov(&positions, &["C"], &["A"], &k);
        let pulled = local_flipped.extract_submatrix(&[(100, "A", "C")]);
        assert!((pulled[(0, 0)] - 4.2).abs() < 1e-12);
    }

    /// Off-diagonal `K[i,j] = G_i·G_j` flips sign when exactly one of the
    /// two variants was stored in flipped orientation. Diagonal is invariant.
    #[test]
    fn extract_submatrix_signs_off_diagonal_on_partial_flip() {
        // Two variants: variant A at pos 100 (canonical = ("A","C")), variant
        // B at pos 200 (canonical = ("G","T")). Store the segment with B in
        // FLIPPED orientation, so the local cov[A,B] is signed against
        // (-G_B), not against +G_B.
        let positions = vec![100u32, 200u32];
        let canonical_k = Mat::<f64>::from_fn(2, 2, |i, j| match (i, j) {
            (0, 0) => 5.0,
            (1, 1) => 3.0,
            _ => 1.5, // off-diagonal
        });
        // Build the LOCAL k that would have been emitted if B were stored
        // flipped: row 1 / col 1 negated except for the diagonal.
        let local_k = Mat::<f64>::from_fn(2, 2, |i, j| {
            let s = if i == 1 || j == 1 { -1.0 } else { 1.0 };
            // diagonal is sign-squared = +1
            let s = if i == j { 1.0 } else { s };
            canonical_k[(i, j)] * s
        });
        assert!((local_k[(0, 1)] - (-1.5)).abs() < 1e-12);

        let seg = segment_cov(&positions, &["A", "T"], &["C", "G"], &local_k);

        // Caller asks in canonical orientation.
        let pulled = seg.extract_submatrix(&[(100, "A", "C"), (200, "G", "T")]);
        assert!((pulled[(0, 0)] - 5.0).abs() < 1e-12);
        assert!((pulled[(1, 1)] - 3.0).abs() < 1e-12);
        assert!(
            (pulled[(0, 1)] - 1.5).abs() < 1e-12,
            "off-diagonal must be sign-corrected: got {}",
            pulled[(0, 1)]
        );
        assert!((pulled[(1, 0)] - 1.5).abs() < 1e-12);
    }

    /// Two studies whose summary statistics algebraically sum to one
    /// pseudo-combined study must yield the same staar_o as direct
    /// `run_staar_from_sumstats` on the combined U/K. This pins the
    /// "split-and-merge equals direct" invariant for MetaSTAAR.
    #[test]
    fn k2_split_recovers_combined_baseline() {
        let m = 4;
        let n_a = 5_000usize;
        let n_b = 5_000usize;
        let n = n_a + n_b;
        let (u, k, mafs, ann, positions) = fixture(m, n);
        let direct = score::run_staar_from_sumstats(&u, &k, &ann, &mafs, n);

        let half_k = Mat::<f64>::from_fn(m, m, |i, j| k[(i, j)] * 0.5);
        let half_u = Mat::<f64>::from_fn(m, 1, |i, _| u[(i, 0)] * 0.5);

        let refs: Vec<&str> = vec!["A"; m];
        let alts: Vec<&str> = vec!["C"; m];
        let variants: Vec<MetaVariant> = (0..m)
            .map(|i| MetaVariant {
                variant: make_variant(positions[i], refs[i], alts[i], mafs[i]),
                u_meta: half_u[(i, 0)] + half_u[(i, 0)],
                mac_total: (2.0 * mafs[i] * n as f64).round() as i64,
                n_total: n as i64,
                study_segments: vec![(0, 10), (1, 20)],
            })
            .collect();

        let group = MaskGroup {
            name: "GENE".into(),
            chromosome: Chromosome::Autosome(22),
            start: positions[0],
            end: *positions.last().unwrap(),
            variant_indices: (0..m).collect(),
        };

        let mut cache = HashMap::new();
        cache.insert((0, 10), segment_cov(&positions, &refs, &alts, &half_k));
        cache.insert((1, 20), segment_cov(&positions, &refs, &alts, &half_k));
        let studies = vec![dummy_study(n_a), dummy_study(n_b)];

        let result = meta_score_gene(&group, &variants, &studies, &cache).expect("score");
        assert!(
            (result.gene.staar.staar_o - direct.staar_o).abs() < 1e-12,
            "K=2 split meta {} vs direct {}",
            result.gene.staar.staar_o,
            direct.staar_o
        );
    }
}
