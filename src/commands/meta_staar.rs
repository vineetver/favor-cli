use std::fs::File;
use std::path::PathBuf;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float64Builder, StringBuilder, UInt32Builder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;
use rayon::prelude::*;
use serde_json::json;

use crate::column::{Col, STAAR_WEIGHTS};
use crate::commands::MetaStaarConfig;
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::staar::masks::MaskGroup;
use crate::staar::{self, GeneResult, MaskCategory, MaskType};
use crate::types::Chromosome;

pub fn build_config(
    studies: Vec<PathBuf>,
    masks: Vec<String>,
    maf_cutoff: f64,
    window_size: u32,
    output_path: Option<PathBuf>,
) -> Result<MetaStaarConfig, CohortError> {
    if studies.is_empty() {
        return Err(CohortError::Input(
            "--studies is required. Provide comma-separated paths to MetaSTAAR summary stat directories.".into(),
        ));
    }

    let mask_categories = crate::commands::parse_mask_categories(&masks)?;
    let output_dir = output_path.unwrap_or_else(|| PathBuf::from("meta_staar_results"));

    Ok(MetaStaarConfig {
        study_dirs: studies,
        mask_categories,
        maf_cutoff,
        window_size,
        output_dir,
    })
}

pub fn run_meta_staar(
    engine: &Engine,
    config: &MetaStaarConfig,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    let studies = staar::meta::load_studies(&config.study_dirs)?;

    if dry_run {
        return emit_dry_run(engine, &studies, config, out);
    }

    let k = studies.len();
    let total_n: usize = studies.iter().map(|s| s.meta.n_samples).sum();
    out.status(&format!("MetaSTAAR: {k} studies, {total_n} total samples"));
    for (i, s) in studies.iter().enumerate() {
        out.status(&format!(
            "  [{i}] {} — N={}, trait={}",
            s.path.display(),
            s.meta.n_samples,
            s.meta.trait_name
        ));
    }
    let resources = engine.resources();
    out.status(&format!(
        "MetaSTAAR: {} memory, {} threads",
        resources.memory_human(),
        resources.threads
    ));

    let df = engine.df();

    std::fs::create_dir_all(&config.output_dir).map_err(|e| {
        CohortError::Resource(format!(
            "Cannot create '{}': {e}",
            config.output_dir.display()
        ))
    })?;

    let results = run_all_chromosomes(&studies, config, df, out)?;
    write_meta_results(&results, &config.output_dir, out)?;
    generate_summary(&studies, &results, config, out);

    out.success(&format!(
        "MetaSTAAR complete -> {}",
        config.output_dir.display()
    ));
    out.result_json(&json!({ "command": "meta-staar", "n_studies": k, "total_samples": total_n }));
    Ok(())
}

fn emit_dry_run(
    engine: &Engine,
    studies: &[staar::meta::StudyHandle],
    config: &MetaStaarConfig,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let k = studies.len();
    let total_n: usize = studies.iter().map(|s| s.meta.n_samples).sum();
    let threads = engine.resources().threads.max(1);
    let runtime_seconds = meta_staar_runtime_seconds(k, total_n, threads);

    let plan = crate::commands::DryRunPlan {
        command: "meta-staar".into(),
        inputs: json!({
            "n_studies": k,
            "total_samples": total_n,
            "trait_type": studies[0].meta.trait_type,
        }),
        memory: crate::commands::MemoryEstimate::default_estimate(),
        runtime: Some(crate::commands::RuntimeEstimate::from_seconds(runtime_seconds)),
        output_path: config.output_dir.to_string_lossy().into(),
    };
    crate::commands::emit(&plan, out);
    Ok(())
}

/// Coarse wall-clock budget for `cohort meta-staar`. Dominated by the
/// per-chromosome SQL merge plus per-segment covariance accumulation,
/// both linear in `(n_studies × total_samples)`.
fn meta_staar_runtime_seconds(n_studies: usize, total_samples: usize, threads: usize) -> u64 {
    const NS_PER_STUDY_SAMPLE: u64 = 5_000;
    const SQL_BASE_S: u64 = 20;

    let merge_ns = (n_studies as u64).saturating_mul(total_samples as u64) * NS_PER_STUDY_SAMPLE;
    let merge_s = merge_ns / 1_000_000_000 / threads as u64;
    SQL_BASE_S.saturating_add(merge_s).max(30)
}

fn run_all_chromosomes(
    studies: &[staar::meta::StudyHandle],
    config: &MetaStaarConfig,
    engine: &DfEngine,
    out: &dyn Output,
) -> Result<Vec<(MaskType, Vec<GeneResult>)>, CohortError> {
    let k = studies.len();
    let chromosomes = discover_chromosomes(&studies[0].path);
    let mut all_results: Vec<(MaskType, Vec<GeneResult>)> = Vec::new();

    for chrom in &chromosomes {
        out.status(&format!(
            "  chr{chrom}: merging variants across {k} studies..."
        ));

        let meta_variants =
            staar::meta::merge_chromosome(engine, studies, chrom, config.maf_cutoff)?;
        if meta_variants.is_empty() {
            out.status(&format!("    chr{chrom}: no variants"));
            continue;
        }
        out.status(&format!("    {} merged variants", meta_variants.len()));

        let annotated: Vec<crate::types::AnnotatedVariant> =
            meta_variants.iter().map(|mv| mv.variant.clone()).collect();
        let chrom_indices: Vec<usize> = (0..annotated.len()).collect();

        let chrom_masks = build_chrom_masks(&annotated, &chrom_indices, chrom, config);

        let segment_cache = load_segment_cache(studies, chrom, out);

        for (mask_type, groups) in &chrom_masks {
            let results: Vec<GeneResult> = groups
                .par_iter()
                .filter_map(|group| {
                    staar::meta::meta_score_gene(
                        group,
                        &meta_variants,
                        studies,
                        &segment_cache,
                    )
                })
                .collect();

            if !results.is_empty() {
                out.status(&format!(
                    "    {}: {} groups, {} with results",
                    mask_type.file_stem(),
                    groups.len(),
                    results.len()
                ));
                if let Some(existing) = all_results.iter_mut().find(|(mt, _)| mt == mask_type) {
                    existing.1.extend(results);
                } else {
                    all_results.push((mask_type.clone(), results));
                }
            }
        }
    }

    Ok(all_results)
}

fn build_chrom_masks(
    annotated: &[crate::types::AnnotatedVariant],
    chrom_indices: &[usize],
    chrom: &str,
    config: &MetaStaarConfig,
) -> Vec<(MaskType, Vec<MaskGroup>)> {
    let chrom_parsed: Chromosome = chrom.parse().unwrap_or(Chromosome::Autosome(1));
    let mut masks = Vec::new();
    for cat in &config.mask_categories {
        match cat {
            MaskCategory::Coding => {
                masks.extend(staar::masks::build_coding_masks(annotated, 2));
            }
            MaskCategory::Noncoding => {
                masks.extend(staar::masks::build_noncoding_masks(annotated, 2));
            }
            MaskCategory::SlidingWindow => {
                let windows = staar::masks::build_sliding_windows(
                    annotated,
                    chrom_indices,
                    chrom_parsed,
                    config.window_size,
                    config.window_size / 2,
                );
                if !windows.is_empty() {
                    masks.push((MaskType::SlidingWindow, windows));
                }
            }
            _ => {}
        }
    }
    masks
}

fn load_segment_cache(
    studies: &[staar::meta::StudyHandle],
    chrom: &str,
    out: &dyn Output,
) -> std::collections::HashMap<(usize, i32), staar::meta::SegmentCov> {
    let mut cache = std::collections::HashMap::new();
    for (study_idx, study) in studies.iter().enumerate() {
        match staar::meta::load_all_segments(study, chrom) {
            Ok(segs) => {
                for (seg_id, seg_cov) in segs {
                    cache.insert((study_idx, seg_id), seg_cov);
                }
            }
            Err(e) => out.warn(&format!(
                "    Failed to load segments for study {study_idx}: {e}"
            )),
        }
    }
    cache
}

fn discover_chromosomes(study_dir: &std::path::Path) -> Vec<String> {
    let mut chroms = Vec::new();
    if let Ok(entries) = std::fs::read_dir(study_dir) {
        for entry in entries.flatten() {
            let name = entry.file_name().to_string_lossy().to_string();
            if let Some(chrom) = name.strip_prefix("chromosome=") {
                if entry.path().join("variants.parquet").exists() {
                    chroms.push(chrom.to_string());
                }
            }
        }
    }
    chroms.sort_by(|a, b| match (a.parse::<u32>(), b.parse::<u32>()) {
        (Ok(a), Ok(b)) => a.cmp(&b),
        _ => a.cmp(b),
    });
    chroms
}

fn generate_summary(
    studies: &[staar::meta::StudyHandle],
    results: &[(MaskType, Vec<GeneResult>)],
    config: &MetaStaarConfig,
    out: &dyn Output,
) {
    let k = studies.len();
    let total_n: usize = studies.iter().map(|s| s.meta.n_samples).sum();
    let n_rare: i64 = results.iter().map(|(_, r)| r.len() as i64).sum();
    let trait_names = vec![studies[0].meta.trait_name.clone()];
    let title = format!("MetaSTAAR Meta-Analysis ({k} studies, N={total_n})");

    match staar::output::generate_report(
        results,
        &trait_names,
        total_n,
        n_rare,
        &config.output_dir,
        &title,
    ) {
        Ok(()) => out.success(&format!(
            "  summary.html -> {}",
            config.output_dir.join("summary.html").display()
        )),
        Err(e) => out.warn(&format!("  Summary report failed: {e}")),
    }
}

fn write_meta_results(
    all_mask_results: &[(MaskType, Vec<GeneResult>)],
    output_dir: &std::path::Path,
    out: &dyn Output,
) -> Result<(), CohortError> {
    out.status("Writing MetaSTAAR results...");
    let channels: Vec<&str> = STAAR_WEIGHTS
        .iter()
        .map(|c| c.weight_display_name().expect("STAAR_WEIGHTS entries have display names"))
        .collect();

    for (mask_type, results) in all_mask_results {
        if results.is_empty() {
            continue;
        }
        write_mask_results(mask_type, results, &channels, output_dir, out)?;
    }
    Ok(())
}

fn write_mask_results(
    mask_type: &MaskType,
    results: &[GeneResult],
    channels: &[&str],
    output_dir: &std::path::Path,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let out_path = output_dir.join(format!("{}.parquet", mask_type.file_stem()));
    let n_channels = channels.len();

    let mut sorted: Vec<&GeneResult> = results.iter().collect();
    sorted.sort_by(|a, b| {
        a.staar
            .staar_o
            .partial_cmp(&b.staar.staar_o)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let (schema, batch) = build_result_batch(&sorted, channels, n_channels)?;

    let file = File::create(&out_path)
        .map_err(|e| CohortError::Resource(format!("Create {}: {e}", out_path.display())))?;
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| CohortError::Resource(format!("Parquet writer: {e}")))?;
    writer
        .write(&batch)
        .map_err(|e| CohortError::Resource(format!("Parquet write: {e}")))?;
    writer
        .close()
        .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;

    let n_sig = results.iter().filter(|r| r.staar.staar_o < 2.5e-6).count();
    out.success(&format!(
        "  {} -> {} genes, {} significant",
        mask_type.file_stem(),
        results.len(),
        n_sig
    ));

    Ok(())
}

fn build_result_batch(
    sorted: &[&GeneResult],
    channels: &[&str],
    n_channels: usize,
) -> Result<(Arc<Schema>, RecordBatch), CohortError> {
    let nr = sorted.len();
    let mut b_ensembl = StringBuilder::with_capacity(nr, nr * 16);
    let mut b_symbol = StringBuilder::with_capacity(nr, nr * 12);
    let mut b_chrom = StringBuilder::with_capacity(nr, nr * 2);
    let mut b_start = UInt32Builder::with_capacity(nr);
    let mut b_end = UInt32Builder::with_capacity(nr);
    let mut b_nvariants = UInt32Builder::with_capacity(nr);
    let mut b_cmac = UInt32Builder::with_capacity(nr);
    let mut b_burden_beta = Float64Builder::with_capacity(nr);
    let mut b_burden_se = Float64Builder::with_capacity(nr);

    let n_pval_cols = 6 + 6 * n_channels + 6 + 2;
    let mut pval_builders: Vec<Float64Builder> = (0..n_pval_cols)
        .map(|_| Float64Builder::with_capacity(nr))
        .collect();

    for r in sorted {
        let s = &r.staar;
        b_ensembl.append_value(&r.ensembl_id);
        b_symbol.append_value(&r.gene_symbol);
        b_chrom.append_value(r.chromosome.label());
        b_start.append_value(r.start);
        b_end.append_value(r.end);
        b_nvariants.append_value(r.n_variants);
        b_cmac.append_value(r.cumulative_mac);
        b_burden_beta.append_value(r.burden_beta);
        b_burden_se.append_value(r.burden_se);

        let mut pi = 0;
        for p in [
            s.burden_1_25,
            s.burden_1_1,
            s.skat_1_25,
            s.skat_1_1,
            s.acat_v_1_25,
            s.acat_v_1_1,
        ] {
            pval_builders[pi].append_value(p);
            pi += 1;
        }
        for ann_p in &s.per_annotation {
            for &v in ann_p {
                pval_builders[pi].append_value(v);
                pi += 1;
            }
        }
        while pi < 6 + 6 * n_channels {
            pval_builders[pi].append_value(f64::NAN);
            pi += 1;
        }
        for p in [
            s.staar_b_1_25,
            s.staar_b_1_1,
            s.staar_s_1_25,
            s.staar_s_1_1,
            s.staar_a_1_25,
            s.staar_a_1_1,
            s.acat_o,
            s.staar_o,
        ] {
            pval_builders[pi].append_value(p);
            pi += 1;
        }
    }

    let test_names = [
        "Burden(1,25)",
        "Burden(1,1)",
        "SKAT(1,25)",
        "SKAT(1,1)",
        "ACAT-V(1,25)",
        "ACAT-V(1,1)",
    ];
    let mut fields = vec![
        Field::new("ensembl_id", DataType::Utf8, false),
        Field::new("gene_symbol", DataType::Utf8, false),
        Field::new(Col::Chromosome.as_str(), DataType::Utf8, false),
        Field::new("start", DataType::UInt32, false),
        Field::new("end", DataType::UInt32, false),
        Field::new("n_variants", DataType::UInt32, false),
        Field::new("cMAC", DataType::UInt32, false),
        Field::new("burden_beta", DataType::Float64, true),
        Field::new("burden_se", DataType::Float64, true),
    ];
    for test in &test_names {
        fields.push(Field::new(*test, DataType::Float64, true));
    }
    for ch in channels {
        for test in &test_names {
            fields.push(Field::new(format!("{test}-{ch}"), DataType::Float64, true));
        }
    }
    for name in [
        "STAAR-B(1,25)",
        "STAAR-B(1,1)",
        "STAAR-S(1,25)",
        "STAAR-S(1,1)",
        "STAAR-A(1,25)",
        "STAAR-A(1,1)",
        "ACAT-O",
        "STAAR-O",
    ] {
        fields.push(Field::new(name, DataType::Float64, true));
    }
    let schema = Arc::new(Schema::new(fields));

    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(b_ensembl.finish()),
        Arc::new(b_symbol.finish()),
        Arc::new(b_chrom.finish()),
        Arc::new(b_start.finish()),
        Arc::new(b_end.finish()),
        Arc::new(b_nvariants.finish()),
        Arc::new(b_cmac.finish()),
        Arc::new(b_burden_beta.finish()),
        Arc::new(b_burden_se.finish()),
    ];
    for b in &mut pval_builders {
        columns.push(Arc::new(b.finish()));
    }

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| CohortError::Resource(format!("Arrow batch: {e}")))?;

    Ok((schema, batch))
}
