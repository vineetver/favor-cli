//! `favor staar` command — thin orchestrator.
//!
//! Validates CLI arguments, constructs a `StaarConfig`, and delegates to
//! `StaarPipeline` for the actual work.

use std::path::PathBuf;

use serde_json::json;

use crate::commands;
use crate::data::VariantSet;
use crate::error::FavorError;
use crate::output::Output;
use crate::staar::store;
use crate::staar::pipeline::{StaarConfig, StaarPipeline};
use crate::staar::{self, MaskCategory};

const GB: u64 = 1024 * 1024 * 1024;

#[allow(clippy::too_many_arguments)]
pub fn run(
    genotypes: PathBuf,
    phenotype: PathBuf,
    trait_names: Vec<String>,
    covariates: Vec<String>,
    annotations: Option<PathBuf>,
    masks: Vec<String>,
    maf_cutoff: f64,
    window_size: u32,
    individual: bool,
    spa: bool,
    scang_lmin: usize,
    scang_lmax: usize,
    scang_step: usize,
    known_loci: Option<PathBuf>,
    emit_sumstats: bool,
    rebuild_store: bool,
    column_map: Vec<String>,
    output_path: Option<PathBuf>,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), FavorError> {
    let config = validate_and_parse(
        genotypes, phenotype, trait_names, covariates, annotations, masks,
        maf_cutoff, window_size, individual, spa,
        scang_lmin, scang_lmax, scang_step, known_loci, emit_sumstats,
        rebuild_store, column_map, output_path,
    )?;

    if dry_run {
        return emit_dry_run(&config, out);
    }

    let pipeline = StaarPipeline::new(config, out)?;
    pipeline.run()
}

#[allow(clippy::too_many_arguments)]
fn validate_and_parse(
    genotypes: PathBuf,
    phenotype: PathBuf,
    trait_names: Vec<String>,
    covariates: Vec<String>,
    annotations: Option<PathBuf>,
    masks: Vec<String>,
    maf_cutoff: f64,
    window_size: u32,
    individual: bool,
    spa: bool,
    scang_lmin: usize,
    scang_lmax: usize,
    scang_step: usize,
    known_loci: Option<PathBuf>,
    emit_sumstats: bool,
    rebuild_store: bool,
    column_map_raw: Vec<String>,
    output_path: Option<PathBuf>,
) -> Result<StaarConfig, FavorError> {
    if !genotypes.exists() {
        return Err(FavorError::Input(format!(
            "Genotype VCF not found: '{}'. Check the path to your multi-sample VCF.",
            genotypes.display()
        )));
    }
    if !phenotype.exists() {
        return Err(FavorError::Input(format!(
            "Phenotype file not found: '{}'. Provide a tab-delimited file with sample IDs, traits, and covariates.",
            phenotype.display()
        )));
    }
    if maf_cutoff <= 0.0 || maf_cutoff >= 0.5 {
        return Err(FavorError::Input(format!(
            "MAF cutoff must be in (0, 0.5), got {maf_cutoff}"
        )));
    }
    let annotations_path = annotations.ok_or_else(|| {
        FavorError::Input("STAAR requires --annotations <path> from `favor annotate`.".into())
    })?;
    if !annotations_path.exists() {
        return Err(FavorError::Input(format!(
            "Annotations not found: {}. Run `favor annotate` first.",
            annotations_path.display()
        )));
    }
    let ann_vs = VariantSet::open(&annotations_path).map_err(|e| {
        FavorError::Input(format!(
            "'{}' is not a valid variant set ({}). Run `favor annotate` to produce one.",
            annotations_path.display(), e,
        ))
    })?;
    ann_vs.require_annotated()?;

    let mask_categories: Vec<MaskCategory> = masks
        .iter()
        .map(|s| {
            s.parse::<MaskCategory>().map_err(|_| {
                FavorError::Input(format!(
                    "Unknown mask '{s}'. Available: coding, noncoding, sliding-window, scang, custom"
                ))
            })
        })
        .collect::<Result<_, _>>()?;

    let output_dir = output_path.unwrap_or_else(|| {
        let stem = genotypes
            .file_stem()
            .unwrap_or_default()
            .to_string_lossy();
        genotypes.with_file_name(format!("{stem}.staar"))
    });

    if let Some(ref loci) = known_loci {
        if !loci.exists() {
            return Err(FavorError::Input(format!(
                "Known loci file not found: {}",
                loci.display()
            )));
        }
    }

    if trait_names.is_empty() {
        return Err(FavorError::Input(
            "At least one --trait-name is required.".into(),
        ));
    }

    let mut column_map = std::collections::HashMap::new();
    for entry in &column_map_raw {
        if let Some((k, v)) = entry.split_once('=') {
            column_map.insert(k.trim().to_string(), v.trim().to_string());
        } else {
            return Err(FavorError::Input(format!(
                "Invalid --column-map entry '{entry}'. Expected key=value."
            )));
        }
    }

    let store_dir = output_dir.join("store");

    Ok(StaarConfig {
        genotypes,
        phenotype,
        annotations: annotations_path,
        trait_names,
        covariates,
        mask_categories,
        maf_cutoff,
        window_size,
        individual,
        spa,
        scang_params: staar::masks::ScangParams {
            lmin: scang_lmin,
            lmax: scang_lmax,
            step: scang_step,
        },
        known_loci,
        emit_sumstats,
        rebuild_store,
        column_map,
        output_dir,
        store_dir,
    })
}

fn emit_dry_run(config: &StaarConfig, out: &dyn Output) -> Result<(), FavorError> {
    let sample_names = staar::genotype::read_sample_names(&config.genotypes)?;
    let n_samples = sample_names.len();

    let ann_rows = parquet_row_count(&config.annotations);
    let est_rare = (ann_rows as f64 * 0.02) as u64;

    let bytes_per_variant = (n_samples as u64) * 8;
    let chr1_geno = (est_rare as f64 * 0.10) as u64 * bytes_per_variant;
    let overhead = 4 * GB;
    let recommended = chr1_geno + overhead;

    let probe_result = store::probe(
        &config.store_dir,
        &config.genotypes,
        &config.annotations,
    );

    let cache_status = if probe_result.hit { "hit" } else { "miss" };
    let store_info = if probe_result.hit {
        // hit == true guarantees manifest is Some
        let m = probe_result.manifest.as_ref().unwrap();
        json!({
            "store_path": probe_result.store_dir.to_string_lossy(),
            "store_variants": m.n_variants,
            "store_samples": m.n_samples,
            "created_at": m.created_at,
        })
    } else {
        json!(null)
    };

    let plan = commands::DryRunPlan {
        command: "staar".into(),
        inputs: json!({
            "genotypes": config.genotypes.to_string_lossy(),
            "genotype_size": commands::file_size(&config.genotypes),
            "annotations": config.annotations.to_string_lossy(),
            "annotation_rows": ann_rows,
            "n_samples": n_samples,
            "estimated_rare_variants": est_rare,
            "trait": config.trait_names[0],
            "maf_cutoff": config.maf_cutoff,
            "cache_status": cache_status,
            "cache": store_info,
        }),
        memory: commands::MemoryEstimate {
            minimum: "4G".into(),
            recommended: commands::human_bytes(if probe_result.hit {
                (16 * GB).max(8 * GB)
            } else {
                recommended.max(8 * GB)
            }),
            minimum_bytes: 4 * GB,
            recommended_bytes: if probe_result.hit {
                16 * GB
            } else {
                recommended.max(8 * GB)
            },
        },
        output_path: config.output_dir.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
}

fn parquet_row_count(path: &std::path::Path) -> i64 {
    let file = match std::fs::File::open(path) {
        Ok(f) => f,
        Err(_) => return 0,
    };
    let reader = match parquet::file::reader::SerializedFileReader::new(file) {
        Ok(r) => r,
        Err(_) => return 0,
    };
    use parquet::file::reader::FileReader;
    reader.metadata().file_metadata().num_rows()
}
