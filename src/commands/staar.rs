//! `cohort staar` command — thin orchestrator.
//!
//! Validates `StaarArgs`, builds a typed `StaarConfig`, and delegates to
//! `StaarPipeline::run`. Replaces the previous 20-positional-arg entry point.

use std::path::PathBuf;

use serde_json::json;

use crate::commands;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::staar;
use crate::staar::masks::ScangParams;
#[cfg(test)]
use crate::staar::MaskCategory;
use crate::staar::pipeline::{StaarConfig, StaarPipeline};
use crate::store::cohort::CohortId;
use crate::store::list::VariantSet;

const GB: u64 = 1024 * 1024 * 1024;

/// Raw, unvalidated CLI arguments for `cohort staar`.
///
/// Mirrors the `Command::Staar` variant in `cli.rs`. Carrying a single
/// struct around lets each layer pass arguments by name, kills the 20-arg
/// signature blob, and gives validation a typed home.
pub struct StaarArgs {
    pub genotypes: PathBuf,
    pub phenotype: PathBuf,
    pub trait_names: Vec<String>,
    pub covariates: Vec<String>,
    pub annotations: Option<PathBuf>,
    pub masks: Vec<String>,
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub individual: bool,
    pub spa: bool,
    pub ancestry_col: Option<String>,
    pub ai_base_tests: usize,
    pub ai_seed: u64,
    pub scang_lmin: usize,
    pub scang_lmax: usize,
    pub scang_step: usize,
    pub kinship: Vec<PathBuf>,
    pub kinship_groups: Option<String>,
    pub known_loci: Option<PathBuf>,
    pub emit_sumstats: bool,
    pub rebuild_store: bool,
    pub column_map: Vec<String>,
    pub output_path: Option<PathBuf>,
    pub cohort_id: Option<String>,
}

pub fn run(
    engine: &Engine,
    args: StaarArgs,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    let config = build_config(args)?;
    if dry_run {
        return emit_dry_run(engine, &config, out);
    }
    StaarPipeline::new(config, engine, out)?.run()
}

/// Validate `StaarArgs` and produce a `StaarConfig`. All input checks live
/// here so the pipeline never sees a half-validated config.
fn build_config(args: StaarArgs) -> Result<StaarConfig, CohortError> {
    if !args.genotypes.exists() {
        return Err(CohortError::Input(format!(
            "Genotype VCF not found: '{}'. Check the path to your multi-sample VCF.",
            args.genotypes.display()
        )));
    }
    if !args.phenotype.exists() {
        return Err(CohortError::Input(format!(
            "Phenotype file not found: '{}'. Provide a tab-delimited file with sample IDs, traits, and covariates.",
            args.phenotype.display()
        )));
    }
    if args.maf_cutoff <= 0.0 || args.maf_cutoff >= 0.5 {
        return Err(CohortError::Input(format!(
            "MAF cutoff must be in (0, 0.5), got {}",
            args.maf_cutoff
        )));
    }
    if args.trait_names.is_empty() {
        return Err(CohortError::Input(
            "At least one --trait-name is required.".into(),
        ));
    }

    let annotations = args.annotations.ok_or_else(|| {
        CohortError::Input("STAAR requires --annotations <path> from `cohort annotate`.".into())
    })?;
    if !annotations.exists() {
        return Err(CohortError::Input(format!(
            "Annotations not found: {}. Run `cohort annotate` first.",
            annotations.display()
        )));
    }
    let ann_vs = VariantSet::open(&annotations).map_err(|e| {
        CohortError::Input(format!(
            "'{}' is not a valid variant set ({}). Run `cohort annotate` to produce one.",
            annotations.display(),
            e,
        ))
    })?;
    ann_vs.require_annotated()?;

    let mask_categories = commands::parse_mask_categories(&args.masks)?;

    if let Some(ref loci) = args.known_loci {
        if !loci.exists() {
            return Err(CohortError::Input(format!(
                "Known loci file not found: {}",
                loci.display()
            )));
        }
    }

    for kp in &args.kinship {
        if !kp.exists() {
            return Err(CohortError::DataMissing(format!(
                "Kinship file not found: '{}'",
                kp.display()
            )));
        }
    }

    let column_map = parse_column_map(&args.column_map)?;
    let output_dir = args.output_path.unwrap_or_else(|| default_output_dir(&args.genotypes));
    let cohort_id = CohortId::new(
        blank_to_none(args.cohort_id).unwrap_or_else(|| derive_cohort_id(&args.genotypes)),
    );

    Ok(StaarConfig {
        genotypes: args.genotypes,
        phenotype: args.phenotype,
        annotations,
        trait_names: args.trait_names,
        covariates: args.covariates,
        mask_categories,
        maf_cutoff: args.maf_cutoff,
        window_size: args.window_size,
        individual: args.individual,
        spa: args.spa,
        ancestry_col: blank_to_none(args.ancestry_col),
        ai_base_tests: args.ai_base_tests,
        ai_seed: args.ai_seed,
        scang_params: ScangParams {
            lmin: args.scang_lmin,
            lmax: args.scang_lmax,
            step: args.scang_step,
        },
        kinship: args.kinship,
        kinship_groups: blank_to_none(args.kinship_groups),
        known_loci: args.known_loci,
        emit_sumstats: args.emit_sumstats,
        rebuild_store: args.rebuild_store,
        column_map,
        output_dir,
        cohort_id,
    })
}

/// Default cohort id derived from the input VCF stem when `--cohort-id`
/// is not given. Strips `.vcf`/`.vcf.gz` extensions and replaces every
/// run of non-alphanumeric chars with `_` so the result is a safe
/// directory name.
fn derive_cohort_id(genotypes: &std::path::Path) -> String {
    let stem = genotypes
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("cohort");
    let stem = stem
        .strip_suffix(".vcf.gz")
        .or_else(|| stem.strip_suffix(".vcf.bgz"))
        .or_else(|| stem.strip_suffix(".vcf"))
        .or_else(|| stem.strip_suffix(".bcf"))
        .unwrap_or(stem);
    let mut out = String::with_capacity(stem.len());
    let mut last_underscore = false;
    for ch in stem.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch);
            last_underscore = false;
        } else if !last_underscore && !out.is_empty() {
            out.push('_');
            last_underscore = true;
        }
    }
    while out.ends_with('_') {
        out.pop();
    }
    if out.is_empty() {
        "cohort".into()
    } else {
        out
    }
}

fn parse_column_map(entries: &[String]) -> Result<std::collections::HashMap<String, String>, CohortError> {
    let mut map = std::collections::HashMap::new();
    for entry in entries {
        let (k, v) = entry.split_once('=').ok_or_else(|| {
            CohortError::Input(format!("Invalid --column-map entry '{entry}'. Expected key=value."))
        })?;
        map.insert(k.trim().to_string(), v.trim().to_string());
    }
    Ok(map)
}

fn default_output_dir(genotypes: &std::path::Path) -> PathBuf {
    let stem = genotypes.file_stem().unwrap_or_default().to_string_lossy();
    genotypes.with_file_name(format!("{stem}.staar"))
}

fn blank_to_none(s: Option<String>) -> Option<String> {
    s.and_then(|s| {
        let trimmed = s.trim().to_string();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed)
        }
    })
}

fn emit_dry_run(
    engine: &Engine,
    config: &StaarConfig,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let sample_names = staar::genotype::read_sample_names(&config.genotypes)?;
    let n_samples = sample_names.len();

    let ann_rows = parquet_row_count(&config.annotations).unwrap_or_else(|e| {
        out.warn(&format!(
            "  dry-run: could not count annotation rows ({e}); estimating with 0"
        ));
        0
    });
    let est_rare = (ann_rows as f64 * 0.02) as u64;

    let bytes_per_variant = (n_samples as u64) * 8;
    let chr1_geno = (est_rare as f64 * 0.10) as u64 * bytes_per_variant;
    let overhead = 4 * GB;
    let recommended = chr1_geno + overhead;

    let cohort = engine.cohort(&config.cohort_id);
    let probe_result = cohort.probe(&config.genotypes, &config.annotations);

    let (cache_status, store_info, recommended_bytes) = match &probe_result.manifest {
        Some(m) => (
            "hit",
            json!({
                "cohort_id": config.cohort_id.as_str(),
                "store_path": probe_result.store_dir.to_string_lossy(),
                "store_variants": m.n_variants,
                "store_samples": m.n_samples,
                "created_at": m.created_at,
            }),
            16 * GB,
        ),
        None => ("miss", json!(null), recommended.max(8 * GB)),
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
            "run_mode": format!("{:?}", config.run_mode()),
        }),
        memory: commands::MemoryEstimate {
            minimum: "4G".into(),
            recommended: commands::human_bytes(recommended_bytes),
            minimum_bytes: 4 * GB,
            recommended_bytes,
        },
        output_path: config.output_dir.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
}

fn parquet_row_count(path: &std::path::Path) -> Result<i64, String> {
    // Annotations may be a directory of parquet files (cohort annotate output)
    // or a single .parquet file. For directories we can't cheaply estimate
    // without scanning all part files, so we report None and let the caller
    // fall back to a safe default.
    if path.is_dir() {
        return Err(format!("{} is a directory; row count not estimated", path.display()));
    }
    let file = std::fs::File::open(path).map_err(|e| format!("open {}: {e}", path.display()))?;
    let reader = parquet::file::reader::SerializedFileReader::new(file)
        .map_err(|e| format!("parse parquet metadata for {}: {e}", path.display()))?;
    use parquet::file::reader::FileReader;
    Ok(reader.metadata().file_metadata().num_rows())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_mask_categories_ok() {
        let m = commands::parse_mask_categories(&["coding".into(), "noncoding".into()]).unwrap();
        assert_eq!(m, vec![MaskCategory::Coding, MaskCategory::Noncoding]);
    }

    #[test]
    fn parse_mask_categories_unknown_rejected() {
        let err = commands::parse_mask_categories(&["bogus".into()]).unwrap_err();
        match err {
            CohortError::Input(msg) => assert!(msg.contains("bogus")),
            _ => panic!("expected Input error"),
        }
    }

    #[test]
    fn parse_column_map_round_trip() {
        let m = parse_column_map(&["id=IID".into(), "trait=BMI".into()]).unwrap();
        assert_eq!(m.get("id"), Some(&"IID".to_string()));
        assert_eq!(m.get("trait"), Some(&"BMI".to_string()));
    }

    #[test]
    fn parse_column_map_rejects_bad_entry() {
        assert!(parse_column_map(&["no_equals_sign".into()]).is_err());
    }

    #[test]
    fn blank_to_none_strips_empty() {
        assert_eq!(blank_to_none(None), None);
        assert_eq!(blank_to_none(Some(String::new())), None);
        assert_eq!(blank_to_none(Some("   ".into())), None);
        assert_eq!(blank_to_none(Some(" foo ".into())), Some("foo".into()));
    }

    #[test]
    fn default_output_dir_appends_staar_suffix() {
        // file_stem only strips the last extension, so .vcf.gz → .vcf → .vcf.staar.
        let p = default_output_dir(std::path::Path::new("/data/cohort.vcf.gz"));
        assert_eq!(p, PathBuf::from("/data/cohort.vcf.staar"));
    }

    #[test]
    fn derive_cohort_id_strips_vcf_extensions_and_sanitizes() {
        assert_eq!(
            derive_cohort_id(std::path::Path::new("/data/ukb_chr22.vcf.gz")),
            "ukb_chr22"
        );
        assert_eq!(
            derive_cohort_id(std::path::Path::new("/data/run-2026-04.vcf")),
            "run_2026_04"
        );
        assert_eq!(
            derive_cohort_id(std::path::Path::new("/data/sample.bcf")),
            "sample"
        );
        assert_eq!(derive_cohort_id(std::path::Path::new("/")), "cohort");
    }
}
