use std::path::PathBuf;

use serde_json::json;

use crate::commands;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::staar;
use crate::staar::masks::ScangParams;
use crate::staar::pipeline::{CohortSource, StaarConfig, StaarPipeline};
#[cfg(test)]
use crate::staar::MaskCategory;
use crate::staar::RunMode;
use crate::store::cohort::CohortId;
use crate::store::list::VariantSet;

const GB: u64 = 1024 * 1024 * 1024;

pub struct StaarArgs {
    pub genotypes: Vec<PathBuf>,
    pub cohort: Option<String>,
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

fn build_config(args: StaarArgs) -> Result<StaarConfig, CohortError> {
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

    // Multi-trait dispatch. Joint MultiSTAAR covers unrelated continuous
    // traits with a shared covariate matrix; everything else (kinship-aware
    // joint null, SPA for binary, AI-STAAR ancestry weighting, sumstats
    // dump) is a separate track. Reject the combinations at config-build
    // time so the pipeline never sees an invalid state.
    let multi_trait = args.trait_names.len() > 1;
    if multi_trait {
        let n_traits = args.trait_names.len();
        let reject = |flag: &str, why: &str| {
            CohortError::Input(format!(
                "--trait-name carries {n_traits} traits (multi-trait joint STAAR), \
                 which is incompatible with {flag}: {why}"
            ))
        };
        if args.spa {
            return Err(reject(
                "--spa",
                "SPA applies to binary traits; joint MultiSTAAR is unrelated continuous only",
            ));
        }
        if blank_to_none(args.ancestry_col.clone()).is_some() {
            return Err(reject(
                "--ancestry-col",
                "AI-STAAR ancestry weighting has not been composed with the joint multi-trait kernel",
            ));
        }
        if !args.kinship.is_empty() {
            return Err(reject(
                "--kinship",
                "kinship-aware joint null (fit_null_glmmkin_multi) is not yet implemented; \
                 the current path assumes unrelated samples",
            ));
        }
        if blank_to_none(args.kinship_groups.clone()).is_some() {
            return Err(reject(
                "--kinship-groups",
                "kinship-aware joint null is not yet implemented",
            ));
        }
        if args.emit_sumstats {
            return Err(reject(
                "--emit-sumstats",
                "MetaSTAAR sumstats are single-trait; run each trait separately to dump U/K",
            ));
        }
    }

    let cohort_name = blank_to_none(args.cohort.clone());
    let (cohort_source, cohort_id, output_dir) = match cohort_name {
        Some(id) => {
            if !args.genotypes.is_empty() {
                return Err(CohortError::Input(
                    "Pass either --cohort <id> (load pre-built) or --genotypes <vcf> \
                     (probe/build), not both."
                        .into(),
                ));
            }
            let cohort_id = CohortId::new(id);
            let output_dir = args
                .output_path
                .clone()
                .unwrap_or_else(|| PathBuf::from(format!("{}.staar", cohort_id.as_str())));
            (CohortSource::Existing, cohort_id, output_dir)
        }
        None => {
            if args.genotypes.is_empty() {
                return Err(CohortError::Input(
                    "STAAR needs either --cohort <id> (pre-built cohort) or --genotypes <vcf>."
                        .into(),
                ));
            }
            let genotypes = crate::commands::ingest::resolve_inputs(args.genotypes)?;
            let annotations = args.annotations.clone().ok_or_else(|| {
                CohortError::Input(
                    "Without --cohort, STAAR requires --annotations <path> from \
                     `favorannotate`."
                        .into(),
                )
            })?;
            if !annotations.exists() {
                return Err(CohortError::Input(format!(
                    "Annotations not found: {}. Run `favorannotate` first.",
                    annotations.display()
                )));
            }
            let ann_vs = VariantSet::open(&annotations).map_err(|e| {
                CohortError::Input(format!(
                    "'{}' is not a valid variant set ({}). Run `favorannotate` to produce one.",
                    annotations.display(),
                    e,
                ))
            })?;
            ann_vs.require_annotated()?;
            let cohort_id = CohortId::new(
                blank_to_none(args.cohort_id.clone())
                    .unwrap_or_else(|| derive_cohort_id(&genotypes[0])),
            );
            let output_dir = args
                .output_path
                .clone()
                .unwrap_or_else(|| default_output_dir(&genotypes[0]));
            (
                CohortSource::Fresh {
                    genotypes,
                    annotations,
                },
                cohort_id,
                output_dir,
            )
        }
    };

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

    Ok(StaarConfig {
        cohort_source,
        phenotype: args.phenotype,
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
        run_mode: if multi_trait {
            RunMode::MultiTrait
        } else if args.emit_sumstats {
            RunMode::EmitSumstats
        } else {
            RunMode::Analyze
        },
        rebuild_store: args.rebuild_store,
        column_map,
        output_dir,
        cohort_id,
    })
}

use super::derive_cohort_id;

fn parse_column_map(
    entries: &[String],
) -> Result<std::collections::HashMap<String, String>, CohortError> {
    let mut map = std::collections::HashMap::new();
    for entry in entries {
        let (k, v) = entry.split_once('=').ok_or_else(|| {
            CohortError::Input(format!(
                "Invalid --column-map entry '{entry}'. Expected key=value."
            ))
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
    let cohort = engine.cohort(&config.cohort_id);
    let (cache_status, store_info, recommended_bytes, n_samples, est_rare, inputs_json) =
        match (config.genotypes(), config.annotations()) {
            // Existing cohort: probe-by-load. Manifest tells us everything.
            (None, None) => {
                let manifest_path = cohort.dir().join("manifest.json");
                if !manifest_path.exists() {
                    return Err(CohortError::DataMissing(format!(
                        "Cohort '{}' not found at {}. Run `favoringest <vcf> \
                         --annotations <set> --cohort-id {}` first.",
                        config.cohort_id.as_str(),
                        cohort.dir().display(),
                        config.cohort_id.as_str()
                    )));
                }
                let txt = std::fs::read_to_string(&manifest_path).map_err(|e| {
                    CohortError::DataMissing(format!(
                        "Cannot read {}: {e}",
                        manifest_path.display()
                    ))
                })?;
                let m: crate::store::cohort::CohortManifest =
                    serde_json::from_str(&txt).map_err(|e| {
                        CohortError::DataMissing(format!("Cohort manifest unparseable: {e}"))
                    })?;
                let inputs = json!({
                    "cohort_id": config.cohort_id.as_str(),
                    "cohort_dir": cohort.dir().to_string_lossy(),
                    "n_samples": m.n_samples,
                    "n_variants": m.n_variants,
                    "trait": config.trait_names[0],
                    "maf_cutoff": config.maf_cutoff,
                    "cache_status": "hit",
                    "run_mode": format!("{:?}", config.run_mode),
                });
                (
                    "hit",
                    json!({
                        "cohort_id": config.cohort_id.as_str(),
                        "store_path": cohort.dir().to_string_lossy(),
                        "store_variants": m.n_variants,
                        "store_samples": m.n_samples,
                        "created_at": m.created_at,
                    }),
                    16 * GB,
                    m.n_samples,
                    0u64,
                    inputs,
                )
            }
            (Some(genotypes), Some(annotations)) => {
                let sample_names = staar::genotype::read_sample_names(&genotypes[0])?;
                let n_samples = sample_names.len();
                let ann_rows = parquet_row_count(annotations).unwrap_or_else(|e| {
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
                let probe_result = cohort.probe(genotypes, annotations);
                let (status, store_info, mem) = match &probe_result.manifest {
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
                let geno_files: Vec<String> = genotypes.iter()
                    .map(|p| p.to_string_lossy().into())
                    .collect();
                let inputs = json!({
                    "genotypes": geno_files,
                    "genotype_file_count": genotypes.len(),
                    "annotations": annotations.to_string_lossy(),
                    "annotation_rows": ann_rows,
                    "n_samples": n_samples,
                    "estimated_rare_variants": est_rare,
                    "trait": config.trait_names[0],
                    "maf_cutoff": config.maf_cutoff,
                    "cache_status": status,
                    "run_mode": format!("{:?}", config.run_mode),
                });
                (status, store_info, mem, n_samples, est_rare, inputs)
            }
            // (Some, None) and (None, Some) are unreachable: build_config
            // rejects mixed input shapes.
            _ => unreachable!("build_config rejects mixed cohort sources"),
        };

    let runtime_seconds = staar_runtime_seconds(
        n_samples,
        est_rare,
        engine.resources().threads.max(1),
        config.has_kinship(),
    );
    let plan = commands::DryRunPlan {
        command: "staar".into(),
        inputs: serde_json::json!({
            "args": inputs_json,
            "cache": store_info,
            "cache_status": cache_status,
        }),
        memory: commands::MemoryEstimate {
            minimum: "4G".into(),
            recommended: commands::human_bytes(recommended_bytes),
            minimum_bytes: 4 * GB,
            recommended_bytes,
        },
        runtime: Some(commands::RuntimeEstimate::from_seconds(runtime_seconds)),
        output_path: config.output_dir.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
}

/// Coarse wall-clock budget for `favorstaar`.
///
/// Calibrated against profiling on the chr22 fixtures: the score-cache fill
/// dominates and scales linearly with `n_samples × est_rare / threads`. Null
/// model fitting adds a fixed overhead that jumps when kinship pulls in REML
/// or PQL. Output is intentionally a budget, not a benchmark — operators use
/// it to size SLURM allocations, not to predict the wall clock to the second.
fn staar_runtime_seconds(n_samples: usize, est_rare: u64, threads: usize, has_kinship: bool) -> u64 {
    const SCORE_NS_PER_SAMPLE_VARIANT: u64 = 1_000;
    const NULL_BASE_S: u64 = 30;
    const NULL_KINSHIP_S: u64 = 270;

    let score_ns = (n_samples as u64).saturating_mul(est_rare) * SCORE_NS_PER_SAMPLE_VARIANT;
    let score_s = score_ns / 1_000_000_000 / threads as u64;
    let null_s = if has_kinship { NULL_KINSHIP_S } else { NULL_BASE_S };
    let total = null_s.saturating_add(score_s);
    total.max(60)
}

fn parquet_row_count(path: &std::path::Path) -> Result<i64, String> {
    // Annotations may be a directory of parquet files (cohort annotate output)
    // or a single .parquet file. For directories we can't cheaply estimate
    // without scanning all part files, so we report None and let the caller
    // fall back to a safe default.
    if path.is_dir() {
        return Err(format!(
            "{} is a directory; row count not estimated",
            path.display()
        ));
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

    #[test]
    fn staar_runtime_grows_with_workload_and_shrinks_with_threads() {
        let small = staar_runtime_seconds(1_000, 100_000, 8, false);
        let big = staar_runtime_seconds(50_000, 1_000_000, 8, false);
        assert!(big > small);
        let serial = staar_runtime_seconds(50_000, 1_000_000, 1, false);
        let parallel = staar_runtime_seconds(50_000, 1_000_000, 16, false);
        assert!(parallel < serial);
    }

    #[test]
    fn staar_runtime_kinship_adds_overhead() {
        let no_kin = staar_runtime_seconds(10_000, 100_000, 8, false);
        let with_kin = staar_runtime_seconds(10_000, 100_000, 8, true);
        assert!(with_kin > no_kin);
    }

    #[test]
    fn staar_runtime_minimum_floor() {
        assert_eq!(staar_runtime_seconds(0, 0, 8, false), 60);
    }

    /// Baseline `StaarArgs` with a real-on-disk phenotype file and a
    /// pre-built cohort source. Covers the minimum set of fields every
    /// multi-trait conflict test wants to override.
    ///
    /// `--cohort <id>` picks the `Existing` path so the fixture does not
    /// need a VCF or annotation directory — the multi-trait conflict block
    /// runs before the cohort branch, which is what we want to exercise.
    fn multi_trait_args(pheno_path: PathBuf, traits: &[&str]) -> StaarArgs {
        StaarArgs {
            genotypes: Vec::new(),
            cohort: Some("dummy".into()),
            phenotype: pheno_path,
            trait_names: traits.iter().map(|s| (*s).into()).collect(),
            covariates: vec!["age".into(), "sex".into()],
            annotations: None,
            masks: vec!["coding".into()],
            maf_cutoff: 0.01,
            window_size: 2000,
            individual: false,
            spa: false,
            ancestry_col: None,
            ai_base_tests: 5,
            ai_seed: 7590,
            scang_lmin: 40,
            scang_lmax: 300,
            scang_step: 10,
            kinship: Vec::new(),
            kinship_groups: None,
            known_loci: None,
            emit_sumstats: false,
            rebuild_store: false,
            column_map: Vec::new(),
            output_path: None,
            cohort_id: None,
        }
    }

    fn pheno_file() -> tempfile::NamedTempFile {
        tempfile::NamedTempFile::new().expect("create temp phenotype file")
    }

    fn assert_input_error_mentions(result: Result<StaarConfig, CohortError>, needle: &str) {
        match result {
            Ok(_) => panic!("expected Input error mentioning {needle:?}, got Ok"),
            Err(CohortError::Input(msg)) => assert!(
                msg.contains(needle),
                "expected error mentioning {needle:?}, got: {msg}"
            ),
            Err(other) => panic!("expected Input error, got {other:?}"),
        }
    }

    #[test]
    fn multi_trait_activates_when_two_traits() {
        let pheno = pheno_file();
        let args = multi_trait_args(pheno.path().to_path_buf(), &["BMI", "HEIGHT"]);
        let cfg = build_config(args).expect("multi-trait alone should be accepted");
        assert_eq!(cfg.run_mode, RunMode::MultiTrait);
        assert_eq!(cfg.trait_names, vec!["BMI".to_string(), "HEIGHT".into()]);
    }

    #[test]
    fn single_trait_stays_analyze() {
        let pheno = pheno_file();
        let args = multi_trait_args(pheno.path().to_path_buf(), &["BMI"]);
        let cfg = build_config(args).expect("single-trait baseline should be accepted");
        assert_eq!(cfg.run_mode, RunMode::Analyze);
    }

    #[test]
    fn multi_trait_rejects_spa() {
        let pheno = pheno_file();
        let mut args = multi_trait_args(pheno.path().to_path_buf(), &["BMI", "HEIGHT"]);
        args.spa = true;
        assert_input_error_mentions(build_config(args), "--spa");
    }

    #[test]
    fn multi_trait_rejects_ancestry_col() {
        let pheno = pheno_file();
        let mut args = multi_trait_args(pheno.path().to_path_buf(), &["BMI", "HEIGHT"]);
        args.ancestry_col = Some("super_population".into());
        assert_input_error_mentions(build_config(args), "--ancestry-col");
    }

    #[test]
    fn multi_trait_rejects_kinship() {
        let pheno = pheno_file();
        let mut args = multi_trait_args(pheno.path().to_path_buf(), &["BMI", "HEIGHT"]);
        args.kinship = vec![PathBuf::from("/no/such/kin.tsv")];
        assert_input_error_mentions(build_config(args), "--kinship");
    }

    #[test]
    fn multi_trait_rejects_kinship_groups() {
        let pheno = pheno_file();
        let mut args = multi_trait_args(pheno.path().to_path_buf(), &["BMI", "HEIGHT"]);
        args.kinship_groups = Some("ancestry".into());
        assert_input_error_mentions(build_config(args), "--kinship-groups");
    }

    #[test]
    fn multi_trait_rejects_emit_sumstats() {
        let pheno = pheno_file();
        let mut args = multi_trait_args(pheno.path().to_path_buf(), &["BMI", "HEIGHT"]);
        args.emit_sumstats = true;
        assert_input_error_mentions(build_config(args), "--emit-sumstats");
    }
}
