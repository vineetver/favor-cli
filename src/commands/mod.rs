pub mod annotate;
pub mod annotation;
pub mod enrich;
pub mod ingest;
pub mod inspect;
pub mod interpret;
pub mod meta_staar;
pub mod staar;
pub mod store;

use std::path::{Path, PathBuf};

use serde::Serialize;

use crate::config::Tier;
use crate::error::CohortError;
use crate::output::Output;
use crate::staar::MaskCategory;

pub fn parse_mask_categories(masks: &[String]) -> Result<Vec<MaskCategory>, CohortError> {
    masks
        .iter()
        .map(|s| {
            s.parse::<MaskCategory>().map_err(|_| {
                CohortError::Input(format!(
                    "Unknown mask '{s}'. Available: coding, noncoding, sliding-window, scang"
                ))
            })
        })
        .collect()
}

pub struct IngestConfig {
    pub inputs: Vec<PathBuf>,
    pub output: PathBuf,
    pub emit_sql: bool,
    pub build_override: Option<crate::cli::GenomeBuild>,
    pub annotations: Option<PathBuf>,
    pub cohort_id: Option<String>,
    pub rebuild: bool,
    pub chromosome_filter: Option<crate::types::ChromosomeSet>,
}

pub fn derive_cohort_id(genotypes: &Path) -> String {
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

pub struct AnnotateConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub tier: Tier,
}

pub struct EnrichConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub tissue_name: String,
}

pub struct MetaStaarConfig {
    pub study_dirs: Vec<PathBuf>,
    pub mask_categories: Vec<MaskCategory>,
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub known_loci: Option<PathBuf>,
    pub conditional_model: crate::cli::ConditionalModel,
    pub output_dir: PathBuf,
}

const GB: u64 = 1024 * 1024 * 1024;

#[derive(Serialize)]
pub struct DryRunPlan {
    pub command: String,
    pub inputs: serde_json::Value,
    pub memory: MemoryEstimate,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub runtime: Option<RuntimeEstimate>,
    pub output_path: String,
}

#[derive(Serialize)]
pub struct MemoryEstimate {
    pub minimum: String,
    pub recommended: String,
    pub minimum_bytes: u64,
    pub recommended_bytes: u64,
}

impl MemoryEstimate {
    pub fn default_estimate() -> Self {
        Self {
            minimum: "4G".into(),
            recommended: "16G".into(),
            minimum_bytes: 4 * GB,
            recommended_bytes: 16 * GB,
        }
    }
}

/// Coarse wall-clock budget for planning, not a benchmark prediction.
/// Callers compute `seconds` from a workload-specific formula and hand the
/// number to `RuntimeEstimate::from_seconds`, which produces the human form.
#[derive(Serialize)]
pub struct RuntimeEstimate {
    pub seconds: u64,
    pub human: String,
}

impl RuntimeEstimate {
    pub fn from_seconds(seconds: u64) -> Self {
        Self {
            seconds,
            human: human_seconds(seconds),
        }
    }
}

fn human_seconds(s: u64) -> String {
    if s < 60 {
        format!("{s}s")
    } else if s < 3600 {
        format!("{}m{:02}s", s / 60, s % 60)
    } else if s < 86_400 {
        format!("{}h{:02}m", s / 3600, (s % 3600) / 60)
    } else {
        format!("{}d{:02}h", s / 86_400, (s % 86_400) / 3600)
    }
}

pub fn emit(plan: &DryRunPlan, out: &dyn Output) {
    out.result_json(&serde_json::to_value(plan).unwrap_or_default());
}

/// Default output path for a command that transforms an input file in
/// place: strip any of `strip_suffixes` from the file name (first match
/// wins), then append `append` and rejoin to the parent directory.
/// Replaces two hand-rolled copies that used to live in `annotate` and
/// `enrich`. `ingest` derives its output via `file_stem` plus a build-
/// tag split and stays inline.
pub fn derive_output_path(input: &Path, strip_suffixes: &[&str], append: &str) -> PathBuf {
    let name = input.file_name().unwrap_or_default().to_string_lossy();
    let stem: &str = strip_suffixes
        .iter()
        .find_map(|suf| name.strip_suffix(*suf))
        .unwrap_or(&name);
    input
        .parent()
        .unwrap_or(input)
        .join(format!("{stem}{append}"))
}

pub use crate::data::transfer::human_size as human_bytes;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn human_seconds_buckets() {
        assert_eq!(human_seconds(0), "0s");
        assert_eq!(human_seconds(45), "45s");
        assert_eq!(human_seconds(60), "1m00s");
        assert_eq!(human_seconds(125), "2m05s");
        assert_eq!(human_seconds(3_600), "1h00m");
        assert_eq!(human_seconds(7_265), "2h01m");
        assert_eq!(human_seconds(90_061), "1d01h");
    }

    #[test]
    fn runtime_estimate_round_trips_through_serde() {
        let est = RuntimeEstimate::from_seconds(125);
        assert_eq!(est.seconds, 125);
        assert_eq!(est.human, "2m05s");
        let json = serde_json::to_string(&est).unwrap();
        assert!(json.contains("\"seconds\":125"));
        assert!(json.contains("\"human\":\"2m05s\""));
    }

    #[test]
    fn dry_run_plan_omits_runtime_when_none() {
        let plan = DryRunPlan {
            command: "annotate".into(),
            inputs: serde_json::json!({}),
            memory: MemoryEstimate::default_estimate(),
            runtime: None,
            output_path: "/tmp/out".into(),
        };
        let json = serde_json::to_string(&plan).unwrap();
        assert!(!json.contains("runtime"));
    }
}
