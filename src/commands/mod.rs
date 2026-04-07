pub mod annotate;
pub mod enrich;
pub mod ingest;
pub mod inspect;
pub mod interpret;
pub mod meta_staar;
pub mod staar;

use std::path::{Path, PathBuf};

use serde::Serialize;

use crate::config::Tier;
use crate::error::CohortError;
use crate::output::Output;
use crate::staar::MaskCategory;

/// Parse `--masks` strings into the typed enum. Shared between
/// `cohort staar` and `cohort meta-staar` so the accepted set stays in sync.
pub fn parse_mask_categories(masks: &[String]) -> Result<Vec<MaskCategory>, CohortError> {
    masks
        .iter()
        .map(|s| {
            s.parse::<MaskCategory>().map_err(|_| {
                CohortError::Input(format!(
                    "Unknown mask '{s}'. Available: coding, noncoding, sliding-window, scang, custom"
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
}

pub struct AnnotateConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub tier: Tier,
    pub data_root: PathBuf,
}

pub struct EnrichConfig {
    pub input: PathBuf,
    pub output: PathBuf,
    pub tissue_name: String,
    pub tissue_dir: PathBuf,
}

pub struct MetaStaarConfig {
    pub study_dirs: Vec<PathBuf>,
    pub mask_categories: Vec<MaskCategory>,
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub output_dir: PathBuf,
}

const GB: u64 = 1024 * 1024 * 1024;

#[derive(Serialize)]
pub struct DryRunPlan {
    pub command: String,
    pub inputs: serde_json::Value,
    pub memory: MemoryEstimate,
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
    /// Default estimate: engine will spill to disk if needed.
    pub fn default_estimate() -> Self {
        Self {
            minimum: "4G".into(),
            recommended: "16G".into(),
            minimum_bytes: 4 * GB,
            recommended_bytes: 16 * GB,
        }
    }
}

pub fn emit(plan: &DryRunPlan, out: &dyn Output) {
    out.result_json(&serde_json::to_value(plan).unwrap_or_default());
}

pub fn file_size(path: &Path) -> u64 {
    std::fs::metadata(path).map(|m| m.len()).unwrap_or(0)
}

pub use crate::data::transfer::human_size as human_bytes;
