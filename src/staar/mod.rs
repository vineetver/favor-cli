pub mod ancestry;
pub mod carrier;
pub mod genotype;
#[cfg(test)]
mod ground_truth_test;
pub mod kinship;
pub mod masks;
pub mod meta;
pub mod model;
pub mod multi;
pub mod output;
pub mod pipeline;
pub mod run_manifest;
pub mod score;
pub mod score_cache;
pub mod scoring;
pub mod sparse_g;
pub mod sparse_g_writer;
pub mod stats;
pub mod store;
#[allow(dead_code)]
pub mod store_validate;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitType {
    Continuous,
    Binary,
}

/// What this `favor staar` invocation actually does end-to-end.
///
/// Replaces the old `emit_sumstats: bool` flag scattered through the
/// pipeline. Each variant routes to a different terminal stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunMode {
    /// Full STAAR analysis: score tests, omnibus, results, report.
    Analyze,
    /// Stop after the null model and dump per-variant U/K for MetaSTAAR.
    EmitSumstats,
}

/// Per-gene scoring backend chosen for this run.
///
/// Picked once at config-time from `(spa, ancestry_col, trait_type)` and
/// then used as a single dispatch value inside the gene scorer instead of
/// re-deriving the boolean combinations at every call site.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringMode {
    /// Cached U/K → summary-stat tests. Default path.
    Standard,
    /// Saddlepoint approximation. Binary traits only.
    Spa,
    /// AI-STAAR ensemble across population groups.
    AiStaar,
}

/// Which null-model fitter to dispatch into.
///
/// Picked from `(trait_type, kinship_present)` at config time. The fitting
/// stage just matches on this — no more `if !kinships.is_empty()` chains.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NullModelKind {
    /// Ordinary least squares for continuous traits.
    Glm,
    /// IRLS logistic regression for binary traits.
    Logistic,
    /// Variance-component REML with kinship for continuous traits.
    KinshipReml,
    /// Penalised quasi-likelihood GLMM for binary traits with kinship.
    KinshipPql,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MaskCategory {
    Coding,
    Noncoding,
    SlidingWindow,
    Scang,
    Custom,
}

impl std::str::FromStr for MaskCategory {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "coding" => Ok(Self::Coding),
            "noncoding" => Ok(Self::Noncoding),
            "sliding-window" | "window" => Ok(Self::SlidingWindow),
            "scang" => Ok(Self::Scang),
            "custom" => Ok(Self::Custom),
            _ => Err(format!("unknown mask: {s}")),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum MaskType {
    PLof,
    Missense,
    DisruptiveMissense,
    PLofMissense,
    Synonymous,
    Ptv,
    PtvDs,
    Upstream,
    Downstream,
    Utr,
    PromoterCage,
    PromoterDhs,
    EnhancerCage,
    EnhancerDhs,
    Ncrna,
    SlidingWindow,
    Scang,
}

impl MaskType {
    pub fn file_stem(&self) -> String {
        match self {
            Self::PLof => "coding_pLoF".into(),
            Self::Missense => "coding_missense".into(),
            Self::DisruptiveMissense => "coding_disruptive_missense".into(),
            Self::PLofMissense => "coding_pLoF_missense".into(),
            Self::Synonymous => "coding_synonymous".into(),
            Self::Ptv => "coding_ptv".into(),
            Self::PtvDs => "coding_ptv_ds".into(),
            Self::Upstream => "noncoding_upstream".into(),
            Self::Downstream => "noncoding_downstream".into(),
            Self::Utr => "noncoding_utr".into(),
            Self::PromoterCage => "noncoding_promoter_CAGE".into(),
            Self::PromoterDhs => "noncoding_promoter_DHS".into(),
            Self::EnhancerCage => "noncoding_enhancer_CAGE".into(),
            Self::EnhancerDhs => "noncoding_enhancer_DHS".into(),
            Self::Ncrna => "noncoding_ncRNA".into(),
            Self::SlidingWindow => "sliding_window".into(),
            Self::Scang => "scang".into(),
        }
    }
}

use crate::types::Chromosome;

#[derive(Debug, Clone)]
pub struct GeneResult {
    pub ensembl_id: String,
    pub gene_symbol: String,
    pub chromosome: Chromosome,
    pub start: u32,
    pub end: u32,
    pub n_variants: u32,
    pub cumulative_mac: u32,
    pub staar: score::StaarResult,
}
