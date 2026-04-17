pub mod ancestry;
pub mod carrier;
pub mod genotype;
pub mod grm;
#[cfg(test)]
mod ground_truth_test;
#[cfg(test)]
mod invariance_test;
pub mod kinship;
pub mod ld_prune;
pub mod masks;
pub mod meta;
pub mod model;
pub mod multi;
pub mod output;
pub mod pipeline;
pub mod run_manifest;
pub mod scang;
pub mod score;
pub mod scoring;
pub mod stats;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TraitType {
    Continuous,
    Binary,
}

/// What this `favorstaar` invocation actually does end-to-end.
///
/// Replaces the old `emit_sumstats: bool` flag scattered through the
/// pipeline. Each variant routes to a different terminal stage.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum RunMode {
    /// Full single-trait STAAR analysis: score tests, omnibus, results, report.
    Analyze,
    /// Stop after the null model and dump per-variant U/K for MetaSTAAR.
    EmitSumstats,
    /// Joint multi-trait STAAR for unrelated continuous traits with shared
    /// covariates. Activated when `--trait-name` carries more than one name.
    /// Rejects `--spa`, `--ancestry-col`, `--kinship`, and `--emit-sumstats`
    /// at config build time; kinship-aware joint nulls are a separate track.
    MultiTrait,
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
}

impl std::str::FromStr for MaskCategory {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "coding" => Ok(Self::Coding),
            "noncoding" => Ok(Self::Noncoding),
            "sliding-window" | "window" => Ok(Self::SlidingWindow),
            "scang" => Ok(Self::Scang),
            _ => Err(format!("unknown mask: {s}")),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum MaskType {
    PLof,
    Missense,
    DisruptiveMissense,
    PLofDs,
    Synonymous,
    Ptv,
    PtvDs,
    CodingInclPtv,
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
    pub fn file_stem(self) -> &'static str {
        match self {
            Self::PLof => "coding_pLoF",
            Self::Missense => "coding_missense",
            Self::DisruptiveMissense => "coding_disruptive_missense",
            Self::PLofDs => "coding_plof_ds",
            Self::Synonymous => "coding_synonymous",
            Self::Ptv => "coding_ptv",
            Self::PtvDs => "coding_ptv_ds",
            Self::CodingInclPtv => "coding_incl_ptv",
            Self::Upstream => "noncoding_upstream",
            Self::Downstream => "noncoding_downstream",
            Self::Utr => "noncoding_utr",
            Self::PromoterCage => "noncoding_promoter_CAGE",
            Self::PromoterDhs => "noncoding_promoter_DHS",
            Self::EnhancerCage => "noncoding_enhancer_CAGE",
            Self::EnhancerDhs => "noncoding_enhancer_DHS",
            Self::Ncrna => "noncoding_ncRNA",
            Self::SlidingWindow => "sliding_window",
            Self::Scang => "scang",
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
    /// SCANG-O empirical −log10(p) threshold at α = 0.05, NaN otherwise.
    /// Emitted alongside per-window p-values so operators can cross-check
    /// `-log10(p) > emthr` matches the R `SCANG_O_res$th0` gate. See
    /// `crate::staar::scang::chrom_threshold` and SCANG R/SCANG.r:181-205.
    pub emthr: f64,
}
