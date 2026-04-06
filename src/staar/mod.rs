pub mod ancestry;
pub mod carrier;
pub mod genotype;
#[cfg(test)]
mod ground_truth_test;
pub mod masks;
pub mod meta;
pub mod model;
pub mod multi;
pub mod output;
pub mod pipeline;
pub mod score;
pub mod score_cache;
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
