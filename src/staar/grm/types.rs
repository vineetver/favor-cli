//! Data types for the FastSparseGRM pipeline.

use faer::Mat;

/// Parsed KING .seg row after degree filtering.
#[derive(Clone, Debug)]
pub struct KingSegEntry {
    pub id1: String,
    pub id2: String,
    pub prop_ibd: f64,
    #[allow(dead_code)]
    pub inf_type: RelatednessType,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum RelatednessType {
    Unrelated,
    Fourth,
    Third,
    Second,
    First,
    Dup,
}

impl RelatednessType {
    pub fn from_king_label(s: &str) -> Self {
        match s.trim() {
            "Dup/MZ" | "Dup" | "MZ" => Self::Dup,
            "PO" | "FS" | "1st" => Self::First,
            "2nd" | "HS" => Self::Second,
            "3rd" => Self::Third,
            "4th" => Self::Fourth,
            _ => Self::Unrelated,
        }
    }

    pub fn degree(self) -> u8 {
        match self {
            Self::Dup => 0,
            Self::First => 1,
            Self::Second => 2,
            Self::Third => 3,
            Self::Fourth => 4,
            Self::Unrelated => 255,
        }
    }
}

/// Connected component of related individuals.
#[derive(Clone, Debug)]
pub struct RelatedComponent {
    pub members: Vec<usize>,
    #[allow(dead_code)]
    pub pairs: Vec<(usize, usize)>,
}

/// Output of the unrelated selection step.
pub struct UnrelatedSubset {
    pub sample_indices: Vec<usize>,
}

/// PCA scores for all samples.
pub struct PcaScores {
    /// (n_samples, n_pcs) column-major faer matrix.
    pub scores: Mat<f64>,
    /// Squared singular values (variance explained per PC).
    pub eigenvalues: Vec<f64>,
}

/// Per-pair kinship accumulator. Numerator and denominator are summed
/// independently across SNP blocks and chromosomes; final kinship is
/// num / den after all blocks are processed.
#[derive(Clone, Debug)]
pub struct KinshipAccum {
    pub idx_i: usize,
    pub idx_j: usize,
    pub numerator: f64,
    pub denominator: f64,
}

/// Final sparse GRM output: symmetric triplets + sample ordering.
pub struct SparseGrm {
    pub triplets: Vec<(usize, usize, f64)>,
    #[allow(dead_code)]
    pub n_samples: usize,
}

/// Combined artifact from a full FastSparseGRM run.
pub struct GrmArtifact {
    pub grm: SparseGrm,
    pub pca: PcaScores,
    pub unrelated: UnrelatedSubset,
    pub sample_ids: Vec<String>,
}

/// Configuration for a GRM build.
#[allow(dead_code)]
pub struct GrmConfig {
    pub cohort_id: crate::store::cohort::CohortId,
    pub king_seg_path: std::path::PathBuf,
    pub degree: u8,
    pub n_pcs: usize,
    pub block_size: usize,
    pub output_dir: Option<std::path::PathBuf>,
}
