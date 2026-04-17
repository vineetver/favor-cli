//! FastSparseGRM: sparse ancestry-adjusted GRM builder.
//!
//! Mirrors the R FastSparseGRM package (Lin & Dey, Nature Genetics 2024).
//! Produces a sparse kinship matrix + PCA scores from a cohort's genotype
//! store and KING IBD segment output.

pub mod cache;
pub mod estimate;
pub mod king;
pub mod pca;
pub mod types;
pub mod unrelated;
