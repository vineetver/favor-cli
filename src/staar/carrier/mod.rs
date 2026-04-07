//! Carrier-indexed sparse-native genotype types.
//!
//! Core types (CarrierList, CarrierEntry) represent non-reference genotypes.
//! SparseG (in staar::sparse_g) owns the mmap'd storage; this module provides
//! the data types and scoring logic that operate on carrier lists.

pub mod encoding;
pub mod reader;
pub mod sparse_score;

pub use reader::{CarrierList, VariantIndex, VariantIndexEntry};
pub use sparse_score::AnalysisVectors;
