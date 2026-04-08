//! Engine-facing newtypes for the cohort layer.
//!
//! Distinct from the underlying `u32` so the resolver cannot mix
//! per-cohort dense indices with raw sample/position counters at the type
//! level. `repr(transparent)` keeps `&[VariantVcf]` cast-equivalent to
//! `&[u32]` for the sparse_g.bin reader, which still speaks `u32`.

use crate::types::Chromosome;

use super::variants::VariantIndexEntry;
use super::CohortId;

#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct VariantVcf(pub u32);

impl VariantVcf {
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

impl From<u32> for VariantVcf {
    #[inline]
    fn from(v: u32) -> Self {
        Self(v)
    }
}

#[repr(transparent)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SampleIdx(pub u32);

impl SampleIdx {
    #[inline]
    pub fn get(self) -> u32 {
        self.0
    }
}

/// `&[VariantVcf]` ↔ `&[u32]` zero-cost cast (both `repr(transparent)` over `u32`).
#[inline]
pub(crate) fn as_u32_slice(vcfs: &[VariantVcf]) -> &[u32] {
    // SAFETY: VariantVcf is `#[repr(transparent)]` over `u32`.
    unsafe { std::slice::from_raw_parts(vcfs.as_ptr() as *const u32, vcfs.len()) }
}

#[inline]
pub(crate) fn from_u32_slice(vcfs: &[u32]) -> &[VariantVcf] {
    // SAFETY: VariantVcf is `#[repr(transparent)]` over `u32`.
    unsafe { std::slice::from_raw_parts(vcfs.as_ptr() as *const VariantVcf, vcfs.len()) }
}

/// Sorted, dedup'd, scoped to one (cohort, chromosome). The cohort+chrom
/// fields are metadata, not enforced at every callsite — `ChromosomeView`
/// methods debug-assert sortedness when they accept a slice.
#[derive(Debug, Clone)]
pub struct SortedVcfs {
    cohort: CohortId,
    chrom: Chromosome,
    vcfs: Vec<VariantVcf>,
}

impl SortedVcfs {
    /// Caller guarantees `vcfs` is sorted ascending and dedup'd. Debug
    /// builds verify; release builds trust.
    pub fn new(cohort: CohortId, chrom: Chromosome, vcfs: Vec<VariantVcf>) -> Self {
        debug_assert!(
            vcfs.windows(2).all(|w| w[0] < w[1]),
            "SortedVcfs requires strictly ascending input"
        );
        Self { cohort, chrom, vcfs }
    }

    /// Take an unsorted Vec, sort+dedup it, then wrap. Use this only when
    /// the input is genuinely arbitrary order.
    pub fn from_unsorted(cohort: CohortId, chrom: Chromosome, mut vcfs: Vec<VariantVcf>) -> Self {
        vcfs.sort_unstable();
        vcfs.dedup();
        Self { cohort, chrom, vcfs }
    }

    pub fn cohort(&self) -> &CohortId {
        &self.cohort
    }
    pub fn chromosome(&self) -> &Chromosome {
        &self.chrom
    }
    pub fn as_slice(&self) -> &[VariantVcf] {
        &self.vcfs
    }
    pub fn len(&self) -> usize {
        self.vcfs.len()
    }
    pub fn is_empty(&self) -> bool {
        self.vcfs.is_empty()
    }
}

/// Per-variant metadata loaded from `variants.parquet`. Same shape as the
/// underlying storage row — the alias keeps the engine API consistent
/// with section 3c of the design doc without duplicating the struct.
pub type VariantMetadata = VariantIndexEntry;

/// Borrowed row yielded by `ChromosomeView::stream`. Holds the dense
/// index plus a borrow of the row metadata.
pub struct VariantRow<'a> {
    pub vcf: VariantVcf,
    pub meta: &'a VariantMetadata,
}

/// Sparse carrier output for a batch of variants. `entries[i]` aligns
/// with the i-th input vcf in the slice that produced this batch.
pub struct CarrierBatch {
    pub entries: Vec<super::variants::CarrierList>,
}
