//! Closed selector enum for variant queries against a `ChromosomeView`.
//!
//! Adding a new arm here is intentionally a touch on this file plus
//! `ChromosomeView::select` — the closed set is small and audited. For
//! open-ended patterns (rsid, custom annotations, third-party indexes)
//! use the `Lookup` extension mechanism in `src/store/lookups/`.
//!
//! `select(...)` returns a `SortedVcfs` scoped to one (cohort, chromosome).
//! For multi-chromosome selectors (`Gene`, `List`), the caller iterates
//! `cohort.chromosomes()` and calls `select` once per chromosome.

pub mod materialize;

use crate::store::cohort::types::VariantVcf;
use crate::store::ids::ListId;
use crate::types::{Chromosome, Consequence};

/// Inclusive 1-based genomic interval scoped to one chromosome.
#[derive(Debug, Clone)]
pub struct Region {
    pub chromosome: Chromosome,
    pub start: u32,
    pub end: u32,
}

/// Closed enum of compile-time selector patterns. The chromosome view
/// resolves each arm against its in-memory index. Arms that name a
/// list (`List`) require the variant-list subsystem; arms that need
/// per-position metadata fall through to the `VariantIndex` walker.
#[derive(Debug, Clone)]
pub enum VariantSelector {
    /// Every variant in the chromosome. O(n).
    All,

    /// All variants tagged with the gene name. O(1) lookup +
    /// O(gene_size) to build the result.
    Gene(String),

    /// Range overlap. O(n) without an interval index, O(log n + matches)
    /// with one (Phase 7).
    Region(Region),

    /// One variant by its dense per-cohort index. O(1).
    Vcf(VariantVcf),

    /// Many variants by index. Pre-sorted by caller. O(len).
    Vcfs(Vec<VariantVcf>),

    /// One variant by `vid` ("chr-pos-ref-alt"). O(1) hash lookup.
    Vid(String),

    /// Many variants by vid. O(len) hash lookups.
    Vids(Vec<String>),

    /// All variants in a previously-imported list. Resolves the list
    /// against the cohort universe (vid join). O(list_size). Returns
    /// `CohortError::DataMissing` until the list subsystem lands
    /// (Phase 5).
    List(ListId),

    /// MAF threshold. Closed-form metadata filter. O(n).
    MafBelow(f64),

    /// Specific consequence (missense, plof, ...). O(n).
    Consequence(Consequence),

    /// Boolean intersect. Resolves both arms then sorted-merge intersect.
    And(Box<VariantSelector>, Box<VariantSelector>),

    /// Boolean union. Sorted-merge union. O(|left| + |right|).
    Or(Box<VariantSelector>, Box<VariantSelector>),
}

/// Sorted-merge intersection of two ascending vcf vectors.
pub(crate) fn intersect_sorted(a: &[VariantVcf], b: &[VariantVcf]) -> Vec<VariantVcf> {
    let mut out = Vec::with_capacity(a.len().min(b.len()));
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

/// Sorted-merge union of two ascending vcf vectors.
pub(crate) fn union_sorted(a: &[VariantVcf], b: &[VariantVcf]) -> Vec<VariantVcf> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => {
                out.push(a[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(b[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&a[i..]);
    out.extend_from_slice(&b[j..]);
    out
}

/// Bridge from a STAAR-style mask predicate to the closest selector.
/// Predicates that test only `consequence` collapse cleanly into
/// `Consequence(c)`; anything else stays as a post-filter applied by
/// the caller after `select` returns.
#[allow(dead_code)]
pub fn bridge_consequence(c: Consequence) -> VariantSelector {
    VariantSelector::Consequence(c)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(xs: &[u32]) -> Vec<VariantVcf> {
        xs.iter().copied().map(VariantVcf).collect()
    }

    #[test]
    fn intersect_keeps_only_common_indices() {
        assert_eq!(intersect_sorted(&v(&[1, 3, 5]), &v(&[3, 4, 5, 6])), v(&[3, 5]));
        assert_eq!(intersect_sorted(&v(&[]), &v(&[1, 2])), v(&[]));
        assert_eq!(intersect_sorted(&v(&[1, 2, 3]), &v(&[1, 2, 3])), v(&[1, 2, 3]));
    }

    #[test]
    fn union_dedups_overlap() {
        assert_eq!(
            union_sorted(&v(&[1, 3, 5]), &v(&[2, 3, 6])),
            v(&[1, 2, 3, 5, 6])
        );
        assert_eq!(union_sorted(&v(&[]), &v(&[1, 2])), v(&[1, 2]));
        assert_eq!(union_sorted(&v(&[1, 2]), &v(&[])), v(&[1, 2]));
    }
}
