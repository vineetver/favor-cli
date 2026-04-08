//! Pre-reserved scratch buffer pool. Hot-path methods (carriers, dense
//! dosages) reuse buffers from the pool instead of allocating per call.
//!
//! Today the pool is a soft cap on total bytes plus a per-call factory
//! that hands back an owned `DenseDosageScratch`. Reuse across calls is
//! the caller's job: hold one scratch across the gene loop and the
//! engine fills it on every iteration.

use crate::store::query::materialize::DenseDosageScratch;

pub struct ScratchPool {
    #[allow(dead_code)]
    budget_bytes: usize,
}

impl ScratchPool {
    pub fn new(budget_bytes: usize) -> Self {
        Self { budget_bytes }
    }

    /// Reserve a fresh dense dosage scratch sized for the given matrix.
    /// The caller owns the returned scratch and reuses it across the
    /// outer loop — the pool exists today only as a place to hang
    /// future cross-call sharing.
    pub fn dense_scratch(&self, n_samples: u32, len: usize) -> DenseDosageScratch {
        let mut s = DenseDosageScratch::new();
        s.reserve(n_samples, len);
        s
    }

    /// Reserve a carrier-list scratch. Today this is just a `Vec`
    /// pre-allocated to the caller's hint; the engine still allocates
    /// the per-variant `CarrierList::entries` Vecs because the carrier
    /// reader produces them. Phase 7 tightens this when the
    /// sample-side index lands.
    pub fn carrier_scratch(&self, n_carriers_hint: usize) -> CarrierScratch {
        CarrierScratch {
            buf: Vec::with_capacity(n_carriers_hint),
        }
    }
}

pub struct CarrierScratch {
    #[allow(dead_code)]
    buf: Vec<crate::store::cohort::variants::CarrierList>,
}
