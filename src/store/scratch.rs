//! Pre-reserved scratch buffer pool. Hot-path methods (carriers, dense
//! dosages) borrow buffers from the pool instead of allocating per call.
//! Bodies are stubs until the query layer wires them in.

pub struct ScratchPool {
    #[allow(dead_code)]
    budget_bytes: usize,
}

impl ScratchPool {
    pub fn new(budget_bytes: usize) -> Self {
        Self { budget_bytes }
    }

    pub fn carrier_scratch(&self, _n_carriers_hint: usize) -> CarrierScratch<'_> {
        CarrierScratch { _marker: std::marker::PhantomData }
    }

    pub fn dense_scratch(&self, _n_samples: u32, _len: usize) -> DenseDosageScratch {
        DenseDosageScratch::new()
    }
}

pub struct CarrierScratch<'a> {
    _marker: std::marker::PhantomData<&'a ()>,
}

pub struct DenseDosageScratch {
    _data: Vec<f32>,
}

impl DenseDosageScratch {
    pub fn new() -> Self {
        Self { _data: Vec::new() }
    }

    pub fn reserve(&mut self, _n_samples: u32, _len: usize) {}
}

impl Default for DenseDosageScratch {
    fn default() -> Self {
        Self::new()
    }
}
