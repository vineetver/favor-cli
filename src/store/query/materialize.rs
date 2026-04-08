//! Output shapes for hot-path engine queries: dense dosages over a
//! caller-owned scratch slice. The engine never allocates inside the
//! inner loop — `DenseDosageBatch` borrows from the scratch.

/// Pre-allocated dense scratch for `dense_dosages_into`. Owned by the
/// caller; reused across iterations of an outer loop.
#[derive(Default)]
pub struct DenseDosageScratch {
    data: Vec<f32>,
    n_samples: u32,
    len: usize,
}

impl DenseDosageScratch {
    pub fn new() -> Self {
        Self::default()
    }

    /// Resize the underlying buffer to fit a [n_samples × len] matrix.
    /// Idempotent: calls with the same size do not reallocate.
    pub fn reserve(&mut self, n_samples: u32, len: usize) {
        let needed = n_samples as usize * len;
        if self.data.len() < needed {
            self.data.resize(needed, 0.0);
        }
        self.n_samples = n_samples;
        self.len = len;
    }

    /// Borrow the live region of the buffer (the `[0..n_samples × len]`
    /// prefix). Engine fills this; caller reads it via
    /// `DenseDosageBatch`.
    pub(crate) fn live_mut(&mut self) -> &mut [f32] {
        let n = self.n_samples as usize * self.len;
        &mut self.data[..n]
    }

    pub fn n_samples(&self) -> u32 {
        self.n_samples
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
}

/// Borrowed dense dosage matrix returned by
/// `ChromosomeView::dense_dosages_into`. Row-major: index
/// `(sample_idx, variant_local) -> data[sample_idx * len + variant_local]`.
pub struct DenseDosageBatch<'a> {
    pub n_samples: u32,
    pub len: usize,
    pub data: &'a mut [f32],
}
