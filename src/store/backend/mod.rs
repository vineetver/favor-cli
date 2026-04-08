//! Filesystem abstraction. Returns parquet readers and byte slices, not
//! typed variants — subsystems do the typed view. One impl today
//! (`LocalFs`); the trait exists so a remote backend can plug in later
//! without rewriting subsystems.

use std::ops::Deref;
use std::path::Path;

use arrow::record_batch::RecordBatchReader;
use memmap2::Mmap;

use crate::error::CohortError;

pub mod local;

pub use local::LocalFs;

/// Boxed `RecordBatchReader` returned by `Backend::open_parquet`. The
/// trait object lets the caller iterate batches without leaking the
/// concrete parquet reader type into subsystem signatures.
pub type BoxedBatchReader = Box<dyn RecordBatchReader + Send>;

/// Owned mmap region. Derefs to `&[u8]` so subsystems can carve typed
/// views without knowing what backend produced it.
pub struct MappedBytes {
    inner: Mmap,
}

impl MappedBytes {
    pub(crate) fn new(inner: Mmap) -> Self {
        Self { inner }
    }
}

impl Deref for MappedBytes {
    type Target = [u8];
    fn deref(&self) -> &[u8] {
        &self.inner
    }
}

impl AsRef<[u8]> for MappedBytes {
    fn as_ref(&self) -> &[u8] {
        &self.inner
    }
}

pub(crate) trait Backend: Send + Sync {
    /// Open a parquet file as a streaming RecordBatch reader.
    fn open_parquet(&self, path: &Path) -> Result<BoxedBatchReader, CohortError>;

    /// Memory-map a binary file (sparse_g.bin, lookup indexes).
    fn mmap(&self, path: &Path) -> Result<MappedBytes, CohortError>;

    /// Atomic write: write to tmp, fsync, rename, fsync parent.
    fn write_atomic(&self, path: &Path, bytes: &[u8]) -> Result<(), CohortError>;

    /// Replace a directory atomically by rename.
    fn swap_dir(&self, staging: &Path, final_path: &Path) -> Result<(), CohortError>;
}
