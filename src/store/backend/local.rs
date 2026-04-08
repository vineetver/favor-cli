//! Local filesystem backend.

use std::fs::{self, File};
use std::path::Path;

use memmap2::MmapOptions;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use crate::error::CohortError;
use crate::store::manifest::{fsync_parent, write_atomic};

use super::{Backend, BoxedBatchReader, MappedBytes};

pub struct LocalFs;

impl LocalFs {
    pub fn new() -> Self {
        Self
    }
}

impl Default for LocalFs {
    fn default() -> Self {
        Self::new()
    }
}

impl Backend for LocalFs {
    fn open_parquet(&self, path: &Path) -> Result<BoxedBatchReader, CohortError> {
        let file = File::open(path)
            .map_err(|e| CohortError::Resource(format!("open {}: {e}", path.display())))?;
        let reader = ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| {
                CohortError::Resource(format!("parquet header {}: {e}", path.display()))
            })?
            .build()
            .map_err(|e| {
                CohortError::Resource(format!("parquet reader {}: {e}", path.display()))
            })?;
        Ok(Box::new(reader))
    }

    fn mmap(&self, path: &Path) -> Result<MappedBytes, CohortError> {
        let file = File::open(path)
            .map_err(|e| CohortError::Resource(format!("open {}: {e}", path.display())))?;
        let mmap = unsafe { MmapOptions::new().map(&file) }
            .map_err(|e| CohortError::Resource(format!("mmap {}: {e}", path.display())))?;
        Ok(MappedBytes::new(mmap))
    }

    fn write_atomic(&self, path: &Path, bytes: &[u8]) -> Result<(), CohortError> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!("mkdir {}: {e}", parent.display()))
            })?;
        }
        write_atomic(path, bytes)
    }

    fn swap_dir(&self, staging: &Path, final_path: &Path) -> Result<(), CohortError> {
        if final_path.exists() {
            fs::remove_dir_all(final_path).map_err(|e| {
                CohortError::Resource(format!("rm {}: {e}", final_path.display()))
            })?;
        } else if let Some(parent) = final_path.parent() {
            fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!("mkdir {}: {e}", parent.display()))
            })?;
        }
        fs::rename(staging, final_path).map_err(|e| {
            CohortError::Resource(format!(
                "rename {} -> {}: {e}",
                staging.display(),
                final_path.display()
            ))
        })?;
        fsync_parent(final_path);
        Ok(())
    }
}
