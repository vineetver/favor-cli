#[allow(dead_code)]
pub mod backend;
pub mod annotation;
pub mod cache;
pub mod cohort;
#[allow(dead_code)]
pub mod config;
#[allow(dead_code)]
pub mod ids;
#[allow(dead_code)]
pub mod layout;
pub mod list;
pub mod manifest;
#[allow(dead_code)]
pub mod scratch;

use std::fs;

use crate::error::CohortError;

use backend::{Backend, LocalFs};
use config::StoreConfig;
use layout::Layout;
use scratch::ScratchPool;

#[allow(dead_code)]
const DEFAULT_SCRATCH_BUDGET_BYTES: usize = 256 * 1024 * 1024;

#[allow(dead_code)]
pub struct Store {
    layout: Layout,
    backend: Box<dyn Backend>,
    scratch: ScratchPool,
}

#[allow(dead_code)]
impl Store {
    pub fn open(config: StoreConfig) -> Result<Self, CohortError> {
        let layout = Layout::new(config.root);
        fs::create_dir_all(layout.root()).map_err(|e| {
            CohortError::Resource(format!(
                "create store root {}: {e}",
                layout.root().display()
            ))
        })?;
        for sub in layout.known_subdirs() {
            fs::create_dir_all(&sub).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", sub.display()))
            })?;
        }
        Ok(Self {
            layout,
            backend: Box::new(LocalFs::new()),
            scratch: ScratchPool::new(DEFAULT_SCRATCH_BUDGET_BYTES),
        })
    }

    pub fn layout(&self) -> &Layout {
        &self.layout
    }

    pub(crate) fn backend(&self) -> &dyn Backend {
        &*self.backend
    }

    #[allow(dead_code)]
    pub fn scratch(&self) -> &ScratchPool {
        &self.scratch
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Read;

    #[test]
    fn open_creates_root_and_subdirs() {
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path().join("store_root");
        let store = Store::open(StoreConfig { root: root.clone() }).unwrap();
        assert!(root.is_dir());
        assert!(store.layout().cohorts_root().is_dir());
        assert!(store.layout().lists_root().is_dir());
        assert!(store.layout().cache_root().is_dir());
    }

    #[test]
    fn backend_write_atomic_round_trips() {
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path().join("store_root");
        let store = Store::open(StoreConfig { root: root.clone() }).unwrap();

        let sentinel = root.join("sentinel.bin");
        store
            .backend()
            .write_atomic(&sentinel, b"phase-1-sentinel")
            .unwrap();

        let reopened = Store::open(StoreConfig { root: root.clone() }).unwrap();
        let mut bytes = Vec::new();
        std::fs::File::open(&sentinel)
            .unwrap()
            .read_to_end(&mut bytes)
            .unwrap();
        assert_eq!(bytes, b"phase-1-sentinel");
        let _ = reopened;
    }

    #[test]
    fn backend_open_parquet_reads_testdata() {
        let path = std::path::Path::new(
            "testdata/ukb_chr22/ingested/chromosome=22/data.parquet",
        );
        if !path.exists() {
            eprintln!("skipping: {} not present", path.display());
            return;
        }
        let tmp = tempfile::tempdir().unwrap();
        let store = Store::open(StoreConfig { root: tmp.path().join("s") }).unwrap();
        let mut reader = store.backend().open_parquet(path).unwrap();
        let batch = reader
            .next()
            .expect("at least one batch")
            .expect("decode ok");
        assert!(batch.num_rows() > 0);
        assert!(batch.num_columns() > 0);
    }
}
