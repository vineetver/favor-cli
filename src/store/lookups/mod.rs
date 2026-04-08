//! Extension point for derived per-cohort indexes.
//!
//! Each lookup owns its own on-disk format under
//! `<store_root>/cache/lookups/<cohort_id>/<name>/`. The framework's job
//! is small: check whether the index file exists, dispatch to the
//! impl's `build()` if missing, then call `query()`. The on-disk
//! invalidation story is the same as score caches: rebuilding a cohort
//! changes its content fingerprint and `CacheStore::prune_cohort` drops
//! the lookup directory along with the score cache.
//!
//! Adding a new lookup is one new file in `src/store/lookups/` plus a
//! `Store::lookup::<L>(...)` call site. No core type touches.

pub mod rsid;
pub mod sample_carriers;

use std::path::{Path, PathBuf};

use crate::error::CohortError;
use crate::store::cohort::CohortHandle;

pub use rsid::RsidLookup;
pub use sample_carriers::SampleCarriersLookup;

/// One on-disk derived index. Implementors own their own file format
/// under the directory the framework hands them. Lifetime is bound by
/// the file path: rebuilding a cohort moves the directory and the next
/// query rebuilds.
pub trait Lookup: Send + Sync + 'static {
    /// Stable name. Used as the on-disk directory stem.
    const NAME: &'static str;

    /// Caller-supplied input key (e.g. an rsid string, a sample index).
    type Key: ?Sized;

    /// Engine-supplied output value (e.g. matching vcfs).
    type Value;

    /// Build the index from the cohort. Called once when the index
    /// directory does not exist or has been invalidated. Implementors
    /// should write atomic files inside `index_dir` and rely on the
    /// framework to mkdir the parent.
    fn build(
        &self,
        cohort: &CohortHandle<'_>,
        index_dir: &Path,
    ) -> Result<(), CohortError>;

    /// Query the on-disk index. Called on every lookup after the index
    /// is materialized. Should mmap files inside `index_dir`.
    fn query(
        &self,
        cohort: &CohortHandle<'_>,
        index_dir: &Path,
        key: &Self::Key,
    ) -> Result<Self::Value, CohortError>;

    /// Sentinel file inside `index_dir` whose presence means the build
    /// completed. The framework checks this; impls must write it last
    /// (atomically) so a crash mid-build doesn't leave a half-baked
    /// index visible to the next query.
    fn ready_marker() -> &'static str {
        "READY"
    }
}

/// Resolve the on-disk directory for a lookup, mkdir its parent, and
/// build the index if no `READY` marker is present yet.
pub(crate) fn ensure_built<L: Lookup>(
    lookup: &L,
    cohort: &CohortHandle<'_>,
) -> Result<PathBuf, CohortError> {
    let store = cohort.store();
    let dir = store
        .layout()
        .lookup_index_dir(cohort.id(), L::NAME);
    if !dir.join(L::ready_marker()).exists() {
        std::fs::create_dir_all(&dir).map_err(|e| {
            CohortError::Resource(format!("create {}: {e}", dir.display()))
        })?;
        lookup.build(cohort, &dir)?;
        // The build is responsible for writing the ready marker. The
        // framework only verifies it landed.
        if !dir.join(L::ready_marker()).exists() {
            return Err(CohortError::Resource(format!(
                "lookup {} did not write {}",
                L::NAME,
                L::ready_marker()
            )));
        }
    }
    Ok(dir)
}
