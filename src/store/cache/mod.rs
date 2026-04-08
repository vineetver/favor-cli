//! Derived caches keyed by cohort id. Today: per-phenotype score cache.
//! Phase 7 lookups will share the `prune_orphans` GC path.

pub mod score_cache;

use std::path::PathBuf;

use crate::error::CohortError;
use crate::store::ids::{CacheKey, CohortId};
use crate::store::layout::Layout;

pub struct CacheStore<'a> {
    layout: &'a Layout,
}

impl<'a> CacheStore<'a> {
    pub(crate) fn new(layout: &'a Layout) -> Self {
        Self { layout }
    }

    /// Directory for one (cohort, cache_key) pair. The on-disk path is
    /// `<store_root>/cache/score_cache/<cohort>/<key>/`. Replaces the
    /// pre-engine `score_cache::cache_dir(store_dir, key)` which
    /// entangled the cache with the cohort dir.
    pub fn score_cache_dir(&self, cohort: &CohortId, key: &CacheKey) -> PathBuf {
        self.layout.score_cache_dir(cohort, key)
    }

    /// Delete the score-cache and lookup-index directories owned by
    /// `cohort`. Called when `--rebuild-store` flips the content
    /// fingerprint and the existing caches no longer match.
    pub fn prune_cohort(&self, cohort: &CohortId) -> Result<PruneSummary, CohortError> {
        let mut summary = PruneSummary::default();
        let score_dir = self
            .layout
            .cache_root()
            .join("score_cache")
            .join(cohort.as_str());
        if score_dir.is_dir() {
            summary.bytes_freed += dir_size(&score_dir);
            std::fs::remove_dir_all(&score_dir).map_err(|e| {
                CohortError::Resource(format!("rm {}: {e}", score_dir.display()))
            })?;
            summary.removed_score_caches = 1;
        }
        let lookup_dir = self
            .layout
            .cache_root()
            .join("lookups")
            .join(cohort.as_str());
        if lookup_dir.is_dir() {
            summary.bytes_freed += dir_size(&lookup_dir);
            std::fs::remove_dir_all(&lookup_dir).map_err(|e| {
                CohortError::Resource(format!("rm {}: {e}", lookup_dir.display()))
            })?;
            summary.removed_lookup_indexes = 1;
        }
        Ok(summary)
    }

    /// Delete every score-cache and lookup-index directory whose cohort
    /// id is not in `live`. Called from the future `cohort store gc`
    /// command.
    pub fn prune_orphans(&self, live: &[CohortId]) -> Result<PruneSummary, CohortError> {
        let mut summary = PruneSummary::default();
        let live_set: std::collections::HashSet<&str> =
            live.iter().map(|c| c.as_str()).collect();

        let score_root = self.layout.cache_root().join("score_cache");
        if score_root.is_dir() {
            for entry in std::fs::read_dir(&score_root).map_err(|e| {
                CohortError::Resource(format!("read {}: {e}", score_root.display()))
            })? {
                let entry = entry
                    .map_err(|e| CohortError::Resource(format!("read entry: {e}")))?;
                let path = entry.path();
                let name = match path.file_name().and_then(|s| s.to_str()) {
                    Some(n) => n,
                    None => continue,
                };
                if live_set.contains(name) {
                    continue;
                }
                summary.bytes_freed += dir_size(&path);
                std::fs::remove_dir_all(&path).map_err(|e| {
                    CohortError::Resource(format!("rm {}: {e}", path.display()))
                })?;
                summary.removed_score_caches += 1;
            }
        }

        let lookups_root = self.layout.cache_root().join("lookups");
        if lookups_root.is_dir() {
            for entry in std::fs::read_dir(&lookups_root).map_err(|e| {
                CohortError::Resource(format!("read {}: {e}", lookups_root.display()))
            })? {
                let entry = entry
                    .map_err(|e| CohortError::Resource(format!("read entry: {e}")))?;
                let path = entry.path();
                let name = match path.file_name().and_then(|s| s.to_str()) {
                    Some(n) => n,
                    None => continue,
                };
                if live_set.contains(name) {
                    continue;
                }
                summary.bytes_freed += dir_size(&path);
                std::fs::remove_dir_all(&path).map_err(|e| {
                    CohortError::Resource(format!("rm {}: {e}", path.display()))
                })?;
                summary.removed_lookup_indexes += 1;
            }
        }

        Ok(summary)
    }
}

#[derive(Debug, Default)]
pub struct PruneSummary {
    pub removed_score_caches: usize,
    pub removed_lookup_indexes: usize,
    pub bytes_freed: u64,
}

fn dir_size(path: &std::path::Path) -> u64 {
    let mut total = 0u64;
    if let Ok(entries) = std::fs::read_dir(path) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.is_file() {
                total += std::fs::metadata(&p).map(|m| m.len()).unwrap_or(0);
            } else if p.is_dir() {
                total += dir_size(&p);
            }
        }
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::store::layout::Layout;

    #[test]
    fn prune_orphans_keeps_live_cohorts() {
        let dir = tempfile::tempdir().unwrap();
        let layout = Layout::new(dir.path().to_path_buf());
        let cache = CacheStore::new(&layout);

        let live_id = CohortId::new("live_cohort");
        let dead_id = CohortId::new("dead_cohort");
        let key = CacheKey::new("k1");

        let live_dir = layout.score_cache_dir(&live_id, &key);
        let dead_dir = layout.score_cache_dir(&dead_id, &key);
        std::fs::create_dir_all(&live_dir).unwrap();
        std::fs::create_dir_all(&dead_dir).unwrap();
        std::fs::write(live_dir.join("scores.bin"), [0u8; 64]).unwrap();
        std::fs::write(dead_dir.join("scores.bin"), [0u8; 128]).unwrap();

        let summary = cache.prune_orphans(&[live_id.clone()]).unwrap();
        assert_eq!(summary.removed_score_caches, 1);
        assert_eq!(summary.removed_lookup_indexes, 0);
        assert_eq!(summary.bytes_freed, 128);

        assert!(live_dir.exists());
        assert!(!dead_dir.parent().unwrap().exists());
    }
}
