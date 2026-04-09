//! Derived caches keyed by cohort id. Today: per-phenotype score cache.

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
    /// `<store_root>/cache/score_cache/<cohort>/<key>/`.
    pub fn score_cache_dir(&self, cohort: &CohortId, key: &CacheKey) -> PathBuf {
        self.layout.score_cache_dir(cohort, key)
    }

    /// Delete the score-cache directory owned by `cohort`. Called when
    /// `--rebuild-store` flips the content fingerprint and the existing
    /// cache no longer matches.
    pub fn prune_cohort(&self, cohort: &CohortId) -> Result<PruneSummary, CohortError> {
        let mut summary = PruneSummary::default();
        for sub in ["score_cache", "lookups"] {
            let dir = self.layout.cache_root().join(sub).join(cohort.as_str());
            if dir.is_dir() {
                summary.bytes_freed += dir_size(&dir);
                std::fs::remove_dir_all(&dir).map_err(|e| {
                    CohortError::Resource(format!("rm {}: {e}", dir.display()))
                })?;
                if sub == "score_cache" {
                    summary.removed_score_caches += 1;
                } else {
                    summary.removed_lookup_indexes += 1;
                }
            }
        }
        Ok(summary)
    }

    /// Walk every per-cohort cache subdirectory and remove the ones whose
    /// cohort id is not in `live`. Used by `cohort store gc`.
    pub fn prune_orphans(&self, live: &[CohortId]) -> Result<PruneSummary, CohortError> {
        let mut summary = PruneSummary::default();
        let alive: std::collections::HashSet<&str> =
            live.iter().map(|c| c.as_str()).collect();
        for sub in ["score_cache", "lookups"] {
            let parent = self.layout.cache_root().join(sub);
            if !parent.is_dir() {
                continue;
            }
            for entry in std::fs::read_dir(&parent).map_err(|e| {
                CohortError::Resource(format!("read {}: {e}", parent.display()))
            })? {
                let entry = entry.map_err(|e| {
                    CohortError::Resource(format!("read {}: {e}", parent.display()))
                })?;
                let p = entry.path();
                if !p.is_dir() {
                    continue;
                }
                let name = match p.file_name().and_then(|n| n.to_str()) {
                    Some(n) => n,
                    None => continue,
                };
                if alive.contains(name) {
                    continue;
                }
                summary.bytes_freed += dir_size(&p);
                std::fs::remove_dir_all(&p).map_err(|e| {
                    CohortError::Resource(format!("rm {}: {e}", p.display()))
                })?;
                if sub == "score_cache" {
                    summary.removed_score_caches += 1;
                } else {
                    summary.removed_lookup_indexes += 1;
                }
            }
        }
        Ok(summary)
    }
}

#[derive(Debug, Default, serde::Serialize)]
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
    fn prune_orphans_removes_only_dead_cohorts() {
        let tmp = tempfile::tempdir().unwrap();
        let layout = Layout::new(tmp.path().to_path_buf());
        let score = layout.cache_root().join("score_cache");
        std::fs::create_dir_all(score.join("alive")).unwrap();
        std::fs::create_dir_all(score.join("dead")).unwrap();
        std::fs::write(score.join("alive/u.bin"), b"x").unwrap();
        std::fs::write(score.join("dead/u.bin"), b"yy").unwrap();

        let cache = CacheStore::new(&layout);
        let summary = cache.prune_orphans(&[CohortId::new("alive")]).unwrap();

        assert!(score.join("alive/u.bin").exists());
        assert!(!score.join("dead").exists());
        assert_eq!(summary.removed_score_caches, 1);
        assert_eq!(summary.removed_lookup_indexes, 0);
        assert!(summary.bytes_freed >= 2);
    }
}
