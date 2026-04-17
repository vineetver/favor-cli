//! Path conventions for the store. Every other module routes path
//! computation through `Layout` so the engine has one source of truth
//! for naming.

use std::path::{Path, PathBuf};

use super::ids::{CacheKey, CohortId};

pub struct Layout {
    root: PathBuf,
}

impl Layout {
    pub fn new(root: PathBuf) -> Self {
        Self { root }
    }

    pub fn root(&self) -> &Path {
        &self.root
    }

    pub fn cohorts_root(&self) -> PathBuf {
        self.root.join("cohorts")
    }

    pub fn cohort_dir(&self, id: &CohortId) -> PathBuf {
        self.cohorts_root().join(id.as_str())
    }

    pub fn lists_root(&self) -> PathBuf {
        self.root.join("lists")
    }

    pub fn annotations_refs(&self) -> PathBuf {
        self.root.join("annotations").join("refs.toml")
    }

    pub fn cache_root(&self) -> PathBuf {
        self.root.join("cache")
    }

    pub fn score_cache_dir(&self, cohort: &CohortId, key: &CacheKey) -> PathBuf {
        self.cache_root()
            .join("score_cache")
            .join(cohort.as_str())
            .join(key.as_str())
    }

    pub fn null_model_cache_dir(&self, cohort: &CohortId, key: &CacheKey) -> PathBuf {
        self.cache_root()
            .join("null_model")
            .join(cohort.as_str())
            .join(key.as_str())
    }

    pub fn grm_cache_dir(&self, cohort: &CohortId, key: &CacheKey) -> PathBuf {
        self.cache_root()
            .join("grm")
            .join(cohort.as_str())
            .join(key.as_str())
    }

    /// Subdirectories that `Store::open` materializes lazily.
    pub(super) fn known_subdirs(&self) -> [PathBuf; 4] {
        [
            self.cohorts_root(),
            self.lists_root(),
            self.cache_root(),
            self.root.join("annotations"),
        ]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paths_are_anchored_at_root() {
        let layout = Layout::new(PathBuf::from("/srv/.cohort"));
        assert_eq!(layout.cohorts_root(), PathBuf::from("/srv/.cohort/cohorts"));
        assert_eq!(layout.lists_root(), PathBuf::from("/srv/.cohort/lists"));
        assert_eq!(layout.cache_root(), PathBuf::from("/srv/.cohort/cache"));
        assert_eq!(
            layout.annotations_refs(),
            PathBuf::from("/srv/.cohort/annotations/refs.toml")
        );
    }

    #[test]
    fn cache_paths_isolate_by_cohort() {
        let layout = Layout::new(PathBuf::from("/r"));
        let id = CohortId::new("c1");
        let key = CacheKey::new("k1");
        assert_eq!(
            layout.score_cache_dir(&id, &key),
            PathBuf::from("/r/cache/score_cache/c1/k1")
        );
    }
}
