//! Path conventions for the store. Every other module routes path
//! computation through `Layout` so the engine has one source of truth
//! for naming.

use std::path::{Path, PathBuf};

use crate::types::Chromosome;

use super::ids::{CacheKey, CohortId, ListId};

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

    pub fn store_toml(&self) -> PathBuf {
        self.root.join("store.toml")
    }

    pub fn cohorts_root(&self) -> PathBuf {
        self.root.join("cohorts")
    }

    pub fn cohort_dir(&self, id: &CohortId) -> PathBuf {
        self.cohorts_root().join(id.as_str())
    }

    pub fn cohort_chromosome_dir(&self, id: &CohortId, chrom: &Chromosome) -> PathBuf {
        self.cohort_dir(id).join(format!("chromosome={}", chrom.label()))
    }

    pub fn lists_root(&self) -> PathBuf {
        self.root.join("lists")
    }

    pub fn list_dir(&self, id: &ListId) -> PathBuf {
        self.lists_root().join(id.as_str())
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

    pub fn lookup_index_dir(&self, cohort: &CohortId, lookup_name: &str) -> PathBuf {
        self.cache_root()
            .join("lookups")
            .join(cohort.as_str())
            .join(lookup_name)
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
        assert_eq!(layout.store_toml(), PathBuf::from("/srv/.cohort/store.toml"));
        assert_eq!(layout.cohorts_root(), PathBuf::from("/srv/.cohort/cohorts"));
        assert_eq!(layout.lists_root(), PathBuf::from("/srv/.cohort/lists"));
        assert_eq!(layout.cache_root(), PathBuf::from("/srv/.cohort/cache"));
        assert_eq!(
            layout.annotations_refs(),
            PathBuf::from("/srv/.cohort/annotations/refs.toml")
        );
    }

    #[test]
    fn cohort_chromosome_dir_uses_chrom_label() {
        let layout = Layout::new(PathBuf::from("/r"));
        let id = CohortId::new("ukb_exome");
        let p = layout.cohort_chromosome_dir(&id, &Chromosome::Autosome(22));
        assert_eq!(p, PathBuf::from("/r/cohorts/ukb_exome/chromosome=22"));
        let p = layout.cohort_chromosome_dir(&id, &Chromosome::X);
        assert_eq!(p, PathBuf::from("/r/cohorts/ukb_exome/chromosome=X"));
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
        assert_eq!(
            layout.lookup_index_dir(&id, "rsid"),
            PathBuf::from("/r/cache/lookups/c1/rsid")
        );
    }
}
