//! Registry of attached external annotation databases.
//!
//! On-disk shape: `<store_root>/annotations/refs.toml`. Each entry names
//! an annotation source by alias and points at the directory holding the
//! data. The registry never owns or copies the bytes — it is only the
//! resolver from alias → typed reader.
//!
//! Default registry seeds from the user's `Config`: `favor-base` and
//! `favor-full` map to the configured FAVOR root + tier; `favor-tissue`
//! maps to the tissue add-on dir. `attach()` writes a refs.toml entry
//! that overrides or extends the defaults; the on-disk file wins when an
//! alias appears in both.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::config::{Config, Tier};
use crate::error::CohortError;
use crate::store::annotation::{AnnotationDb, TissueDb};
use crate::store::manifest::write_atomic;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum AnnotationKind {
    /// FAVOR base or full tier — chromosome-sharded `sorted.parquet`.
    Favor { tier: Tier },
    /// FAVOR tissue add-on — per-table parquet directories.
    Tissue,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationRef {
    pub name: String,
    #[serde(flatten)]
    pub kind: AnnotationKind,
    pub path: PathBuf,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
struct RefsFile {
    #[serde(default)]
    refs: Vec<AnnotationRef>,
}

pub struct AnnotationRegistry {
    refs_path: PathBuf,
    entries: BTreeMap<String, AnnotationRef>,
}

impl AnnotationRegistry {
    /// Build a registry by merging the on-disk `refs.toml` (if present)
    /// over the defaults derived from the user's `Config`. Entries in
    /// the file shadow defaults of the same alias.
    pub fn load(refs_path: PathBuf, config: &Config) -> Result<Self, CohortError> {
        let mut entries: BTreeMap<String, AnnotationRef> = BTreeMap::new();
        for entry in default_entries(config) {
            entries.insert(entry.name.clone(), entry);
        }
        if refs_path.exists() {
            let body = std::fs::read_to_string(&refs_path).map_err(|e| {
                CohortError::Resource(format!("read {}: {e}", refs_path.display()))
            })?;
            let parsed: RefsFile = toml::from_str(&body).map_err(|e| {
                CohortError::Input(format!("parse {}: {e}", refs_path.display()))
            })?;
            for entry in parsed.refs {
                entries.insert(entry.name.clone(), entry);
            }
        }
        Ok(Self {
            refs_path,
            entries,
        })
    }

    pub fn get(&self, name: &str) -> Option<&AnnotationRef> {
        self.entries.get(name)
    }

    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.entries.keys().map(|s| s.as_str())
    }

    /// Open the FAVOR-style annotation DB at `name`. Errors if the
    /// alias is unknown or points at something that is not a Favor ref.
    pub fn open_db(&self, name: &str) -> Result<AnnotationDb, CohortError> {
        let r = self.require(name)?;
        match r.kind {
            AnnotationKind::Favor { tier } => {
                let parent = r.path.parent().ok_or_else(|| {
                    CohortError::Input(format!(
                        "annotation ref {name}: path {} has no parent",
                        r.path.display()
                    ))
                })?;
                AnnotationDb::open_tier(tier, parent)
            }
            AnnotationKind::Tissue => Err(CohortError::Input(format!(
                "annotation ref {name} is a tissue ref; use open_tissue"
            ))),
        }
    }

    /// Open the tissue add-on DB at `name`.
    pub fn open_tissue(&self, name: &str) -> Result<TissueDb, CohortError> {
        let r = self.require(name)?;
        match r.kind {
            AnnotationKind::Tissue => TissueDb::open(&r.path),
            AnnotationKind::Favor { .. } => Err(CohortError::Input(format!(
                "annotation ref {name} is a favor ref; use open_db"
            ))),
        }
    }

    /// Add or replace an entry, then persist `refs.toml` atomically.
    /// Default entries (those derived from `Config`) are not stored on
    /// disk; only `attach`'d entries appear in the file so a fresh
    /// install picks up new defaults transparently.
    pub fn attach(
        &mut self,
        name: String,
        kind: AnnotationKind,
        path: PathBuf,
    ) -> Result<(), CohortError> {
        let entry = AnnotationRef {
            name: name.clone(),
            kind,
            path,
        };
        self.entries.insert(name, entry);
        self.persist()
    }

    fn persist(&self) -> Result<(), CohortError> {
        let user_refs: Vec<AnnotationRef> = self
            .entries
            .values()
            .filter(|r| !is_default_alias(&r.name))
            .cloned()
            .collect();
        let body = toml::to_string_pretty(&RefsFile { refs: user_refs })
            .map_err(|e| CohortError::Resource(format!("serialize refs.toml: {e}")))?;
        write_atomic(&self.refs_path, body.as_bytes())
    }

    fn require(&self, name: &str) -> Result<&AnnotationRef, CohortError> {
        self.entries.get(name).ok_or_else(|| {
            CohortError::Input(format!(
                "annotation ref {name} not registered; known: {}",
                self.entries
                    .keys()
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ")
            ))
        })
    }
}

fn default_entries(config: &Config) -> Vec<AnnotationRef> {
    let mut out = Vec::new();
    let root = config.root_dir();
    if !root.as_os_str().is_empty() {
        out.push(AnnotationRef {
            name: "favor-base".into(),
            kind: AnnotationKind::Favor { tier: Tier::Base },
            path: root.join(Tier::Base.as_str()),
        });
        out.push(AnnotationRef {
            name: "favor-full".into(),
            kind: AnnotationKind::Favor { tier: Tier::Full },
            path: root.join(Tier::Full.as_str()),
        });
        out.push(AnnotationRef {
            name: "favor-tissue".into(),
            kind: AnnotationKind::Tissue,
            path: config.tissue_dir(),
        });
    }
    out
}

/// Default aliases are stable names the registry always seeds from the
/// `Config`; they should never be persisted into refs.toml because the
/// disk file would shadow user changes to the `Config`.
fn is_default_alias(name: &str) -> bool {
    matches!(name, "favor-base" | "favor-full" | "favor-tissue")
}

/// Convenience: pick the default alias for a tier.
pub fn favor_alias_for(tier: Tier) -> &'static str {
    match tier {
        Tier::Base => "favor-base",
        Tier::Full => "favor-full",
    }
}

/// Convenience: open the FAVOR DB for `tier` from the standard
/// configured root, bypassing the registry. Used by command sites that
/// take a `Tier` directly (`commands/annotate`, `commands/inspect`).
pub fn open_favor_tier(
    config: &Config,
    refs_path: &Path,
    tier: Tier,
) -> Result<AnnotationDb, CohortError> {
    let registry = AnnotationRegistry::load(refs_path.to_path_buf(), config)?;
    registry.open_db(favor_alias_for(tier))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::DataConfig;

    fn test_config(root: &Path) -> Config {
        let mut c = Config::default();
        c.data = DataConfig {
            root_dir: root.to_string_lossy().into_owned(),
            tier: Tier::Base,
            packs: Vec::new(),
        };
        c
    }

    #[test]
    fn defaults_seed_from_config() {
        let dir = tempfile::tempdir().unwrap();
        let cfg = test_config(dir.path());
        let registry = AnnotationRegistry::load(dir.path().join("refs.toml"), &cfg).unwrap();
        assert!(registry.get("favor-base").is_some());
        assert!(registry.get("favor-full").is_some());
        assert!(registry.get("favor-tissue").is_some());
    }

    #[test]
    fn attach_persists_only_user_refs() {
        let dir = tempfile::tempdir().unwrap();
        let cfg = test_config(dir.path());
        let refs_path = dir.path().join("refs.toml");
        let mut registry = AnnotationRegistry::load(refs_path.clone(), &cfg).unwrap();
        registry
            .attach(
                "ukb-v3".into(),
                AnnotationKind::Favor { tier: Tier::Full },
                dir.path().join("ukb_v3"),
            )
            .unwrap();

        let body = std::fs::read_to_string(&refs_path).unwrap();
        assert!(body.contains("ukb-v3"));
        assert!(!body.contains("favor-base"));

        let reloaded = AnnotationRegistry::load(refs_path, &cfg).unwrap();
        assert!(reloaded.get("ukb-v3").is_some());
        assert!(reloaded.get("favor-base").is_some());
    }
}
