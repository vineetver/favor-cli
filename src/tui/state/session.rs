use std::fs;
use std::io;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::cli::GenomeBuild;
use crate::config::Tier;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum BuildTag {
    Hg38,
    Hg19,
}

impl BuildTag {
    pub fn from_build(b: &GenomeBuild) -> Self {
        match b {
            GenomeBuild::Hg38 => BuildTag::Hg38,
            GenomeBuild::Hg19 => BuildTag::Hg19,
        }
    }
    #[allow(dead_code)]
    pub fn to_build(self) -> GenomeBuild {
        match self {
            BuildTag::Hg38 => GenomeBuild::Hg38,
            BuildTag::Hg19 => GenomeBuild::Hg19,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SessionState {
    pub cwd: PathBuf,
    pub last_artifact: Option<PathBuf>,
    pub transform: Option<TransformSnapshot>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum TransformSnapshot {
    Ingest {
        inputs: Vec<PathBuf>,
        output: Option<PathBuf>,
        emit_sql: bool,
        build: Option<BuildTag>,
    },
    Annotate {
        input: PathBuf,
        output: Option<PathBuf>,
        tier: Tier,
        data_root: PathBuf,
    },
}

#[derive(Debug, Clone)]
pub struct SessionId(String);

impl SessionId {
    pub fn for_cwd(cwd: &Path) -> Self {
        let mut h = Sha256::new();
        h.update(cwd.as_os_str().as_encoded_bytes());
        let digest = h.finalize();
        let mut s = String::with_capacity(16);
        for byte in &digest[..8] {
            s.push_str(&format!("{byte:02x}"));
        }
        SessionId(s)
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }
}

pub struct SessionStore {
    root: PathBuf,
}

impl SessionStore {
    pub fn from_home() -> Option<Self> {
        let home = std::env::var_os("HOME").map(PathBuf::from)?;
        Some(Self {
            root: home.join(".cohort").join("sessions"),
        })
    }

    #[cfg(test)]
    pub fn with_root(root: PathBuf) -> Self {
        Self { root }
    }

    fn dir_for(&self, id: &SessionId) -> PathBuf {
        self.root.join(id.as_str())
    }

    fn path_for(&self, id: &SessionId) -> PathBuf {
        self.dir_for(id).join("state.json")
    }

    pub fn load(&self, id: &SessionId) -> io::Result<Option<SessionState>> {
        let path = self.path_for(id);
        match fs::read(&path) {
            Ok(bytes) => {
                let state: SessionState = serde_json::from_slice(&bytes)
                    .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                Ok(Some(state))
            }
            Err(e) if e.kind() == io::ErrorKind::NotFound => Ok(None),
            Err(e) => Err(e),
        }
    }

    pub fn save(&self, id: &SessionId, state: &SessionState) -> io::Result<()> {
        let dir = self.dir_for(id);
        fs::create_dir_all(&dir)?;
        let final_path = dir.join("state.json");
        let tmp_path = dir.join("state.json.tmp");
        let bytes = serde_json::to_vec_pretty(state)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        fs::write(&tmp_path, &bytes)?;
        fs::rename(&tmp_path, &final_path)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn session_id_stable_for_cwd() {
        let a = SessionId::for_cwd(Path::new("/tmp/foo"));
        let b = SessionId::for_cwd(Path::new("/tmp/foo"));
        let c = SessionId::for_cwd(Path::new("/tmp/bar"));
        assert_eq!(a.as_str(), b.as_str());
        assert_ne!(a.as_str(), c.as_str());
        assert_eq!(a.as_str().len(), 16);
    }

    #[test]
    fn round_trip_save_load() {
        let tmp = std::env::temp_dir().join(format!("cohort-session-test-{}", std::process::id()));
        let store = SessionStore::with_root(tmp.clone());
        let id = SessionId::for_cwd(Path::new("/some/cwd"));
        let state = SessionState {
            cwd: PathBuf::from("/some/cwd"),
            last_artifact: Some(PathBuf::from("/some/cwd/x.parquet")),
            transform: Some(TransformSnapshot::Ingest {
                inputs: vec![PathBuf::from("a.vcf.gz")],
                output: None,
                emit_sql: true,
                build: Some(BuildTag::Hg38),
            }),
        };
        store.save(&id, &state).unwrap();
        let loaded = store.load(&id).unwrap().unwrap();
        assert_eq!(loaded.cwd, state.cwd);
        assert_eq!(loaded.last_artifact, state.last_artifact);
        match loaded.transform {
            Some(TransformSnapshot::Ingest { emit_sql, .. }) => assert!(emit_sql),
            _ => panic!("expected ingest snapshot"),
        }
        let _ = fs::remove_dir_all(&tmp);
    }

    #[test]
    fn load_missing_returns_none() {
        let tmp = std::env::temp_dir().join(format!("cohort-session-missing-{}", std::process::id()));
        let store = SessionStore::with_root(tmp);
        let id = SessionId::for_cwd(Path::new("/never/written"));
        assert!(store.load(&id).unwrap().is_none());
    }

    #[test]
    fn atomic_write_no_partial_on_repeat() {
        let tmp = std::env::temp_dir().join(format!("cohort-session-atomic-{}", std::process::id()));
        let store = SessionStore::with_root(tmp.clone());
        let id = SessionId::for_cwd(Path::new("/atomic"));
        let s1 = SessionState {
            cwd: PathBuf::from("/atomic"),
            last_artifact: None,
            transform: None,
        };
        store.save(&id, &s1).unwrap();
        let s2 = SessionState {
            cwd: PathBuf::from("/atomic"),
            last_artifact: Some(PathBuf::from("/atomic/y")),
            transform: None,
        };
        store.save(&id, &s2).unwrap();
        let loaded = store.load(&id).unwrap().unwrap();
        assert_eq!(loaded.last_artifact, Some(PathBuf::from("/atomic/y")));
        assert!(!store.dir_for(&id).join("state.json.tmp").exists());
        let _ = fs::remove_dir_all(&tmp);
    }
}
