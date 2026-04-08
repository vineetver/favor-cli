use std::fs::{self, File};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::CohortError;

#[allow(dead_code)]
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct StoreManifest {
    pub version: u32,
}

#[allow(dead_code)]
impl StoreManifest {
    pub const CURRENT_VERSION: u32 = 1;

    pub fn current() -> Self {
        Self { version: Self::CURRENT_VERSION }
    }
}

pub fn tmp_path(path: &Path) -> PathBuf {
    path.with_extension(format!(
        "{}.tmp",
        path.extension().and_then(|s| s.to_str()).unwrap_or("")
    ))
}

/// fsync the parent dir so the rename survives a crash. Best-effort: some
/// network mounts refuse to open a directory; skip rather than fail.
pub fn fsync_parent(path: &Path) {
    if let Some(parent) = path.parent() {
        if let Ok(dir) = File::open(parent) {
            let _ = dir.sync_all();
        }
    }
}

/// Write `bytes` to `path` atomically: write to `path.tmp`, fsync the
/// file, rename, then fsync the parent directory so the rename survives
/// a crash.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), CohortError> {
    use std::io::Write;
    let tmp = tmp_path(path);
    {
        let mut f = File::create(&tmp)
            .map_err(|e| CohortError::Resource(format!("create {}: {e}", tmp.display())))?;
        f.write_all(bytes)
            .map_err(|e| CohortError::Resource(format!("write {}: {e}", tmp.display())))?;
        f.sync_all()
            .map_err(|e| CohortError::Resource(format!("fsync {}: {e}", tmp.display())))?;
    }
    fs::rename(&tmp, path).map_err(|e| {
        CohortError::Resource(format!(
            "rename {} -> {}: {e}",
            tmp.display(),
            path.display()
        ))
    })?;
    fsync_parent(path);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_atomic_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("manifest.json");
        write_atomic(&path, b"hello world").unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"hello world");
        let entries: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(Result::ok)
            .collect();
        assert_eq!(entries.len(), 1);
    }

    #[test]
    fn tmp_path_appends_tmp_extension() {
        let p = tmp_path(Path::new("/foo/bar.json"));
        assert_eq!(p, PathBuf::from("/foo/bar.json.tmp"));
    }
}
