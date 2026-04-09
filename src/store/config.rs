//! Store configuration and path resolution.
//!
//! One resolver, one precedence order:
//!
//! 1. Explicit `--store-path <PATH>` CLI argument, if provided.
//! 2. `FAVOR_STORE` environment variable, if set and non-empty.
//! 3. Walk parent directories from cwd looking for an existing
//!    `.cohort/`. Stops at `$HOME` or filesystem root.
//! 4. Fall back to `<cwd>/.cohort/` (created lazily by `Store::open`).

use std::path::{Path, PathBuf};

use crate::error::CohortError;

const STORE_DIR_NAME: &str = ".cohort";
const ENV_VAR: &str = "FAVOR_STORE";

#[derive(Debug, Clone)]
pub struct StoreConfig {
    pub root: PathBuf,
}

impl StoreConfig {
    pub fn resolve(cli_arg: Option<PathBuf>) -> Result<Self, CohortError> {
        let cwd = std::env::current_dir().map_err(|e| {
            CohortError::Resource(format!("cannot determine current directory: {e}"))
        })?;
        let env = std::env::var_os(ENV_VAR);
        let root = resolve_store_path_with(cli_arg, env.as_deref().map(Path::new), &cwd)?;
        Ok(Self { root })
    }
}

/// Pure resolver — takes the cwd and env var explicitly so tests can
/// drive every arm without mutating process state.
pub(crate) fn resolve_store_path_with(
    cli_arg: Option<PathBuf>,
    env: Option<&Path>,
    cwd: &Path,
) -> Result<PathBuf, CohortError> {
    if let Some(p) = cli_arg {
        return Ok(absolutize(p, cwd));
    }
    if let Some(p) = env {
        if !p.as_os_str().is_empty() {
            return Ok(absolutize(p.to_path_buf(), cwd));
        }
    }
    if let Some(found) = walk_parents_for_store(cwd) {
        return Ok(found);
    }
    Ok(cwd.join(STORE_DIR_NAME))
}

fn absolutize(p: PathBuf, cwd: &Path) -> PathBuf {
    if p.is_absolute() {
        p
    } else {
        cwd.join(p)
    }
}

/// Walk from `start` upward looking for an existing `.cohort/` dir.
/// Stops at `$HOME` (inclusive) or filesystem root. Returns the
/// `.cohort/` path itself, not the parent.
fn walk_parents_for_store(start: &Path) -> Option<PathBuf> {
    let home = dirs::home_dir();
    let mut cur: Option<&Path> = Some(start);
    while let Some(dir) = cur {
        let candidate = dir.join(STORE_DIR_NAME);
        if candidate.is_dir() {
            return Some(candidate);
        }
        if let Some(h) = home.as_deref() {
            if dir == h {
                return None;
            }
        }
        cur = dir.parent();
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;

    #[test]
    fn cli_arg_wins_over_env_and_parent_walk() {
        let tmp = tempfile::tempdir().unwrap();
        let cwd = tmp.path().join("nested/deep");
        fs::create_dir_all(&cwd).unwrap();
        fs::create_dir_all(tmp.path().join(STORE_DIR_NAME)).unwrap();
        let env_path = tmp.path().join("env_store");
        let cli_path = tmp.path().join("cli_store");

        let resolved = resolve_store_path_with(
            Some(cli_path.clone()),
            Some(env_path.as_path()),
            &cwd,
        )
        .unwrap();
        assert_eq!(resolved, cli_path);
    }

    #[test]
    fn env_wins_over_parent_walk_and_fallback() {
        let tmp = tempfile::tempdir().unwrap();
        let cwd = tmp.path().join("nested");
        fs::create_dir_all(&cwd).unwrap();
        fs::create_dir_all(tmp.path().join(STORE_DIR_NAME)).unwrap();
        let env_path = tmp.path().join("env_store");

        let resolved =
            resolve_store_path_with(None, Some(env_path.as_path()), &cwd).unwrap();
        assert_eq!(resolved, env_path);
    }

    #[test]
    fn parent_walk_finds_existing_cohort() {
        let tmp = tempfile::tempdir().unwrap();
        let project = tmp.path().join("proj");
        let nested = project.join("a/b/c");
        fs::create_dir_all(&nested).unwrap();
        let store_dir = project.join(STORE_DIR_NAME);
        fs::create_dir_all(&store_dir).unwrap();

        let resolved = resolve_store_path_with(None, None, &nested).unwrap();
        assert_eq!(resolved, store_dir);
    }

    #[test]
    fn fallback_returns_cwd_cohort_when_nothing_found() {
        let tmp = tempfile::tempdir().unwrap();
        let cwd = tmp.path().join("isolated");
        fs::create_dir_all(&cwd).unwrap();

        let resolved = resolve_store_path_with(None, None, &cwd).unwrap();
        assert_eq!(resolved, cwd.join(STORE_DIR_NAME));
        assert!(!resolved.exists());
    }

    #[test]
    fn empty_env_falls_through_to_walk() {
        let tmp = tempfile::tempdir().unwrap();
        let cwd = tmp.path().join("c");
        fs::create_dir_all(&cwd).unwrap();
        fs::create_dir_all(tmp.path().join(STORE_DIR_NAME)).unwrap();

        let resolved =
            resolve_store_path_with(None, Some(Path::new("")), &cwd).unwrap();
        assert_eq!(resolved, tmp.path().join(STORE_DIR_NAME));
    }

    #[test]
    fn relative_cli_arg_resolved_against_cwd() {
        let tmp = tempfile::tempdir().unwrap();
        let cwd = tmp.path().join("c");
        fs::create_dir_all(&cwd).unwrap();

        let resolved =
            resolve_store_path_with(Some(PathBuf::from("rel/store")), None, &cwd).unwrap();
        assert_eq!(resolved, cwd.join("rel/store"));
    }
}
