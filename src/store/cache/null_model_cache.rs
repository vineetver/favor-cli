//! Null model disk cache: persist fitted NullModel so re-runs with the
//! same cohort/trait/covariates skip the fitting step.
//!
//! Disk format (`null_model.bin`):
//!   Header (64 bytes): magic, version, n_samples, n_covariates, sigma2, flags
//!   Body: residuals[n], x_matrix[n*k], xtx_inv[k*k] as raw f64 LE
//!   If has_fitted: fitted_values[n]
//!   If has_working_weights: working_weights[n]
//!
//! Kinship-aware models are NOT cached in v1 — the KinshipState contains a
//! sparse Cholesky factor that is not trivially serializable. A warning is
//! logged and the null model is refit.

use std::io::{Read as IoRead, Write as IoWrite};
use std::path::{Path, PathBuf};

use faer::Mat;
use sha2::{Digest, Sha256};

use crate::error::CohortError;
use crate::staar::model::NullModel;

const MAGIC: &[u8; 8] = b"FVNULLM1";
const VERSION: u16 = 1;
const HEADER_SIZE: usize = 64;

const FLAG_HAS_FITTED: u8 = 1;
const FLAG_HAS_WORKING_WEIGHTS: u8 = 2;

/// Cache key for the null model, following the same hashing pattern as
/// `score_cache::cache_key`. Keyed by (store manifest key, trait, sorted
/// covariates, known-loci content, kinship paths, kinship groups).
pub fn cache_key(
    store_key: &str,
    trait_name: &str,
    covariates: &[String],
    known_loci: Option<&Path>,
    kinship: &[PathBuf],
    kinship_groups: Option<&str>,
) -> String {
    let mut hasher = Sha256::new();
    hasher.update(b"null_model|");
    hasher.update(store_key.as_bytes());
    hasher.update(b"|");
    hasher.update(trait_name.as_bytes());
    hasher.update(b"|");
    let mut sorted = covariates.to_vec();
    sorted.sort();
    for cov in &sorted {
        hasher.update(cov.as_bytes());
        hasher.update(b",");
    }
    hasher.update(b"|");
    if let Some(p) = known_loci {
        hasher.update(b"loci=");
        match std::fs::read(p) {
            Ok(bytes) => hasher.update(&bytes),
            Err(_) => hasher.update(p.to_string_lossy().as_bytes()),
        }
    }
    hasher.update(b"|");
    let mut canonical: Vec<String> = kinship
        .iter()
        .map(|p| {
            std::fs::canonicalize(p)
                .map(|c| c.to_string_lossy().into_owned())
                .unwrap_or_else(|_| p.to_string_lossy().into_owned())
        })
        .collect();
    canonical.sort();
    for c in &canonical {
        hasher.update(c.as_bytes());
        hasher.update(b",");
    }
    hasher.update(b"|");
    if let Some(g) = kinship_groups {
        hasher.update(g.as_bytes());
    }
    format!("{:x}", hasher.finalize())
}

/// Check whether a cached null model exists and has a valid header.
pub fn probe(dir: &Path) -> bool {
    let path = dir.join("null_model.bin");
    let Ok(meta) = std::fs::metadata(&path) else {
        return false;
    };
    if meta.len() < HEADER_SIZE as u64 {
        return false;
    }
    let Ok(mut f) = std::fs::File::open(&path) else {
        return false;
    };
    let mut header = [0u8; HEADER_SIZE];
    if f.read_exact(&mut header).is_err() {
        return false;
    }
    &header[0..8] == MAGIC && u16::from_le_bytes([header[8], header[9]]) == VERSION
}

/// Serialize a NullModel to disk. Atomic: writes .tmp then renames.
pub fn save(dir: &Path, model: &NullModel) -> Result<(), CohortError> {
    std::fs::create_dir_all(dir)
        .map_err(|e| CohortError::Resource(format!("create {}: {e}", dir.display())))?;

    let path = dir.join("null_model.bin");
    let tmp = dir.join("null_model.bin.tmp");

    let n = model.n_samples;
    let k = model.x_matrix.ncols();
    let mut flags: u8 = 0;
    if model.fitted_values.is_some() {
        flags |= FLAG_HAS_FITTED;
    }
    if model.working_weights.is_some() {
        flags |= FLAG_HAS_WORKING_WEIGHTS;
    }

    let mut w = std::io::BufWriter::new(
        std::fs::File::create(&tmp)
            .map_err(|e| CohortError::Resource(format!("create {}: {e}", tmp.display())))?,
    );

    // Header
    let mut header = [0u8; HEADER_SIZE];
    header[0..8].copy_from_slice(MAGIC);
    header[8..10].copy_from_slice(&VERSION.to_le_bytes());
    header[10..14].copy_from_slice(&(n as u32).to_le_bytes());
    header[14..18].copy_from_slice(&(k as u32).to_le_bytes());
    header[18..26].copy_from_slice(&model.sigma2.to_le_bytes());
    header[26] = flags;
    w.write_all(&header)
        .map_err(|e| CohortError::Resource(format!("write header: {e}")))?;

    // residuals: n × 1 column-major (faer stores column-major)
    write_mat(&mut w, &model.residuals)?;
    // x_matrix: n × k
    write_mat(&mut w, &model.x_matrix)?;
    // xtx_inv: k × k
    write_mat(&mut w, &model.xtx_inv)?;

    if let Some(ref fv) = model.fitted_values {
        write_vec(&mut w, fv)?;
    }
    if let Some(ref ww) = model.working_weights {
        write_vec(&mut w, ww)?;
    }

    w.flush()
        .map_err(|e| CohortError::Resource(format!("flush: {e}")))?;
    drop(w);

    // fsync + rename for atomicity
    #[cfg(unix)]
    {
        use std::os::unix::io::AsRawFd;
        let f = std::fs::File::open(&tmp)
            .map_err(|e| CohortError::Resource(format!("reopen for fsync: {e}")))?;
        unsafe { libc::fsync(f.as_raw_fd()) };
    }

    std::fs::rename(&tmp, &path)
        .map_err(|e| CohortError::Resource(format!("rename {}: {e}", tmp.display())))?;

    Ok(())
}

/// Deserialize a NullModel from disk.
pub fn load(dir: &Path) -> Result<NullModel, CohortError> {
    let path = dir.join("null_model.bin");
    let mut f = std::io::BufReader::new(
        std::fs::File::open(&path)
            .map_err(|e| CohortError::DataMissing(format!("open {}: {e}", path.display())))?,
    );

    let mut header = [0u8; HEADER_SIZE];
    f.read_exact(&mut header)
        .map_err(|e| CohortError::DataMissing(format!("read header: {e}")))?;

    if &header[0..8] != MAGIC {
        return Err(CohortError::DataMissing(format!(
            "{}: bad magic",
            path.display()
        )));
    }
    let version = u16::from_le_bytes([header[8], header[9]]);
    if version != VERSION {
        return Err(CohortError::DataMissing(format!(
            "{}: version {version}, expected {VERSION}",
            path.display()
        )));
    }

    let n = u32::from_le_bytes(header[10..14].try_into().unwrap()) as usize;
    let k = u32::from_le_bytes(header[14..18].try_into().unwrap()) as usize;
    let sigma2 = f64::from_le_bytes(header[18..26].try_into().unwrap());
    let flags = header[26];

    let residuals = read_mat(&mut f, n, 1)?;
    let x_matrix = read_mat(&mut f, n, k)?;
    let xtx_inv = read_mat(&mut f, k, k)?;

    let fitted_values = if flags & FLAG_HAS_FITTED != 0 {
        Some(read_vec(&mut f, n)?)
    } else {
        None
    };

    let working_weights = if flags & FLAG_HAS_WORKING_WEIGHTS != 0 {
        Some(read_vec(&mut f, n)?)
    } else {
        None
    };

    Ok(NullModel {
        residuals,
        x_matrix,
        xtx_inv,
        sigma2,
        n_samples: n,
        fitted_values,
        working_weights,
        kinship: None,
    })
}

fn write_mat(w: &mut impl IoWrite, m: &Mat<f64>) -> Result<(), CohortError> {
    for j in 0..m.ncols() {
        for i in 0..m.nrows() {
            w.write_all(&m[(i, j)].to_le_bytes())
                .map_err(|e| CohortError::Resource(format!("write mat: {e}")))?;
        }
    }
    Ok(())
}

fn write_vec(w: &mut impl IoWrite, v: &[f64]) -> Result<(), CohortError> {
    for val in v {
        w.write_all(&val.to_le_bytes())
            .map_err(|e| CohortError::Resource(format!("write vec: {e}")))?;
    }
    Ok(())
}

fn read_mat(r: &mut impl IoRead, nrows: usize, ncols: usize) -> Result<Mat<f64>, CohortError> {
    let mut m = Mat::zeros(nrows, ncols);
    let mut buf = [0u8; 8];
    for j in 0..ncols {
        for i in 0..nrows {
            r.read_exact(&mut buf)
                .map_err(|e| CohortError::DataMissing(format!("read mat: {e}")))?;
            m[(i, j)] = f64::from_le_bytes(buf);
        }
    }
    Ok(m)
}

fn read_vec(r: &mut impl IoRead, n: usize) -> Result<Vec<f64>, CohortError> {
    let mut v = Vec::with_capacity(n);
    let mut buf = [0u8; 8];
    for _ in 0..n {
        r.read_exact(&mut buf)
            .map_err(|e| CohortError::DataMissing(format!("read vec: {e}")))?;
        v.push(f64::from_le_bytes(buf));
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_model(n: usize, k: usize) -> NullModel {
        let mut residuals = Mat::zeros(n, 1);
        for i in 0..n {
            residuals[(i, 0)] = (i as f64) * 0.1;
        }
        let mut x = Mat::zeros(n, k);
        for i in 0..n {
            for j in 0..k {
                x[(i, j)] = (i * k + j) as f64;
            }
        }
        let xtx_inv = Mat::identity(k, k);
        NullModel {
            residuals,
            x_matrix: x,
            xtx_inv,
            sigma2: 1.234,
            n_samples: n,
            fitted_values: None,
            working_weights: None,
            kinship: None,
        }
    }

    #[test]
    fn round_trip_continuous() {
        let dir = tempfile::tempdir().unwrap();
        let model = dummy_model(100, 3);
        save(dir.path(), &model).unwrap();
        assert!(probe(dir.path()));

        let loaded = load(dir.path()).unwrap();
        assert_eq!(loaded.n_samples, 100);
        assert_eq!(loaded.x_matrix.ncols(), 3);
        assert!((loaded.sigma2 - 1.234).abs() < 1e-12);
        assert!(loaded.fitted_values.is_none());
        assert!(loaded.working_weights.is_none());
        assert!(loaded.kinship.is_none());

        for i in 0..100 {
            assert!((loaded.residuals[(i, 0)] - (i as f64) * 0.1).abs() < 1e-12);
        }
    }

    #[test]
    fn round_trip_binary() {
        let dir = tempfile::tempdir().unwrap();
        let mut model = dummy_model(50, 2);
        model.fitted_values = Some(vec![0.5; 50]);
        model.working_weights = Some(vec![0.25; 50]);
        save(dir.path(), &model).unwrap();

        let loaded = load(dir.path()).unwrap();
        assert_eq!(loaded.fitted_values.as_ref().unwrap().len(), 50);
        assert_eq!(loaded.working_weights.as_ref().unwrap().len(), 50);
        assert!((loaded.fitted_values.unwrap()[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn probe_missing_returns_false() {
        let dir = tempfile::tempdir().unwrap();
        assert!(!probe(dir.path()));
    }

    #[test]
    fn cache_key_differs_by_trait() {
        let k1 = cache_key("store", "BMI", &[], None, &[], None);
        let k2 = cache_key("store", "HEIGHT", &[], None, &[], None);
        assert_ne!(k1, k2);
    }

    #[test]
    fn cache_key_covariate_order_irrelevant() {
        let k1 = cache_key("s", "t", &["age".into(), "sex".into()], None, &[], None);
        let k2 = cache_key("s", "t", &["sex".into(), "age".into()], None, &[], None);
        assert_eq!(k1, k2);
    }
}
