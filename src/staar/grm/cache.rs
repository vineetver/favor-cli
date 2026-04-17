//! GRM cache: fingerprint, probe, save, load.
//!
//! Layout under the cohort store:
//!   .cohort/cache/grm/<cohort_id>/<fingerprint>/
//!     grm.tsv           sample_i  sample_j  kinship
//!     pca_scores.tsv    sample_id  PC1  PC2  ...
//!     unrelated.txt     one sample_id per line
//!     manifest.json     { fingerprint, degree, n_pcs, created_at }

use std::io::Write;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::error::CohortError;

use super::types::{GrmArtifact, PcaScores, SparseGrm, UnrelatedSubset};

#[derive(Serialize, Deserialize)]
pub struct GrmManifest {
    pub fingerprint: String,
    pub degree: u8,
    pub n_pcs: usize,
    pub n_samples: usize,
    pub n_kinship_pairs: usize,
    pub n_unrelated: usize,
    pub created_at: String,
}

pub fn fingerprint(
    cohort_key: &str,
    king_seg_path: &Path,
    degree: u8,
    n_pcs: usize,
) -> Result<String, CohortError> {
    let seg_fp = crate::store::cohort::file_content_fingerprint(king_seg_path)?;
    let seg_hex = seg_fp.iter().map(|b| format!("{b:02x}")).collect::<String>();
    let input = format!("{cohort_key}|{seg_hex}|{degree}|{n_pcs}");
    Ok(crate::store::cohort::sha256_str(&input))
}

pub fn probe(dir: &Path) -> bool {
    dir.join("manifest.json").exists()
        && dir.join("grm.tsv").exists()
        && dir.join("pca_scores.tsv").exists()
}

pub fn save(
    dir: &Path,
    artifact: &GrmArtifact,
    fp: &str,
    degree: u8,
    n_pcs: usize,
) -> Result<(), CohortError> {
    std::fs::create_dir_all(dir)
        .map_err(|e| CohortError::Resource(format!("create {}: {e}", dir.display())))?;

    write_grm_tsv(&dir.join("grm.tsv"), &artifact.grm, &artifact.sample_ids)?;
    write_pca_tsv(
        &dir.join("pca_scores.tsv"),
        &artifact.pca,
        &artifact.sample_ids,
    )?;
    write_unrelated(&dir.join("unrelated.txt"), &artifact.unrelated, &artifact.sample_ids)?;

    let manifest = GrmManifest {
        fingerprint: fp.to_string(),
        degree,
        n_pcs,
        n_samples: artifact.sample_ids.len(),
        n_kinship_pairs: artifact.grm.triplets.len(),
        n_unrelated: artifact.unrelated.sample_indices.len(),
        created_at: chrono_now(),
    };
    let json = serde_json::to_string_pretty(&manifest)
        .map_err(|e| CohortError::Resource(format!("serialize manifest: {e}")))?;
    std::fs::write(dir.join("manifest.json"), json)
        .map_err(|e| CohortError::Resource(format!("write manifest: {e}")))?;

    Ok(())
}

fn write_grm_tsv(
    path: &Path,
    grm: &SparseGrm,
    sample_ids: &[String],
) -> Result<(), CohortError> {
    let mut f = std::fs::File::create(path)
        .map_err(|e| CohortError::Resource(format!("create {}: {e}", path.display())))?;
    writeln!(f, "ID1\tID2\tKinship")
        .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
    for &(i, j, k) in &grm.triplets {
        writeln!(f, "{}\t{}\t{k:.6}", sample_ids[i], sample_ids[j])
            .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
    }
    Ok(())
}

fn write_pca_tsv(
    path: &Path,
    pca: &PcaScores,
    sample_ids: &[String],
) -> Result<(), CohortError> {
    let n = pca.scores.nrows();
    let k = pca.scores.ncols();
    let mut f = std::fs::File::create(path)
        .map_err(|e| CohortError::Resource(format!("create {}: {e}", path.display())))?;
    let mut header = String::from("sample_id");
    for c in 0..k {
        header.push_str(&format!("\tPC{}", c + 1));
    }
    writeln!(f, "{header}")
        .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
    for (i, sid) in sample_ids.iter().enumerate().take(n) {
        write!(f, "{sid}")
            .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
        for c in 0..k {
            write!(f, "\t{:.6}", pca.scores[(i, c)])
                .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
        }
        writeln!(f)
            .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
    }
    Ok(())
}

fn write_unrelated(
    path: &Path,
    unrel: &UnrelatedSubset,
    sample_ids: &[String],
) -> Result<(), CohortError> {
    let mut f = std::fs::File::create(path)
        .map_err(|e| CohortError::Resource(format!("create {}: {e}", path.display())))?;
    for &idx in &unrel.sample_indices {
        writeln!(f, "{}", sample_ids[idx])
            .map_err(|e| CohortError::Resource(format!("write {}: {e}", path.display())))?;
    }
    Ok(())
}

pub fn grm_tsv_path(cache_dir: &Path) -> PathBuf {
    cache_dir.join("grm.tsv")
}

pub fn pca_tsv_path(cache_dir: &Path) -> PathBuf {
    cache_dir.join("pca_scores.tsv")
}

fn chrono_now() -> String {
    let dur = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default();
    format!("{}", dur.as_secs())
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::Mat;

    #[test]
    fn round_trip_grm_tsv() {
        let dir = tempfile::tempdir().unwrap();
        let grm = SparseGrm {
            triplets: vec![(0, 1, 0.25), (0, 0, 1.0), (1, 1, 1.0)],
            n_samples: 2,
        };
        let ids = vec!["s1".into(), "s2".into()];
        write_grm_tsv(&dir.path().join("grm.tsv"), &grm, &ids).unwrap();
        let content = std::fs::read_to_string(dir.path().join("grm.tsv")).unwrap();
        assert!(content.contains("s1\ts2\t0.250000"));
        assert!(content.starts_with("ID1\tID2\tKinship"));
    }

    #[test]
    fn round_trip_pca_tsv() {
        let dir = tempfile::tempdir().unwrap();
        let scores = Mat::from_fn(2, 3, |i, j| (i * 3 + j) as f64);
        let pca = PcaScores {
            scores,
            eigenvalues: vec![1.0, 0.5, 0.1],
        };
        let ids = vec!["s1".into(), "s2".into()];
        write_pca_tsv(&dir.path().join("pca.tsv"), &pca, &ids).unwrap();
        let content = std::fs::read_to_string(dir.path().join("pca.tsv")).unwrap();
        assert!(content.starts_with("sample_id\tPC1\tPC2\tPC3"));
        assert!(content.contains("s1\t0.000000\t1.000000\t2.000000"));
    }

    #[test]
    fn probe_returns_false_on_empty_dir() {
        let dir = tempfile::tempdir().unwrap();
        assert!(!probe(dir.path()));
    }
}
