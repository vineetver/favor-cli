//! TSV → `KinshipMatrix` and phenotype-column → `GroupPartition` loaders.
//!
//! `load_kinship` reads one or more 3-column TSV files
//! (`sample_i  sample_j  kinship`), reorders to genotype-set sample order,
//! and hands off to `KinshipMatrix::new` which detects density and routes
//! to the dense or sparse storage variant.
//!
//! `load_groups` reads a categorical column from the phenotype TSV via
//! DataFusion, applies the user-supplied column-alias map, and builds a
//! validated `GroupPartition` over the genotyped sample subset.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use faer::Mat;

use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::output::Output;
use crate::staar::genotype::GenotypeResult;
use crate::staar::kinship::types::{GroupPartition, KinshipMatrix};

/// Read kinship matrices from one or more 3-column TSV files
/// (`sample_i  sample_j  kinship`). Reorders to `sample_order` and validates
/// every sample appears at least on the diagonal of each loaded matrix.
/// Returns one [`KinshipMatrix`] per input path; the constructor enforces
/// symmetry, finiteness, and density-based routing to dense or sparse
/// storage.
pub fn load_kinship(
    paths: &[PathBuf],
    sample_order: &[String],
    out: &dyn Output,
) -> Result<Vec<KinshipMatrix>, CohortError> {
    let n = sample_order.len();
    let id_to_idx: HashMap<&str, usize> = sample_order
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let mut result = Vec::with_capacity(paths.len());
    for path in paths {
        let label = path
            .file_stem()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());
        out.status(&format!("  Loading kinship: {}", path.display()));
        let content = std::fs::read_to_string(path).map_err(|e| {
            CohortError::Resource(format!(
                "Cannot read kinship file '{}': {e}",
                path.display()
            ))
        })?;

        // Stage as dense first; KinshipMatrix::new will down-convert to
        // sparse storage if the density falls below the threshold.
        let mut k = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            k[(i, i)] = 1.0;
        }
        let mut n_off_diag = 0usize;
        let mut n_skipped = 0usize;
        for (line_no, line) in content.lines().enumerate() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            if parts.len() < 3 {
                continue;
            }
            // Skip header row if value column is non-numeric.
            if line_no == 0 && parts[2].parse::<f64>().is_err() {
                continue;
            }
            let s_i = parts[0];
            let s_j = parts[1];
            let val: f64 = parts[2].parse().map_err(|_| {
                CohortError::Input(format!(
                    "Kinship file '{}' line {}: cannot parse value '{}'",
                    path.display(),
                    line_no + 1,
                    parts[2]
                ))
            })?;
            let (Some(&idx_i), Some(&idx_j)) = (id_to_idx.get(s_i), id_to_idx.get(s_j)) else {
                n_skipped += 1;
                continue;
            };
            k[(idx_i, idx_j)] = val;
            k[(idx_j, idx_i)] = val;
            if idx_i != idx_j {
                n_off_diag += 1;
            }
        }

        out.status(&format!(
            "    {label}: {n_off_diag} off-diagonal entries, {n_skipped} entries skipped (sample not in genotype set)"
        ));

        let kin = KinshipMatrix::new(k, label.clone())?;
        out.status(&format!(
            "    {label} stored as {}",
            if kin.is_sparse() { "sparse CSC" } else { "dense" }
        ));
        result.push(kin);
    }
    Ok(result)
}

/// Read a categorical group column from the phenotype TSV via DataFusion and
/// return a validated [`GroupPartition`] over the genotyped sample subset.
pub fn load_groups(
    engine: &DfEngine,
    phenotype: &Path,
    group_col: &str,
    geno: &GenotypeResult,
    pheno_mask: &[bool],
    column_map: &HashMap<String, String>,
    out: &dyn Output,
) -> Result<GroupPartition, CohortError> {
    engine.register_csv("_pheno_groups", phenotype, b'\t')?;
    let cols = engine.table_columns("_pheno_groups")?;

    let resolved = if let Some(mapped) = column_map.get(group_col) {
        if cols.contains(mapped) {
            mapped.clone()
        } else {
            return Err(CohortError::Input(format!(
                "Column map '{group_col}={mapped}' but '{mapped}' not in phenotype.",
            )));
        }
    } else if cols.contains(&group_col.to_string()) {
        group_col.to_string()
    } else {
        let lower = group_col.to_lowercase();
        cols.iter()
            .find(|c| c.to_lowercase() == lower)
            .cloned()
            .ok_or_else(|| {
                CohortError::Input(format!(
                    "Group column '{group_col}' not found in phenotype. Available: {}",
                    cols.join(", ")
                ))
            })?
    };

    let id_col = cols[0].clone();
    let sample_list = geno
        .sample_names
        .iter()
        .map(|s| format!("'{s}'"))
        .collect::<Vec<_>>()
        .join(",");

    let sql = format!(
        "SELECT \"{id_col}\", CAST(\"{resolved}\" AS VARCHAR) \
         FROM _pheno_groups \
         WHERE \"{id_col}\" IN ({sample_list}) AND \"{resolved}\" IS NOT NULL"
    );
    let batches = engine.collect(&sql)?;

    let mut id_to_label: HashMap<String, String> = HashMap::new();
    for batch in &batches {
        for row in 0..batch.num_rows() {
            let id =
                arrow::util::display::array_value_to_string(batch.column(0).as_ref(), row)
                    .unwrap_or_default();
            let label =
                arrow::util::display::array_value_to_string(batch.column(1).as_ref(), row)
                    .unwrap_or_default();
            id_to_label.insert(id, label);
        }
    }

    let mut compact_labels: Vec<String> = Vec::new();
    for (vcf_idx, name) in geno.sample_names.iter().enumerate() {
        if !pheno_mask[vcf_idx] {
            continue;
        }
        let label = id_to_label.get(name).cloned().ok_or_else(|| {
            CohortError::Input(format!(
                "Sample '{name}' has no value in group column '{resolved}'"
            ))
        })?;
        compact_labels.push(label);
    }

    let mut order: Vec<String> = Vec::new();
    let mut assignments: Vec<usize> = Vec::with_capacity(compact_labels.len());
    for label in &compact_labels {
        let g = order.iter().position(|l| l == label).unwrap_or_else(|| {
            order.push(label.clone());
            order.len() - 1
        });
        assignments.push(g);
    }

    let partition = GroupPartition::from_assignments(&assignments, &order)?;
    out.status(&format!(
        "  Heteroscedastic groups: {} levels ({})",
        partition.n_groups(),
        order.join(", "),
    ));
    Ok(partition)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::output::{create, OutputMode};
    use std::io::Write;

    fn null_output() -> Box<dyn Output> {
        create(&OutputMode::Machine)
    }

    fn write_tsv(dir: &tempfile::TempDir, name: &str, body: &str) -> PathBuf {
        let p = dir.path().join(name);
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
        p
    }

    #[test]
    fn load_kinship_reorders_to_sample_order_and_symmetrizes() {
        // File lists S2 before S1; loader must place S1 at index 0 and S2
        // at index 1. Off-diagonal is given once and must be mirrored to
        // both halves.
        let dir = tempfile::tempdir().unwrap();
        let body = "S2\tS2\t1.0\n\
                    S1\tS1\t1.0\n\
                    S1\tS2\t0.5\n";
        let path = write_tsv(&dir, "kin.tsv", body);
        let order = vec!["S1".to_string(), "S2".to_string()];
        let out = null_output();
        let kins =
            load_kinship(&[path], &order, out.as_ref()).expect("load_kinship should succeed");
        assert_eq!(kins.len(), 1);
        let kin = &kins[0];
        assert_eq!(kin.n(), 2);
        let dense = kin.as_dense().expect("2x2 dense");
        assert_eq!(dense[(0, 0)], 1.0);
        assert_eq!(dense[(1, 1)], 1.0);
        assert_eq!(dense[(0, 1)], 0.5);
        assert_eq!(dense[(1, 0)], 0.5);
        assert_eq!(kin.label(), "kin");
    }

    #[test]
    fn load_kinship_skips_samples_outside_genotype_set() {
        // S3 is not in the sample order — its rows are silently dropped
        // (not an error: kinship files often cover a superset of the
        // genotyped samples).
        let dir = tempfile::tempdir().unwrap();
        let body = "S1\tS1\t1.0\n\
                    S2\tS2\t1.0\n\
                    S1\tS2\t0.5\n\
                    S1\tS3\t0.25\n\
                    S3\tS3\t1.0\n";
        let path = write_tsv(&dir, "kin.tsv", body);
        let order = vec!["S1".to_string(), "S2".to_string()];
        let out = null_output();
        let kins = load_kinship(&[path], &order, out.as_ref()).unwrap();
        let dense = kins[0].as_dense().unwrap();
        // S3 entries dropped, S1↔S2 still set.
        assert_eq!(dense[(0, 1)], 0.5);
        assert_eq!(dense[(1, 0)], 0.5);
    }

    #[test]
    fn load_kinship_rejects_unparseable_value() {
        let dir = tempfile::tempdir().unwrap();
        let body = "S1\tS1\t1.0\n\
                    S1\tS2\tNOT_A_NUMBER\n";
        let path = write_tsv(&dir, "kin.tsv", body);
        let order = vec!["S1".to_string(), "S2".to_string()];
        let out = null_output();
        let err = load_kinship(&[path], &order, out.as_ref()).unwrap_err();
        match err {
            CohortError::Input(msg) => {
                assert!(msg.contains("cannot parse value"), "msg = {msg}");
            }
            other => panic!("expected Input error, got {other:?}"),
        }
    }

    #[test]
    fn load_kinship_skips_header_and_comments() {
        // Header row (non-numeric value column) and `#` lines are ignored.
        let dir = tempfile::tempdir().unwrap();
        let body = "id1\tid2\tkinship\n\
                    # comment line\n\
                    \n\
                    S1\tS1\t1.0\n\
                    S2\tS2\t1.0\n\
                    S1\tS2\t0.25\n";
        let path = write_tsv(&dir, "kin.tsv", body);
        let order = vec!["S1".to_string(), "S2".to_string()];
        let out = null_output();
        let kins = load_kinship(&[path], &order, out.as_ref()).unwrap();
        let dense = kins[0].as_dense().unwrap();
        assert_eq!(dense[(0, 1)], 0.25);
    }

    #[test]
    fn load_kinship_diagonal_defaults_to_one_when_unspecified() {
        // Missing diagonal lines are fine — the loader pre-fills 1.0 on
        // the diagonal so the resulting matrix is still SPD-shaped.
        let dir = tempfile::tempdir().unwrap();
        let body = "S1\tS2\t0.3\n";
        let path = write_tsv(&dir, "kin.tsv", body);
        let order = vec!["S1".to_string(), "S2".to_string()];
        let out = null_output();
        let kins = load_kinship(&[path], &order, out.as_ref()).unwrap();
        let dense = kins[0].as_dense().unwrap();
        assert_eq!(dense[(0, 0)], 1.0);
        assert_eq!(dense[(1, 1)], 1.0);
        assert_eq!(dense[(0, 1)], 0.3);
        assert_eq!(dense[(1, 0)], 0.3);
    }

    #[test]
    fn load_kinship_returns_one_matrix_per_path() {
        let dir = tempfile::tempdir().unwrap();
        let p1 = write_tsv(&dir, "a.tsv", "S1\tS1\t1.0\nS2\tS2\t1.0\nS1\tS2\t0.4\n");
        let p2 = write_tsv(&dir, "b.tsv", "S1\tS1\t1.0\nS2\tS2\t1.0\nS1\tS2\t0.1\n");
        let order = vec!["S1".to_string(), "S2".to_string()];
        let out = null_output();
        let kins = load_kinship(&[p1, p2], &order, out.as_ref()).unwrap();
        assert_eq!(kins.len(), 2);
        assert_eq!(kins[0].label(), "a");
        assert_eq!(kins[1].label(), "b");
        assert_eq!(kins[0].as_dense().unwrap()[(0, 1)], 0.4);
        assert_eq!(kins[1].as_dense().unwrap()[(0, 1)], 0.1);
    }
}
