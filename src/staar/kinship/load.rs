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

use faer::sparse::Triplet;
use faer::Mat;

use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::output::Output;
use crate::staar::genotype::GenotypeResult;
use crate::staar::kinship::types::{
    GroupPartition, KinshipMatrix, DENSE_DENSITY_THRESHOLD,
};

/// Parsed rows from one kinship TSV: symmetric triplet list, the count
/// of distinct off-diagonal pairs, and how many rows referenced samples
/// outside the genotype set. Shared by the sparse-direct and dense-
/// replay branches of `load_kinship`.
struct ParsedKinship {
    triplets: Vec<Triplet<u32, u32, f64>>,
    n_off_diag: usize,
    n_skipped: usize,
}

/// Parse a 3-column kinship TSV into a symmetric triplet list plus the
/// default-to-1 diagonal fill. Header detection, comment skip, and
/// sample-not-in-set handling match the original dense loader; only
/// the storage target (triplets vs dense `Mat`) is new.
fn parse_kinship_triplets(
    content: &str,
    path: &Path,
    id_to_idx: &HashMap<&str, usize>,
    n: usize,
) -> Result<ParsedKinship, CohortError> {
    let mut triplets: Vec<Triplet<u32, u32, f64>> = Vec::new();
    let mut diagonal_set = vec![false; n];
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
        let val: f64 = parts[2].parse().map_err(|_| {
            CohortError::Input(format!(
                "Kinship file '{}' line {}: cannot parse value '{}'",
                path.display(),
                line_no + 1,
                parts[2]
            ))
        })?;
        let (Some(&idx_i), Some(&idx_j)) = (id_to_idx.get(parts[0]), id_to_idx.get(parts[1]))
        else {
            n_skipped += 1;
            continue;
        };
        let ri = idx_i as u32;
        let rj = idx_j as u32;
        if idx_i == idx_j {
            diagonal_set[idx_i] = true;
            triplets.push(Triplet::new(ri, ri, val));
        } else {
            // Symmetric entries: push both halves so the triplet list is
            // a full-storage representation of the symmetric matrix.
            triplets.push(Triplet::new(ri, rj, val));
            triplets.push(Triplet::new(rj, ri, val));
            n_off_diag += 1;
        }
    }

    // Default-fill any diagonal the TSV didn't set. Matches the old
    // dense loader which pre-zeroed then wrote 1.0 on every diagonal
    // before parsing the file.
    for (i, &set) in diagonal_set.iter().enumerate() {
        if !set {
            triplets.push(Triplet::new(i as u32, i as u32, 1.0));
        }
    }

    Ok(ParsedKinship {
        triplets,
        n_off_diag,
        n_skipped,
    })
}

/// Read kinship matrices from one or more 3-column TSV files
/// (`sample_i  sample_j  kinship`). Reorders to `sample_order` and
/// validates every sample appears at least on the diagonal of each
/// loaded matrix.
///
/// Parses each file once into a symmetric triplet list, then picks the
/// storage variant from the measured density:
///
/// * density < `DENSE_DENSITY_THRESHOLD` → `KinshipMatrix::from_triplets`
///   builds a sparse CSC directly. Never touches a dense `Mat`, so a
///   pedigree with n = 50 000 samples and ~5 off-diagonal entries per
///   individual occupies ~3 MB instead of the ~20 GB staging buffer the
///   old dense-first path tried to allocate before the budget check.
/// * density ≥ threshold → allocate `Mat::<f64>::zeros(n, n)`, replay
///   the triplets, and call `KinshipMatrix::new`. The dense working set
///   is gated by the AI-REML memory budget in `kinship::budget`.
///
/// Both branches share the same parse pass so the header / comment /
/// default-diagonal / sample-not-in-set edge cases cannot drift.
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

        let parsed = parse_kinship_triplets(&content, path, &id_to_idx, n)?;

        out.status(&format!(
            "    {}: {} off-diagonal entries, {} entries skipped (sample not in genotype set)",
            label, parsed.n_off_diag, parsed.n_skipped,
        ));

        // Density accounts for both halves of every off-diagonal plus
        // the n diagonal entries we always push.
        let nnz = n + 2 * parsed.n_off_diag;
        let density = nnz as f64 / (n as f64 * n as f64);

        let kin = if density < DENSE_DENSITY_THRESHOLD {
            // Sparse-direct: never materialize n×n.
            KinshipMatrix::from_triplets(n, parsed.triplets, label.clone())?
        } else {
            // Dense replay: budget-bound, will be checked by the REML
            // dispatcher before the path is actually taken.
            let mut k = Mat::<f64>::zeros(n, n);
            for t in &parsed.triplets {
                k[(t.row as usize, t.col as usize)] = t.val;
            }
            KinshipMatrix::new(k, label.clone())?
        };
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

    let resolved = crate::staar::model::resolve_column(group_col, &cols, column_map, "Group")?;
    let id_col = crate::staar::model::resolve_id_column(&cols, column_map);
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
            let id = crate::staar::model::arrow_str(batch.column(0).as_ref(), row);
            let label = crate::staar::model::arrow_str(batch.column(1).as_ref(), row);
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

    let (order, assignments) =
        crate::staar::model::intern_labels(compact_labels.iter().map(|s| s.as_str()));
    let partition = GroupPartition::from_assignments(&assignments, &order)?;
    out.status(&format!(
        "  Heteroscedastic groups: {} levels ({})",
        partition.n_groups(),
        order.join(", "),
    ));
    Ok(partition)
}

/// Load the per-sample time vector for the random-slope LMM. Ports the
/// `time.var <- data[idx, random.slope]` extraction at upstream
/// `GMMAT/R/glmmkin.R:113`. Rows are returned in the same compacted
/// order `load_phenotype` produces (i.e., filtered by `pheno_mask` and
/// aligned with `geno.sample_names`), so `time_var[i]` aligns with the
/// i-th row of `y` / `x`.
pub fn load_random_slope(
    engine: &DfEngine,
    phenotype: &Path,
    time_col: &str,
    geno: &GenotypeResult,
    pheno_mask: &[bool],
    column_map: &HashMap<String, String>,
    out: &dyn Output,
) -> Result<Vec<f64>, CohortError> {
    engine.register_csv("_pheno_rs", phenotype, b'\t')?;
    let cols = engine.table_columns("_pheno_rs")?;

    let resolved =
        crate::staar::model::resolve_column(time_col, &cols, column_map, "RandomSlope")?;
    let id_col = crate::staar::model::resolve_id_column(&cols, column_map);
    let sample_list = geno
        .sample_names
        .iter()
        .map(|s| format!("'{s}'"))
        .collect::<Vec<_>>()
        .join(",");

    let sql = format!(
        "SELECT \"{id_col}\", CAST(\"{resolved}\" AS DOUBLE) \
         FROM _pheno_rs \
         WHERE \"{id_col}\" IN ({sample_list}) AND \"{resolved}\" IS NOT NULL"
    );
    let batches = engine.collect(&sql)?;

    let mut id_to_value: HashMap<String, f64> = HashMap::new();
    for batch in &batches {
        for row in 0..batch.num_rows() {
            let id = crate::staar::model::arrow_str(batch.column(0).as_ref(), row);
            let v = crate::staar::model::arrow_f64(batch.column(1).as_ref(), row);
            if !v.is_finite() {
                return Err(CohortError::Input(format!(
                    "random.slope column '{resolved}' has non-finite value for sample '{id}'"
                )));
            }
            id_to_value.insert(id, v);
        }
    }

    let mut out_vec: Vec<f64> = Vec::new();
    for (vcf_idx, name) in geno.sample_names.iter().enumerate() {
        if !pheno_mask[vcf_idx] {
            continue;
        }
        let v = *id_to_value.get(name).ok_or_else(|| {
            CohortError::Input(format!(
                "Sample '{name}' has no value in random.slope column '{resolved}'"
            ))
        })?;
        out_vec.push(v);
    }
    out.status(&format!(
        "  Random-slope time values: n = {}, range = [{}, {}]",
        out_vec.len(),
        out_vec.iter().cloned().fold(f64::INFINITY, f64::min),
        out_vec.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    ));
    Ok(out_vec)
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

    #[test]
    fn load_kinship_sparse_pedigree_takes_triplet_path() {
        // Pedigree fixture: 500 families × 2 sibs (n = 1000), each sib
        // pair carrying a single 0.5 off-diagonal entry. Overall density
        // is (1000 + 2·500) / 1_000_000 = 0.002, well under the 0.30
        // sparse threshold. The triplet-first loader must route this
        // straight into the sparse variant without ever allocating the
        // n×n staging matrix. This test would have driven the old loader
        // through an 8 MB scratch alloc; it now takes the sparse-direct
        // branch and never builds a Mat::zeros(1000, 1000).
        let dir = tempfile::tempdir().unwrap();
        let n_fam = 500;
        let n = n_fam * 2;
        let mut body = String::with_capacity(n_fam * 40);
        let mut order = Vec::with_capacity(n);
        for f in 0..n_fam {
            let a = format!("F{f}_A");
            let b = format!("F{f}_B");
            body.push_str(&format!("{a}\t{b}\t0.5\n"));
            order.push(a);
            order.push(b);
        }
        let path = write_tsv(&dir, "ped.tsv", &body);
        let out = null_output();
        let kins = load_kinship(&[path], &order, out.as_ref()).unwrap();
        assert_eq!(kins.len(), 1);
        let kin = &kins[0];
        assert!(
            kin.is_sparse(),
            "sparse pedigree must take the sparse-direct load path"
        );
        assert_eq!(kin.n(), n);
        let sparse = kin.as_sparse().unwrap();
        // Each family contributes 2 diagonal + 2 off-diagonal entries.
        let nnz = sparse.val().len();
        assert_eq!(nnz, 4 * n_fam, "nnz = {nnz}, expected {}", 4 * n_fam);
    }
}
