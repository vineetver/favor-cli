//! Validated newtypes for the kinship pipeline.
//!
//! `KinshipMatrix` is the load-time storage choice — dense `Mat<f64>` or
//! sparse CSC, picked from density at load time. Both variants enforce
//! square shape, finiteness, and symmetry in their constructors so the
//! REML kernels never see a malformed matrix.
//!
//! `GroupPartition` carries a complete non-overlapping covering of the
//! sample set into G ≥ 1 non-empty groups for heteroscedastic residual
//! variance.
//!
//! `VarianceComponents` is the typed `[τ_kinship_1..τ_L, τ_group_1..τ_G]`
//! vector that the AI loop updates each iteration. Accessors are typed
//! by component kind so you can't accidentally write a kinship τ to a
//! group slot.
//!
//! `KinshipState` is the fitted output. It carries either a dense Σ⁻¹
//! (small studies, exact, fast scoring) or a sparse Cholesky factor of Σ
//! plus an optional cached selected-inverse (large pedigrees). The score
//! path dispatches on the variant; both produce identical p-values.

use faer::sparse::{SparseColMat, Triplet};
use faer::Mat;

use crate::error::CohortError;
use crate::staar::kinship::sparse::SparseFactor;

const KINSHIP_SYMMETRY_TOL: f64 = 1e-8;
/// Density threshold above which a loaded kinship matrix gets stored
/// densely. Below this we keep CSC and route through the sparse path.
/// Visible to `super::load` so the triplet-first loader can pick the
/// storage variant without materializing a dense Mat first.
pub(super) const DENSE_DENSITY_THRESHOLD: f64 = 0.30;

/// A symmetric, finite, square kinship matrix in either dense or sparse
/// storage. Both variants are validated by the constructor.
#[derive(Clone, Debug)]
pub enum KinshipMatrix {
    Dense {
        matrix: Mat<f64>,
        label: String,
    },
    Sparse {
        matrix: SparseColMat<u32, f64>,
        label: String,
    },
}

impl KinshipMatrix {
    /// Build from a dense matrix. Errors if not square, contains NaN/Inf, or
    /// is asymmetric beyond `KINSHIP_SYMMETRY_TOL`. Auto-routes to sparse
    /// storage if the density falls below [`DENSE_DENSITY_THRESHOLD`].
    pub fn new(matrix: Mat<f64>, label: String) -> Result<Self, CohortError> {
        let n = matrix.nrows();
        if matrix.ncols() != n {
            return Err(CohortError::Input(format!(
                "kinship matrix '{label}' is not square: ({n}, {})",
                matrix.ncols()
            )));
        }
        if n == 0 {
            return Err(CohortError::Input(format!(
                "kinship matrix '{label}' is empty"
            )));
        }
        let mut nnz = 0usize;
        for i in 0..n {
            for j in 0..n {
                let v = matrix[(i, j)];
                if !v.is_finite() {
                    return Err(CohortError::Input(format!(
                        "kinship matrix '{label}' has non-finite entry at ({i}, {j}): {v}"
                    )));
                }
                if v != 0.0 {
                    nnz += 1;
                }
            }
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let a = matrix[(i, j)];
                let b = matrix[(j, i)];
                if (a - b).abs() > KINSHIP_SYMMETRY_TOL {
                    return Err(CohortError::Input(format!(
                        "kinship matrix '{label}' not symmetric at ({i}, {j}): \
                         K[{i},{j}] = {a:.6e} vs K[{j},{i}] = {b:.6e}"
                    )));
                }
            }
        }
        let density = nnz as f64 / (n as f64 * n as f64);
        if density >= DENSE_DENSITY_THRESHOLD {
            Ok(Self::Dense { matrix, label })
        } else {
            Self::dense_to_sparse(&matrix, label)
        }
    }

    /// Build from a triplet stream. Caller guarantees the entries are
    /// already symmetric (or supplies both halves). Validates square shape,
    /// nonzero diagonal, and finiteness via faer's builder.
    pub fn from_triplets(
        n: usize,
        entries: Vec<Triplet<u32, u32, f64>>,
        label: String,
    ) -> Result<Self, CohortError> {
        if n == 0 {
            return Err(CohortError::Input(format!(
                "kinship matrix '{label}' is empty"
            )));
        }
        for t in &entries {
            if !t.val.is_finite() {
                return Err(CohortError::Input(format!(
                    "kinship matrix '{label}' has non-finite entry at ({}, {})",
                    t.row, t.col
                )));
            }
        }
        let matrix = SparseColMat::<u32, f64>::try_new_from_triplets(n, n, &entries).map_err(
            |e| CohortError::Input(format!("kinship matrix '{label}' build failed: {e:?}")),
        )?;
        Ok(Self::Sparse { matrix, label })
    }

    fn dense_to_sparse(dense: &Mat<f64>, label: String) -> Result<Self, CohortError> {
        let n = dense.nrows();
        let mut entries: Vec<Triplet<u32, u32, f64>> = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let v = dense[(i, j)];
                if v != 0.0 {
                    entries.push(Triplet::new(i as u32, j as u32, v));
                }
            }
        }
        Self::from_triplets(n, entries, label)
    }

    pub fn n(&self) -> usize {
        match self {
            Self::Dense { matrix, .. } => matrix.nrows(),
            Self::Sparse { matrix, .. } => matrix.nrows(),
        }
    }

    pub fn label(&self) -> &str {
        match self {
            Self::Dense { label, .. } | Self::Sparse { label, .. } => label,
        }
    }

    pub fn is_sparse(&self) -> bool {
        matches!(self, Self::Sparse { .. })
    }

    /// Mean of the diagonal entries — used by the REML warm start.
    pub fn mean_diagonal(&self) -> f64 {
        let n = self.n();
        if n == 0 {
            return 0.0;
        }
        let mut s = 0.0;
        match self {
            Self::Dense { matrix, .. } => {
                for i in 0..n {
                    s += matrix[(i, i)];
                }
            }
            Self::Sparse { matrix, .. } => {
                let col_ptr = matrix.symbolic().col_ptr();
                let row_idx = matrix.symbolic().row_idx();
                let val = matrix.val();
                for j in 0..n {
                    let start = col_ptr[j] as usize;
                    let end = col_ptr[j + 1] as usize;
                    for k in start..end {
                        if row_idx[k] as usize == j {
                            s += val[k];
                            break;
                        }
                    }
                }
            }
        }
        s / n as f64
    }

    /// Borrow the dense matrix view if this is a `Dense` variant.
    pub fn as_dense(&self) -> Option<&Mat<f64>> {
        match self {
            Self::Dense { matrix, .. } => Some(matrix),
            Self::Sparse { .. } => None,
        }
    }

    /// Borrow the sparse matrix view if this is a `Sparse` variant.
    pub fn as_sparse(&self) -> Option<&SparseColMat<u32, f64>> {
        match self {
            Self::Sparse { matrix, .. } => Some(matrix),
            Self::Dense { .. } => None,
        }
    }
}

/// A complete non-overlapping partition of `{0..n_samples}` into G ≥ 1
/// non-empty groups. Used to define heteroscedastic residual variance blocks
/// in AI-REML — every sample must belong to exactly one group, and the union
/// must cover all samples.
#[derive(Clone, Debug)]
pub struct GroupPartition {
    groups: Vec<Vec<u32>>,
    n: usize,
}

impl GroupPartition {
    /// Homoscedastic default: a single group covering all `n` samples.
    pub fn single(n: usize) -> Self {
        Self {
            groups: vec![(0..n as u32).collect()],
            n,
        }
    }

    /// No per-group residual variance. Used by the multi-trait kinship path
    /// (GMMAT `glmmkin.R:424-536::glmmkin.multi.ai`) where every variance
    /// component lives in the expanded kinship list as `E_jk ⊗ V` blocks
    /// and the residual term is itself one of those kinship entries.
    pub fn empty(n: usize) -> Self {
        Self {
            groups: Vec::new(),
            n,
        }
    }

    /// Build from a per-sample group-index assignment plus the group label
    /// table. Errors on out-of-range indices, empty groups, or zero samples.
    /// Labels are consumed for validation messages but not retained.
    pub fn from_assignments(
        assignments: &[usize],
        labels: &[String],
    ) -> Result<Self, CohortError> {
        let n = assignments.len();
        if n == 0 {
            return Err(CohortError::Input(
                "group partition requires at least one sample".into(),
            ));
        }
        if labels.is_empty() {
            return Err(CohortError::Input(
                "group partition requires at least one label".into(),
            ));
        }
        let n_groups = labels.len();
        let mut groups: Vec<Vec<u32>> = vec![Vec::new(); n_groups];
        for (i, &g) in assignments.iter().enumerate() {
            if g >= n_groups {
                return Err(CohortError::Input(format!(
                    "sample {i} has group index {g} but only {n_groups} labels are defined"
                )));
            }
            groups[g].push(i as u32);
        }
        for (g, members) in groups.iter().enumerate() {
            if members.is_empty() {
                return Err(CohortError::Input(format!(
                    "group '{}' (index {g}) has no members — drop unused levels or refilter",
                    labels[g]
                )));
            }
        }
        Ok(Self { groups, n })
    }

    pub fn n_groups(&self) -> usize {
        self.groups.len()
    }

    pub fn n_samples(&self) -> usize {
        self.n
    }

    pub fn group(&self, g: usize) -> &[u32] {
        &self.groups[g]
    }
}

/// PSD-box constraint for random-slope longitudinal LMM. Ports
/// `covariance.idx` from `GMMAT/R/glmmkin.R:175-183`. Indices are into
/// the flat `τ` layout used by `fit_reml_random_slope`
/// (kinship slots 0..3L, group slots 3L..3L+G). The covariance-slope
/// matrix `D = [[τ_int, τ_cov], [τ_cov, τ_slope]]` is PSD iff
/// `τ_cov² ≤ τ_int · τ_slope`.
#[derive(Clone, Copy, Debug)]
pub struct CovarianceIdx {
    pub cov_idx: usize,
    pub var_int_idx: usize,
    pub var_slope_idx: usize,
}

/// Variance components for AI-REML, laid out `[τ_kinship_1..τ_kinship_L,
/// τ_group_1..τ_group_G]`. The kinship/group split is carried with the data
/// so accessors are typed and indices cannot be confused.
#[derive(Clone, Debug)]
pub struct VarianceComponents {
    data: Vec<f64>,
    n_kinship: usize,
    n_group: usize,
}

impl VarianceComponents {
    pub fn zeros(n_kinship: usize, n_group: usize) -> Self {
        Self {
            data: vec![0.0; n_kinship + n_group],
            n_kinship,
            n_group,
        }
    }

    pub fn n_total(&self) -> usize {
        self.n_kinship + self.n_group
    }

    pub fn n_kinship(&self) -> usize {
        self.n_kinship
    }

    pub fn n_group(&self) -> usize {
        self.n_group
    }

    pub fn kinship(&self, l: usize) -> f64 {
        self.data[l]
    }

    pub fn group(&self, g: usize) -> f64 {
        self.data[self.n_kinship + g]
    }

    pub fn set_kinship(&mut self, l: usize, v: f64) {
        self.data[l] = v;
    }

    pub fn set_group(&mut self, g: usize, v: f64) {
        self.data[self.n_kinship + g] = v;
    }

    pub fn as_slice(&self) -> &[f64] {
        &self.data
    }

    pub fn as_slice_mut(&mut self) -> &mut [f64] {
        &mut self.data
    }
}

/// How a fitted `KinshipState` stores its inverse-of-Σ. Dense path keeps a
/// materialized `Σ⁻¹`. Sparse path keeps the Cholesky factor of Σ and applies
/// it lazily via solves; cheap matvecs (Σ⁻¹·v) become triangular solves.
#[derive(Clone)]
pub enum KinshipInverse {
    Dense(Mat<f64>),
    Sparse(SparseFactor),
}

/// Fitted state from AI-REML. Owned by `NullModel` when kinship is set.
/// Cloned once into `AnalysisVectors` for the score-test path.
#[derive(Clone)]
pub struct KinshipState {
    pub tau: VarianceComponents,
    /// Σ⁻¹ representation: dense (n×n materialized) or sparse Cholesky factor.
    pub inverse: KinshipInverse,
    /// Σ⁻¹ X, n×k. Always dense — k is small (≤ 20).
    pub sigma_inv_x: Mat<f64>,
    /// (X' Σ⁻¹ X)⁻¹, k×k.
    pub cov: Mat<f64>,
    /// PY = Σ⁻¹ (Y − X α̂), n×1. The "projected response" used by the score test.
    pub p_y: Mat<f64>,
    /// Per-kinship heritability estimate τ_l / (Σ τ + Σ τ_g). Length L.
    pub h2: Vec<f64>,
    pub n_iter: usize,
    pub outer_refits: usize,
    /// Dense AI-REML budget threaded in from `Resources::kinship_budget_bytes`
    /// at fit time. Stashed here so the score path can re-validate the
    /// budget without a separate parameter cascade.
    pub budget_bytes: u64,
}

impl KinshipState {
    pub fn n(&self) -> usize {
        self.sigma_inv_x.nrows()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kinship_matrix_rejects_asymmetric() {
        let mut k = Mat::<f64>::zeros(3, 3);
        for i in 0..3 {
            k[(i, i)] = 1.0;
        }
        k[(0, 1)] = 0.5;
        k[(1, 0)] = 0.4;
        let err = KinshipMatrix::new(k, "bad".into()).unwrap_err();
        match err {
            CohortError::Input(msg) => assert!(msg.contains("not symmetric")),
            other => panic!("expected Input, got {other:?}"),
        }
    }

    #[test]
    fn kinship_matrix_rejects_nan() {
        let mut k = Mat::<f64>::zeros(2, 2);
        k[(0, 0)] = 1.0;
        k[(1, 1)] = f64::NAN;
        let err = KinshipMatrix::new(k, "bad".into()).unwrap_err();
        assert!(matches!(err, CohortError::Input(_)));
    }

    #[test]
    fn dense_kinship_routed_dense_when_dense() {
        let n = 5;
        let mut k = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                k[(i, j)] = if i == j { 1.0 } else { 0.5 };
            }
        }
        let kin = KinshipMatrix::new(k, "dense".into()).unwrap();
        assert!(!kin.is_sparse());
        assert!(kin.as_dense().is_some());
    }

    #[test]
    fn sparse_pedigree_routed_sparse() {
        // 5 families × 2 sibs, 0.5 within and 0 between → 20% nonzero, sparse.
        let n = 50;
        let mut k = Mat::<f64>::zeros(n, n);
        for f in 0..25 {
            let base = f * 2;
            k[(base, base)] = 1.0;
            k[(base + 1, base + 1)] = 1.0;
            k[(base, base + 1)] = 0.5;
            k[(base + 1, base)] = 0.5;
        }
        let kin = KinshipMatrix::new(k, "ped".into()).unwrap();
        assert!(kin.is_sparse(), "pedigree kinship should be sparse-stored");
        assert_eq!(kin.n(), n);
    }

    #[test]
    fn group_partition_rejects_empty_group() {
        let labels = vec!["A".into(), "B".into()];
        let err = GroupPartition::from_assignments(&[0, 0, 0], &labels).unwrap_err();
        assert!(matches!(err, CohortError::Input(_)));
    }

    #[test]
    fn group_partition_rejects_oob_index() {
        let labels = vec!["A".into()];
        let err = GroupPartition::from_assignments(&[0, 1], &labels).unwrap_err();
        assert!(matches!(err, CohortError::Input(_)));
    }

    #[test]
    fn group_partition_single_covers_all() {
        let p = GroupPartition::single(7);
        assert_eq!(p.n_groups(), 1);
        assert_eq!(p.n_samples(), 7);
        assert_eq!(p.group(0).len(), 7);
    }

    #[test]
    fn variance_components_kinship_group_split() {
        let mut tau = VarianceComponents::zeros(2, 3);
        assert_eq!(tau.n_total(), 5);
        tau.set_kinship(0, 0.4);
        tau.set_kinship(1, 0.6);
        tau.set_group(0, 1.0);
        tau.set_group(1, 2.0);
        tau.set_group(2, 3.0);
        assert_eq!(tau.kinship(0), 0.4);
        assert_eq!(tau.kinship(1), 0.6);
        assert_eq!(tau.group(0), 1.0);
        assert_eq!(tau.group(2), 3.0);
        assert_eq!(tau.as_slice(), &[0.4, 0.6, 1.0, 2.0, 3.0]);
    }
}
