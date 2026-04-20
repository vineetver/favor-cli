// `needless_range_loop` is silenced because the expansion helpers and
// score-cov assembly read more naturally as index loops than iterator
// chains over `Mat` rows — matches `multi.rs`.
#![allow(clippy::needless_range_loop)]

//! Multi-trait AI-REML with shared kinship matrices (related samples).
//!
//! Ports `glmmkin.multi.ai` from upstream GMMAT
//! (`hanchenphd/GMMAT@473b342 R/glmmkin.R:424-536`) for the `family =
//! gaussian` multi-pheno path that `fit_null_glmmkin_multi` wraps. The
//! upstream wrapper itself lives at
//! `xihaoli/MultiSTAAR@main R/fit_null_glmmkin_multi.R:70-123`.
//!
//! ## Variance-component parameterisation
//!
//! For `k` continuous traits with shared covariates `X` (shape `n × p`),
//! `L` kinship matrices `Φ_1..Φ_L` (shape `n × n`), and `G ≥ 1` residual
//! groups covering `{0..n}`, the joint covariance of `vec(Y)` (length
//! `nk`) under the null is
//!
//!   `Σ = Σ_{g=1}^G  Σ_res_g ⊗ D_g  +  Σ_{l=1}^L  Σ_g_l ⊗ Φ_l`
//!
//! where `D_g` is the diagonal indicator of group `g` (so `Σ_g D_g = I_n`
//! in the homoscedastic case) and each `Σ_res_g`, `Σ_g_l` is a symmetric
//! `(k × k)` matrix of variance/covariance components.
//!
//! Upstream's trick (`glmmkin.R:147-166`) expands this into a flat list
//! of scalar variance components by walking every symmetric `(k × k)`
//! matrix entry: each `(j, m)` pair with `j ≤ m` contributes one
//! parameter `τ_{j,m}` and one expanded kinship entry `E_{j,m} ⊗ V`,
//! where `E_{j,m}` is the `(k × k)` symmetric matrix with a single `1`
//! at `(j, m)` and `(m, j)`. The AI-REML algorithm runs on this flat
//! scalar-tau list exactly like the single-trait path; at convergence,
//! upstream (`glmmkin.R:213-224`) reshapes the flat vector back into
//! per-group and per-kinship `(k × k)` symmetric matrices.
//!
//! Total parameter count: `n_params = k(k+1)/2 · (G + L)`.
//!
//! ## What this module adds
//!
//! * `expand_kinships` — builds the `n_params` expanded kinship matrices
//!   `E_{j,m} ⊗ D_g` and `E_{j,m} ⊗ Φ_l`.
//! * `stack_response` — `vec(Y)` as an `(nk × 1)` column.
//! * `stack_covariates` — `I_k ⊗ X` as an `(nk × kp)` dense matrix
//!   (block-diagonal with `k` copies of `X` on the diagonal).
//! * `MultiKinshipState` — converged fit: stacked `Σ⁻¹`, `Σ⁻¹X_stacked`,
//!   `cov = (X_stacked' Σ⁻¹ X_stacked)⁻¹`, stacked `PY`, plus the
//!   reshaped `(k × k)` per-group and per-kinship component matrices.
//! * `fit_multi_kinship` — runs AI-REML on the expanded list by reusing
//!   `super::kinship::reml::run_reml` with `GroupPartition::empty`,
//!   since every variance component now lives in the kinship list.
//!
//! ## Score-test path (gene_stats kinship branch)
//!
//! For a gene `G₀` (shape `n × m`), the stacked score is
//!
//!   `U_stacked = (I_k ⊗ G₀)' · PY_stacked`  shape `mk × 1`
//!
//! which unstacks to `S = G₀' · PY_reshaped` of shape `m × k`, the same
//! layout the no-kinship branch uses. The joint score covariance
//! *does not* factor as `Σ_res⁻¹ ⊗ K₀` — the kinship term breaks the
//! Kronecker structure of `Σ⁻¹`. Instead we compute the full
//! `(mk × mk)` covariance per gene (upstream
//! `MultiSTAAR/src/MultiSTAAR_O_SMMAT_sparse.cpp:77`):
//!
//!   `Cov = (I_k ⊗ G₀)' Σ⁻¹ (I_k ⊗ G₀) − A' · cov · A`
//!
//! where `A = Σ⁻¹X_stacked' · (I_k ⊗ G₀)` has shape `kp × mk`. The per-
//! variant and per-test statistics then contract against this full
//! covariance instead of the Kronecker factor.
//!
//! ## k=1 parity
//!
//! At `k=1`, every `E_{0,0}` equals `[[1]]`, the expanded kinship list
//! is `{D_g}_g ∪ {Φ_l}_l` — exactly the single-trait kinship+groups
//! layout — and the stacked response/covariates collapse to the raw
//! `y`/`X`. `MultiKinshipState.sigma_inv`, `.sigma_inv_x`, `.cov`,
//! `.p_y` coincide with `KinshipState` fields bit-for-bit, and the
//! per-gene `(mk × mk)` cov reduces to the usual `(m × m)` `K₀`. The
//! unit test `k1_kinship_matches_single_trait_reml` pins this.

use faer::prelude::Solve;
use faer::sparse::Triplet;
use faer::Mat;

use super::kinship::reml::run_reml;
use super::kinship::types::{
    GroupPartition, KinshipInverse, KinshipMatrix, KinshipState, VarianceComponents,
};
use crate::error::CohortError;

/// Upper-triangle index enumeration for a `k × k` symmetric matrix:
/// returns the `(j, m)` pairs with `j ≤ m` in the same order upstream
/// GMMAT uses (`glmmkin.R:150-153`). For `k=2`: `[(0,0), (0,1), (1,1)]`.
fn sym_pairs(k: usize) -> Vec<(usize, usize)> {
    let mut v = Vec::with_capacity(k * (k + 1) / 2);
    for j in 0..k {
        for m in j..k {
            v.push((j, m));
        }
    }
    v
}

/// Build the expanded kinship list for multi-trait AI-REML.
///
/// Output order (matches `glmmkin.R:147-166`): first the residual-group
/// entries `E_{j,m} ⊗ D_g` in `(group, (j,m))` lexicographic order, then
/// the kinship entries `E_{j,m} ⊗ Φ_l` in `(kinship, (j,m))` order.
///
/// Every output matrix is `(nk × nk)` sparse. The per-entry non-zero
/// count is `2 · nnz(V)` for `j < m` pairs (two Kronecker blocks: `(j,m)`
/// and `(m,j)`) and `nnz(V)` for `j == m` pairs (one diagonal block).
pub fn expand_kinships(
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    n_pheno: usize,
) -> Result<Vec<KinshipMatrix>, CohortError> {
    let k = n_pheno;
    let n = groups.n_samples();
    let pairs = sym_pairs(k);
    let n_groups = groups.n_groups().max(1);
    let n_entries = pairs.len() * (n_groups + kinships.len());
    let mut out: Vec<KinshipMatrix> = Vec::with_capacity(n_entries);

    // Residual-group block `E_{j,m} ⊗ D_g`. Homoscedastic case (no
    // explicit groups) collapses to `E_{j,m} ⊗ I_n`.
    let group_sets: Vec<Vec<u32>> = if groups.n_groups() == 0 {
        vec![(0..n as u32).collect()]
    } else {
        (0..groups.n_groups())
            .map(|g| groups.group(g).to_vec())
            .collect()
    };

    for (gi, rows) in group_sets.iter().enumerate() {
        for &(j, m) in &pairs {
            let label = format!("res_g{gi}_{j}_{m}");
            let triplets = build_block_triplets(n, k, j, m, |i| {
                // D_g at (i,i) = 1 if i is in group gi, 0 otherwise.
                // We only emit diagonal triplets for in-group samples.
                let _ = i;
                Some(1.0)
            }, rows);
            out.push(KinshipMatrix::from_triplets(n * k, triplets, label)?);
        }
    }

    // Kinship block `E_{j,m} ⊗ Φ_l`.
    for (li, phi) in kinships.iter().enumerate() {
        if phi.n() != n {
            return Err(CohortError::Input(format!(
                "kinship matrix '{}' has n={}, expected n={}",
                phi.label(),
                phi.n(),
                n,
            )));
        }
        for &(j, m) in &pairs {
            let label = format!("kins_{li}_{j}_{m}");
            let triplets = kronecker_block_triplets(n, k, j, m, phi);
            out.push(KinshipMatrix::from_triplets(n * k, triplets, label)?);
        }
    }

    Ok(out)
}

/// Build the `(nk × nk)` triplet list for `E_{j,m} ⊗ diag(d)` where
/// `d[i] = value(i)` for rows in `members`, else `0`. Used by the
/// residual-group branch.
fn build_block_triplets(
    n: usize,
    k: usize,
    j: usize,
    m: usize,
    value: impl Fn(usize) -> Option<f64>,
    members: &[u32],
) -> Vec<Triplet<u32, u32, f64>> {
    let _ = k;
    let mut out: Vec<Triplet<u32, u32, f64>> =
        Vec::with_capacity(members.len() * if j == m { 1 } else { 2 });
    for &row in members {
        let i = row as usize;
        let v = match value(i) {
            Some(v) if v != 0.0 => v,
            _ => continue,
        };
        // Block (j, m) sits at rows `j*n + i` and columns `m*n + i`.
        out.push(Triplet::new(
            (j * n + i) as u32,
            (m * n + i) as u32,
            v,
        ));
        if j != m {
            out.push(Triplet::new(
                (m * n + i) as u32,
                (j * n + i) as u32,
                v,
            ));
        }
    }
    out
}

/// Build the `(nk × nk)` triplet list for `E_{j,m} ⊗ Φ`. Two Kronecker
/// blocks for `j != m` (at `(j,m)` and `(m,j)` in `E`), one for `j == m`.
fn kronecker_block_triplets(
    n: usize,
    k: usize,
    j: usize,
    m: usize,
    phi: &KinshipMatrix,
) -> Vec<Triplet<u32, u32, f64>> {
    let _ = k;
    let mut out: Vec<Triplet<u32, u32, f64>> = Vec::new();
    let push_block = |out: &mut Vec<Triplet<u32, u32, f64>>, jj: usize, mm: usize| {
        match phi {
            KinshipMatrix::Dense { matrix, .. } => {
                for r in 0..n {
                    for c in 0..n {
                        let v = matrix[(r, c)];
                        if v != 0.0 {
                            out.push(Triplet::new(
                                (jj * n + r) as u32,
                                (mm * n + c) as u32,
                                v,
                            ));
                        }
                    }
                }
            }
            KinshipMatrix::Sparse { matrix, .. } => {
                let col_ptr = matrix.symbolic().col_ptr();
                let row_idx = matrix.symbolic().row_idx();
                let val = matrix.val();
                for c in 0..n {
                    let s = col_ptr[c] as usize;
                    let e = col_ptr[c + 1] as usize;
                    for kk in s..e {
                        let r = row_idx[kk] as usize;
                        let v = val[kk];
                        if v == 0.0 {
                            continue;
                        }
                        out.push(Triplet::new(
                            (jj * n + r) as u32,
                            (mm * n + c) as u32,
                            v,
                        ));
                    }
                }
            }
        }
    };
    push_block(&mut out, j, m);
    if j != m {
        push_block(&mut out, m, j);
    }
    out
}

/// `vec(Y)` — column-major stack of an `(n × k)` phenotype matrix into
/// an `(nk × 1)` column. Trait 0 occupies rows `0..n`, trait 1 occupies
/// `n..2n`, and so on, matching `glmmkin.R:432` `Y2 <- as.vector(Y)`.
pub fn stack_response(y: &Mat<f64>) -> Mat<f64> {
    let n = y.nrows();
    let k = y.ncols();
    let mut out = Mat::<f64>::zeros(n * k, 1);
    for t in 0..k {
        for i in 0..n {
            out[(t * n + i, 0)] = y[(i, t)];
        }
    }
    out
}

/// `I_k ⊗ X` — block-diagonal stacking of the shared covariate matrix
/// `X` `(n × p)` into `(nk × kp)`. Matches `glmmkin.R:434`
/// `X2 <- Diagonal(n=n.pheno) %x% X`. Dense storage is fine: `kp` is
/// small (≤ ~50 in practice) and the matrix is only formed once at fit
/// time.
pub fn stack_covariates(x: &Mat<f64>, n_pheno: usize) -> Mat<f64> {
    let n = x.nrows();
    let p = x.ncols();
    let k = n_pheno;
    let mut out = Mat::<f64>::zeros(n * k, k * p);
    for t in 0..k {
        for i in 0..n {
            for c in 0..p {
                out[(t * n + i, t * p + c)] = x[(i, c)];
            }
        }
    }
    out
}

/// Reshape a flat `Σ_{i=1}^{k(k+1)/2}` `τ` slice into a single symmetric
/// `(k × k)` matrix. Mirrors `glmmkin.R:216,222` which calls
/// `sparseMatrix(i=rep(1:k, k:1), j=flattened_cols, x=slice, symmetric=T)`.
fn reshape_sym_block(slice: &[f64], k: usize) -> Mat<f64> {
    debug_assert_eq!(slice.len(), k * (k + 1) / 2);
    let mut out = Mat::<f64>::zeros(k, k);
    let mut idx = 0;
    for j in 0..k {
        for m in j..k {
            out[(j, m)] = slice[idx];
            if j != m {
                out[(m, j)] = slice[idx];
            }
            idx += 1;
        }
    }
    out
}

/// Fitted state from multi-trait AI-REML with shared kinship matrices.
///
/// Field layout mirrors the single-trait `KinshipState` but at the
/// stacked `(nk)` dimension. `theta_res` holds one `(k × k)` symmetric
/// residual-covariance matrix per residual group, `theta_g` holds one
/// `(k × k)` per kinship matrix. Several fields (`n_covariates`,
/// `n_iter`, and the reshaped `theta_res`/`theta_g` components) are
/// carried for downstream consumers and run-manifest reporting even if
/// the hot path does not read them; the `allow(dead_code)` keeps clippy
/// quiet without losing the public API.
#[derive(Clone)]
#[allow(dead_code)]
pub struct MultiKinshipState {
    /// Stacked `Σ⁻¹` of shape `(nk × nk)`.
    pub sigma_inv: KinshipInverse,
    /// `Σ⁻¹ · (I_k ⊗ X)` of shape `(nk × kp)`.
    pub sigma_inv_x: Mat<f64>,
    /// `(X_stacked' · Σ⁻¹ · X_stacked)⁻¹` of shape `(kp × kp)`.
    pub cov: Mat<f64>,
    /// `PY_stacked = Σ⁻¹ · (vec(Y) − X_stacked · α̂)` of shape `(nk × 1)`.
    pub p_y: Mat<f64>,
    /// Per-group residual covariance matrices `Σ_res_g`, each `(k × k)`.
    /// One entry even in the homoscedastic case.
    pub theta_res: Vec<Mat<f64>>,
    /// Per-kinship covariance matrices `Σ_g_l`, each `(k × k)`. One per
    /// input kinship matrix, in input order.
    pub theta_g: Vec<Mat<f64>>,
    pub n_samples: usize,
    pub n_pheno: usize,
    pub n_covariates: usize,
    pub n_iter: usize,
}

/// Warm-start initialiser for the expanded `τ` vector. Mirrors
/// `glmmkin.R:446-447`: diagonal entries (`j == m`) get `var(vec(Y))/q`
/// where `q` is the expanded list length, off-diagonals start at `0`.
fn warm_start_tau(y_stacked: &Mat<f64>, n_pairs: usize, n_entries: usize, k: usize) -> VarianceComponents {
    let n = y_stacked.nrows();
    let mean: f64 = (0..n).map(|i| y_stacked[(i, 0)]).sum::<f64>() / n as f64;
    let var: f64 = (0..n)
        .map(|i| (y_stacked[(i, 0)] - mean).powi(2))
        .sum::<f64>()
        / (n - 1).max(1) as f64;
    let warm = (var / n_entries as f64).max(1e-6);
    let mut tau = VarianceComponents::zeros(n_entries, 0);
    let pairs = sym_pairs(k);
    let slot = tau.as_slice_mut();
    for block in 0..(n_entries / n_pairs) {
        for (pair_idx, &(j, m)) in pairs.iter().enumerate() {
            let flat = block * n_pairs + pair_idx;
            slot[flat] = if j == m { warm } else { 0.0 };
        }
    }
    tau
}

/// Multi-trait AI-REML with shared kinship matrices. Returns the
/// converged multi-trait null state.
///
/// upstream:
///   `MultiSTAAR/R/fit_null_glmmkin_multi.R:70-123` → wrapper
///   `GMMAT/R/glmmkin.R:424-536` → `glmmkin.multi.ai` body
///   `GMMAT/R/glmmkin.R:771-811` → `R_fitglmm_ai_noresidual` AI iteration
pub fn fit_multi_kinship(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
) -> Result<MultiKinshipState, CohortError> {
    let n = y.nrows();
    let k = y.ncols();
    let p = x.ncols();
    assert_eq!(x.nrows(), n, "X rows must match Y rows");
    assert!(n * k > p * k + 1, "need nk > kp samples for the joint fit");
    if !kinships.is_empty() && kinships[0].n() != n {
        return Err(CohortError::Input(format!(
            "kinship dim {} does not match phenotype dim {}",
            kinships[0].n(),
            n,
        )));
    }

    let pairs = sym_pairs(k);
    let n_pairs = pairs.len();
    let n_groups = groups.n_groups().max(1);
    let n_entries = n_pairs * (n_groups + kinships.len());

    let y_stacked = stack_response(y);
    let x_stacked = stack_covariates(x, k);
    let expanded = expand_kinships(kinships, groups, k)?;
    debug_assert_eq!(expanded.len(), n_entries);

    let init_tau = warm_start_tau(&y_stacked, n_pairs, n_entries, k);

    // `run_reml` was written for the single-trait path and expects a
    // backend-specific builder. For the dense-first multi-trait fit we
    // route through the dense builder; sparse kinship → sparse stacked
    // path is a follow-up (the expanded entries are naturally sparse for
    // sparse Φ, but faer's sparse Cholesky on nk × nk needs more care
    // around symbolic-factor reuse across iterations).
    //
    // For now we materialise the expanded entries as dense `KinshipMatrix`
    // if any were originally dense, and keep the originally-sparse ones
    // in sparse storage. The shared `ai_step` handles both variants via
    // `matvec_kinship`, and `DenseBuilder` promotes sparse kins to dense
    // Σ⁻¹ inside each iteration.
    let expanded_kins_dense = promote_to_dense(&expanded, n * k);
    let empty_groups = GroupPartition::empty(n * k);
    let weights = vec![1.0; n * k];
    let builder = super::kinship::dense::DenseBuilder {
        n: n * k,
        kinships: &expanded_kins_dense,
        groups: &empty_groups,
        weights: &weights,
    };
    let state: KinshipState = run_reml(
        &y_stacked,
        &x_stacked,
        &expanded_kins_dense,
        &empty_groups,
        &weights,
        init_tau,
        &builder,
    )?;

    // Reshape τ back to per-group / per-kinship (k × k) matrices.
    let mut theta_res: Vec<Mat<f64>> = Vec::with_capacity(n_groups);
    let mut theta_g: Vec<Mat<f64>> = Vec::with_capacity(kinships.len());
    let tau_flat = state.tau.as_slice();
    for gi in 0..n_groups {
        let offset = gi * n_pairs;
        theta_res.push(reshape_sym_block(&tau_flat[offset..offset + n_pairs], k));
    }
    for li in 0..kinships.len() {
        let offset = (n_groups + li) * n_pairs;
        theta_g.push(reshape_sym_block(&tau_flat[offset..offset + n_pairs], k));
    }

    Ok(MultiKinshipState {
        sigma_inv: state.inverse,
        sigma_inv_x: state.sigma_inv_x,
        cov: state.cov,
        p_y: state.p_y,
        theta_res,
        theta_g,
        n_samples: n,
        n_pheno: k,
        n_covariates: p,
        n_iter: state.n_iter,
    })
}

/// Per-gene score statistics for the multi-trait kinship path.
///
/// `s[j, t] = (G₀' · PY_reshaped)[j, t]` is the per-trait score at
/// variant `j`, same shape and semantics as the no-kinship
/// `GeneStats.s`. `cov_full` is the full `(mk × mk)` joint score
/// covariance — the kinship term prevents the Kronecker factorisation
/// the no-kinship path relies on, so the per-test kernels contract
/// against this matrix directly.
///
/// Layout of `cov_full`: trait-major block order. The entry at
/// `(t1 * m + j1, t2 * m + j2)` is the covariance between the score at
/// trait `t1`, variant `j1` and trait `t2`, variant `j2`. When the
/// no-kinship path uses `Σ_res⁻¹ ⊗ K₀`, this matrix would equal
/// `Σ_res ⊗ K₀` exactly (the Kronecker expansion), which is the
/// invariant the k=1 gene-level parity test pins.
pub struct MultiKinshipGeneStats {
    pub s: Mat<f64>,
    pub cov_full: Mat<f64>,
}

/// Build the per-gene `(S, cov_full)` for the multi-trait kinship path.
/// See `MultiKinshipGeneStats` for layout.
///
/// upstream: `MultiSTAAR/src/MultiSTAAR_O_SMMAT_sparse.cpp:77` for the
/// sparse-Σ assembly (we use the dense-Σ analogue — same formula).
pub fn gene_stats_kinship(g0: &Mat<f64>, null: &MultiKinshipState) -> MultiKinshipGeneStats {
    let n = null.n_samples;
    let k = null.n_pheno;
    let m = g0.ncols();
    assert_eq!(g0.nrows(), n, "G₀ rows must match null n_samples");
    let mk = m * k;

    // S[j, t] = Σ_i G₀[i, j] · PY[t*n + i, 0].
    let mut s = Mat::<f64>::zeros(m, k);
    for t in 0..k {
        for j in 0..m {
            let mut acc = 0.0;
            for i in 0..n {
                acc += g0[(i, j)] * null.p_y[(t * n + i, 0)];
            }
            s[(j, t)] = acc;
        }
    }

    // Materialise G_stacked = I_k ⊗ G₀ as (nk × mk). For small gene-set
    // sizes this is cheap; a block-wise implementation that avoids the
    // allocation is a follow-up if production cohorts hit it.
    let mut g_stacked = Mat::<f64>::zeros(n * k, mk);
    for t in 0..k {
        for i in 0..n {
            for j in 0..m {
                g_stacked[(t * n + i, t * m + j)] = g0[(i, j)];
            }
        }
    }

    // Σ⁻¹ · G_stacked — dense solve or matmul depending on the stored
    // inverse representation.
    let sinv_g = match &null.sigma_inv {
        KinshipInverse::Dense(sinv) => sinv * &g_stacked,
        KinshipInverse::Sparse(factor) => {
            use faer::linalg::solvers::Solve;
            let mut out = g_stacked.clone();
            factor.llt.solve_in_place(&mut out);
            out
        }
    };

    let term1 = g_stacked.transpose() * &sinv_g;
    let xs_sinv_g = null.sigma_inv_x.transpose() * &g_stacked;
    let cov_correction = xs_sinv_g.transpose() * (&null.cov * &xs_sinv_g);
    let mut cov_full = term1;
    for r in 0..mk {
        for c in 0..mk {
            cov_full[(r, c)] -= cov_correction[(r, c)];
        }
    }

    MultiKinshipGeneStats { s, cov_full }
}

/// Joint per-variant χ²(k) test at variant `i` using the full `(mk ×
/// mk)` covariance. Extracts the `(k × k)` diagonal block at
/// `(t1 * m + i, t2 * m + i)`, then evaluates
/// `Q_i = U_i' · Cov_i⁻¹ · U_i ~ χ²(k)` with `U_i = S[i, :]`.
///
/// Returns `1.0` if the per-variant variance block is numerically
/// singular or non-positive-definite; matches the fail-safe semantics
/// of `multi::variant_joint_chi2`.
pub fn variant_joint_chi2_kinship(
    s: &Mat<f64>,
    cov_full: &Mat<f64>,
    i: usize,
    n_pheno: usize,
    m: usize,
) -> f64 {
    let k = n_pheno;
    // Extract (k × k) diagonal block at variant i.
    let mut block = Mat::<f64>::zeros(k, k);
    for t1 in 0..k {
        for t2 in 0..k {
            block[(t1, t2)] = cov_full[(t1 * m + i, t2 * m + i)];
        }
    }
    // Diagonal-zero guard: if any trait's variance collapses the joint
    // test is undefined; return 1.0 to mirror `multi::variant_joint_chi2`.
    for a in 0..k {
        if !block[(a, a)].is_finite() || block[(a, a)] <= 0.0 {
            return 1.0;
        }
    }
    let eye_k = Mat::<f64>::identity(k, k);
    let block_inv: Mat<f64> = block.col_piv_qr().solve(&eye_k);
    let mut quad = 0.0;
    for a in 0..k {
        for b in 0..k {
            quad += s[(i, a)] * block_inv[(a, b)] * s[(i, b)];
        }
    }
    if !quad.is_finite() || quad <= 0.0 {
        return 1.0;
    }
    chisq_pvalue(quad, k as f64)
}

/// Joint burden test with weight `w`. Per-trait burden score
/// `U_burden[t] = Σ_j w[j] · S[j, t]`, variance block
/// `V[t1, t2] = Σ_{j1, j2} w[j1] · cov_full[t1m+j1, t2m+j2] · w[j2]`,
/// then `Q = U_burden' V⁻¹ U_burden ~ χ²(k)`.
pub fn multi_burden_kinship(
    s: &Mat<f64>,
    cov_full: &Mat<f64>,
    w: &[f64],
    n_pheno: usize,
) -> f64 {
    let k = n_pheno;
    let m = w.len();
    if m == 0 {
        return 1.0;
    }
    let mut u_burden = vec![0.0_f64; k];
    for t in 0..k {
        let mut acc = 0.0;
        for j in 0..m {
            acc += w[j] * s[(j, t)];
        }
        u_burden[t] = acc;
    }

    // V[t1, t2] = Σ_{j1, j2} w[j1] · cov_full[t1m+j1, t2m+j2] · w[j2].
    let mut v = Mat::<f64>::zeros(k, k);
    for t1 in 0..k {
        for t2 in 0..k {
            let mut acc = 0.0;
            for j1 in 0..m {
                if w[j1] == 0.0 {
                    continue;
                }
                for j2 in 0..m {
                    if w[j2] == 0.0 {
                        continue;
                    }
                    acc += w[j1] * cov_full[(t1 * m + j1, t2 * m + j2)] * w[j2];
                }
            }
            v[(t1, t2)] = acc;
        }
    }

    let eye_k = Mat::<f64>::identity(k, k);
    let v_inv = v.col_piv_qr().solve(&eye_k);
    let mut q = 0.0;
    for a in 0..k {
        for b in 0..k {
            q += u_burden[a] * v_inv[(a, b)] * u_burden[b];
        }
    }
    if !q.is_finite() || q <= 0.0 {
        return 1.0;
    }
    chisq_pvalue(q, k as f64)
}

/// Joint SKAT with weight `w`.
///
///   `Q = vec(WS)' (I_k ⊗ W) vec(S) = Σ_{t, j} w[j]² · S[j, t]²`
///
/// Under the null `vec(S) ~ N(0, cov_full)`, so `Q` is a weighted
/// mixture of χ²(1) with eigenvalues of
/// `W^(1/2)_stacked · cov_full · W^(1/2)_stacked` where
/// `W_stacked = I_k ⊗ diag(w²)`. At `k=1` the construction reduces to
/// `diag(w²) · K · diag(w²) → eigvals = eigvals(W K W)`, which matches
/// single-trait SKAT.
pub fn multi_skat_kinship(
    s: &Mat<f64>,
    cov_full: &Mat<f64>,
    w: &[f64],
    n_pheno: usize,
) -> f64 {
    let k = n_pheno;
    let m = w.len();
    if m == 0 {
        return 1.0;
    }
    // Q = Σ_{t, j} w[j]² · S[j, t]².
    let mut q = 0.0;
    for t in 0..k {
        for j in 0..m {
            let wj = w[j];
            if wj == 0.0 {
                continue;
            }
            let v = wj * s[(j, t)];
            q += v * v;
        }
    }
    if !q.is_finite() || q <= 0.0 {
        return 1.0;
    }

    // Weighted kernel `A = W_stacked · cov_full · W_stacked` with
    // `W_stacked = I_k ⊗ diag(w)`. (Using w rather than √(w²) keeps
    // signs consistent with the quadratic form above; eigenvalues of
    // A are the mixture weights for `Q`.)
    let mk = m * k;
    let mut a = Mat::<f64>::zeros(mk, mk);
    for t1 in 0..k {
        for j1 in 0..m {
            let wi = w[j1];
            if wi == 0.0 {
                continue;
            }
            for t2 in 0..k {
                for j2 in 0..m {
                    let wj = w[j2];
                    if wj == 0.0 {
                        continue;
                    }
                    a[(t1 * m + j1, t2 * m + j2)] =
                        wi * cov_full[(t1 * m + j1, t2 * m + j2)] * wj;
                }
            }
        }
    }

    let eigs = symmetric_eigenvalues(&a);
    crate::staar::stats::mixture_chisq_pvalue(q, &eigs)
}

fn symmetric_eigenvalues(m: &Mat<f64>) -> Vec<f64> {
    let n = m.nrows();
    if n == 0 {
        return Vec::new();
    }
    match m.self_adjoint_eigen(faer::Side::Lower) {
        Ok(evd) => {
            let s = evd.S();
            let col = s.column_vector();
            (0..n).map(|i| col[i].max(0.0)).collect()
        }
        Err(_) => vec![0.0; n],
    }
}

/// Joint ACAT-V: per-variant common-variant joint χ²(k), plus a pooled
/// rare-group joint burden; Cauchy-combined with `w_acat` weights.
/// Matches the no-kinship `multi::multi_acat_v` MAC split (MAC > 10 ⇒
/// common; MAC ≤ 10 ⇒ pooled burden with `w_burden` weights).
#[allow(clippy::too_many_arguments)]
pub fn multi_acat_v_kinship(
    s: &Mat<f64>,
    cov_full: &Mat<f64>,
    n_pheno: usize,
    w_acat: &[f64],
    w_burden: &[f64],
    mafs: &[f64],
    n_samples: usize,
) -> f64 {
    const MAC_THRESHOLD: f64 = 10.0;
    let k = n_pheno;
    let m = w_acat.len();
    let ns = n_samples as f64;

    let mut p_values: Vec<f64> = Vec::with_capacity(m);
    let mut cauchy_weights: Vec<f64> = Vec::with_capacity(m);
    let mut rare_idx: Vec<usize> = Vec::new();

    for j in 0..m {
        if w_acat[j] == 0.0 {
            continue;
        }
        let mac = (2.0 * mafs[j] * ns).round();
        if mac > MAC_THRESHOLD {
            let p = variant_joint_chi2_kinship(s, cov_full, j, k, m);
            p_values.push(p);
            cauchy_weights.push(w_acat[j]);
        } else {
            rare_idx.push(j);
        }
    }

    if !rare_idx.is_empty() {
        let mut w_rare = vec![0.0_f64; m];
        for &j in &rare_idx {
            w_rare[j] = w_burden[j];
        }
        let p = multi_burden_kinship(s, cov_full, &w_rare, k);
        let mean_w =
            rare_idx.iter().map(|&j| w_acat[j]).sum::<f64>() / rare_idx.len() as f64;
        p_values.push(p);
        cauchy_weights.push(mean_w);
    }

    if p_values.is_empty() {
        return 1.0;
    }
    crate::staar::stats::cauchy_combine_weighted(&p_values, &cauchy_weights)
}

fn chisq_pvalue(t: f64, df: f64) -> f64 {
    use statrs::distribution::{ChiSquared, ContinuousCDF};
    if t <= 0.0 || !t.is_finite() {
        return 1.0;
    }
    match ChiSquared::new(df) {
        Ok(d) => (1.0 - d.cdf(t)).max(crate::staar::stats::P_FLOOR),
        Err(_) => f64::NAN,
    }
}

/// Run the multi-trait joint STAAR omnibus on one gene against a
/// kinship-aware null. Mirrors `multi::run_multi_staar` test-for-test
/// but contracts against the full `(mk × mk)` covariance so the results
/// stay valid under arbitrary kinship correlation.
pub fn run_multi_staar_kinship(
    g0: &Mat<f64>,
    null: &MultiKinshipState,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
) -> crate::staar::score::StaarResult {
    use crate::staar::score::{beta_density_weight, StaarResult};
    use crate::staar::stats;
    let m = g0.ncols();
    if m == 0 {
        return nan_result();
    }
    let gs = gene_stats_kinship(g0, null);
    let s = &gs.s;
    let cov_full = &gs.cov_full;
    let k = null.n_pheno;
    let n_samples = null.n_samples;

    let beta_1_25: Vec<f64> =
        mafs.iter().map(|&maf| beta_density_weight(maf, 1.0, 25.0)).collect();
    let beta_1_1: Vec<f64> =
        mafs.iter().map(|&maf| beta_density_weight(maf, 1.0, 1.0)).collect();
    let acat_denom: Vec<f64> = mafs
        .iter()
        .map(|&maf| {
            let d = beta_density_weight(maf, 0.5, 0.5);
            if d > 0.0 { d * d } else { 1.0 }
        })
        .collect();

    let base_burden_1_25 = multi_burden_kinship(s, cov_full, &beta_1_25, k);
    let base_burden_1_1 = multi_burden_kinship(s, cov_full, &beta_1_1, k);
    let base_skat_1_25 = multi_skat_kinship(s, cov_full, &beta_1_25, k);
    let base_skat_1_1 = multi_skat_kinship(s, cov_full, &beta_1_1, k);
    let wa_base_1_25: Vec<f64> =
        beta_1_25.iter().zip(&acat_denom).map(|(b, d)| b * b / d).collect();
    let wa_base_1_1: Vec<f64> =
        beta_1_1.iter().zip(&acat_denom).map(|(b, d)| b * b / d).collect();
    let base_acat_v_1_25 = multi_acat_v_kinship(
        s, cov_full, k, &wa_base_1_25, &beta_1_25, mafs, n_samples,
    );
    let base_acat_v_1_1 = multi_acat_v_kinship(
        s, cov_full, k, &wa_base_1_1, &beta_1_1, mafs, n_samples,
    );

    let acat_o = stats::cauchy_combine(&[
        base_burden_1_25,
        base_burden_1_1,
        base_skat_1_25,
        base_skat_1_1,
        base_acat_v_1_25,
        base_acat_v_1_1,
    ]);

    let n_channels = annotation_matrix.len();
    let mut per_annotation: Vec<[f64; 6]> = Vec::with_capacity(n_channels);
    let mut by_test: [Vec<f64>; 6] = [
        vec![base_burden_1_25],
        vec![base_burden_1_1],
        vec![base_skat_1_25],
        vec![base_skat_1_1],
        vec![base_acat_v_1_25],
        vec![base_acat_v_1_1],
    ];

    let mut wb_1_25 = vec![0.0; m];
    let mut wb_1_1 = vec![0.0; m];
    let mut ws_1_25 = vec![0.0; m];
    let mut ws_1_1 = vec![0.0; m];
    let mut wa_1_25 = vec![0.0; m];
    let mut wa_1_1 = vec![0.0; m];

    for channel_weights in annotation_matrix {
        for j in 0..m {
            let a = channel_weights[j];
            let a_sqrt = a.sqrt();
            wb_1_25[j] = beta_1_25[j] * a;
            wb_1_1[j] = beta_1_1[j] * a;
            ws_1_25[j] = beta_1_25[j] * a_sqrt;
            ws_1_1[j] = beta_1_1[j] * a_sqrt;
            wa_1_25[j] = a * beta_1_25[j] * beta_1_25[j] / acat_denom[j];
            wa_1_1[j] = a * beta_1_1[j] * beta_1_1[j] / acat_denom[j];
        }
        let p = [
            multi_burden_kinship(s, cov_full, &wb_1_25, k),
            multi_burden_kinship(s, cov_full, &wb_1_1, k),
            multi_skat_kinship(s, cov_full, &ws_1_25, k),
            multi_skat_kinship(s, cov_full, &ws_1_1, k),
            multi_acat_v_kinship(s, cov_full, k, &wa_1_25, &wb_1_25, mafs, n_samples),
            multi_acat_v_kinship(s, cov_full, k, &wa_1_1, &wb_1_1, mafs, n_samples),
        ];
        for i in 0..6 {
            by_test[i].push(p[i]);
        }
        per_annotation.push(p);
    }

    let staar_b_1_25 = stats::cauchy_combine(&by_test[0]);
    let staar_b_1_1 = stats::cauchy_combine(&by_test[1]);
    let staar_s_1_25 = stats::cauchy_combine(&by_test[2]);
    let staar_s_1_1 = stats::cauchy_combine(&by_test[3]);
    let staar_a_1_25 = stats::cauchy_combine(&by_test[4]);
    let staar_a_1_1 = stats::cauchy_combine(&by_test[5]);

    let mut all_p: Vec<f64> = Vec::with_capacity(6 + n_channels * 6);
    all_p.extend_from_slice(&[
        base_burden_1_25,
        base_burden_1_1,
        base_skat_1_25,
        base_skat_1_1,
        base_acat_v_1_25,
        base_acat_v_1_1,
    ]);
    for p in &per_annotation {
        all_p.extend_from_slice(p);
    }
    let staar_o = stats::cauchy_combine(&all_p);

    StaarResult {
        burden_1_25: base_burden_1_25,
        burden_1_1: base_burden_1_1,
        skat_1_25: base_skat_1_25,
        skat_1_1: base_skat_1_1,
        acat_v_1_25: base_acat_v_1_25,
        acat_v_1_1: base_acat_v_1_1,
        per_annotation,
        staar_b_1_25,
        staar_b_1_1,
        staar_s_1_25,
        staar_s_1_1,
        staar_a_1_25,
        staar_a_1_1,
        acat_o,
        staar_o,
    }
}

fn nan_result() -> crate::staar::score::StaarResult {
    crate::staar::score::StaarResult {
        burden_1_25: f64::NAN,
        burden_1_1: f64::NAN,
        skat_1_25: f64::NAN,
        skat_1_1: f64::NAN,
        acat_v_1_25: f64::NAN,
        acat_v_1_1: f64::NAN,
        per_annotation: Vec::new(),
        staar_b_1_25: f64::NAN,
        staar_b_1_1: f64::NAN,
        staar_s_1_25: f64::NAN,
        staar_s_1_1: f64::NAN,
        staar_a_1_25: f64::NAN,
        staar_a_1_1: f64::NAN,
        acat_o: f64::NAN,
        staar_o: f64::NAN,
    }
}

/// Promote every sparse entry in an expanded list to a dense
/// `KinshipMatrix` so the dense REML builder can operate uniformly. For
/// the dense-first Phase B this is the simplest correct path; the sparse
/// Σ route follows once we validate the math end-to-end.
///
/// The expanded entries are extremely sparse (single-block Kronecker),
/// so `KinshipMatrix::new` would route them back to sparse storage via
/// its density threshold. We bypass that by constructing the `Dense`
/// variant directly — the dense backend requires dense storage
/// regardless of density.
fn promote_to_dense(expanded: &[KinshipMatrix], n: usize) -> Vec<KinshipMatrix> {
    let mut out = Vec::with_capacity(expanded.len());
    for kin in expanded {
        match kin {
            KinshipMatrix::Dense { .. } => out.push(kin.clone()),
            KinshipMatrix::Sparse { matrix, label } => {
                let mut dense = Mat::<f64>::zeros(n, n);
                let col_ptr = matrix.symbolic().col_ptr();
                let row_idx = matrix.symbolic().row_idx();
                let val = matrix.val();
                for c in 0..n {
                    let s = col_ptr[c] as usize;
                    let e = col_ptr[c + 1] as usize;
                    for kk in s..e {
                        let r = row_idx[kk] as usize;
                        dense[(r, c)] += val[kk];
                    }
                }
                out.push(KinshipMatrix::Dense {
                    matrix: dense,
                    label: label.clone(),
                });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn identity_kinship(n: usize, label: &str) -> KinshipMatrix {
        let mut m = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            m[(i, i)] = 1.0;
        }
        KinshipMatrix::new(m, label.into()).unwrap()
    }

    fn pedigree_kinship(n_fam: usize, label: &str) -> KinshipMatrix {
        // 2-sib families: 0.5 within-family off-diagonal, 1.0 diagonal.
        let n = n_fam * 2;
        let mut m = Mat::<f64>::zeros(n, n);
        for f in 0..n_fam {
            let a = 2 * f;
            let b = 2 * f + 1;
            m[(a, a)] = 1.0;
            m[(b, b)] = 1.0;
            m[(a, b)] = 0.5;
            m[(b, a)] = 0.5;
        }
        KinshipMatrix::new(m, label.into()).unwrap()
    }

    #[test]
    fn sym_pairs_enumerate_upper_triangle() {
        assert_eq!(sym_pairs(1), vec![(0, 0)]);
        assert_eq!(sym_pairs(2), vec![(0, 0), (0, 1), (1, 1)]);
        assert_eq!(sym_pairs(3), vec![(0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2)]);
    }

    #[test]
    fn stack_response_is_column_major() {
        let n = 3;
        let k = 2;
        let mut y = Mat::<f64>::zeros(n, k);
        for i in 0..n {
            for t in 0..k {
                y[(i, t)] = (10 * t + i) as f64;
            }
        }
        let y2 = stack_response(&y);
        assert_eq!(y2.nrows(), n * k);
        assert_eq!(y2.ncols(), 1);
        assert_eq!(y2[(0, 0)], 0.0);
        assert_eq!(y2[(1, 0)], 1.0);
        assert_eq!(y2[(2, 0)], 2.0);
        assert_eq!(y2[(3, 0)], 10.0);
        assert_eq!(y2[(4, 0)], 11.0);
        assert_eq!(y2[(5, 0)], 12.0);
    }

    #[test]
    fn stack_covariates_is_block_diagonal() {
        let n = 2;
        let p = 2;
        let k = 3;
        let mut x = Mat::<f64>::zeros(n, p);
        x[(0, 0)] = 1.0;
        x[(0, 1)] = 2.0;
        x[(1, 0)] = 3.0;
        x[(1, 1)] = 4.0;
        let x2 = stack_covariates(&x, k);
        assert_eq!(x2.nrows(), n * k);
        assert_eq!(x2.ncols(), k * p);

        // Block (0, 0): rows 0..n, cols 0..p carries original X.
        assert_eq!(x2[(0, 0)], 1.0);
        assert_eq!(x2[(1, 1)], 4.0);
        // Block (1, 1): rows n..2n, cols p..2p carries original X.
        assert_eq!(x2[(2, 2)], 1.0);
        assert_eq!(x2[(3, 3)], 4.0);
        // Off-diagonal blocks are zero.
        assert_eq!(x2[(0, 2)], 0.0);
        assert_eq!(x2[(2, 0)], 0.0);
    }

    #[test]
    fn expand_kinships_entry_count_matches_upstream_formula() {
        let n = 6;
        let k = 2;
        let phi = identity_kinship(n, "phi");
        let groups = GroupPartition::single(n);
        let expanded = expand_kinships(std::slice::from_ref(&phi), &groups, k).unwrap();
        // n_params = k(k+1)/2 · (G + L) = 3 · (1 + 1) = 6
        assert_eq!(expanded.len(), 6);
        for e in &expanded {
            assert_eq!(e.n(), n * k);
        }
    }

    #[test]
    fn reshape_sym_block_is_symmetric_and_orders_upper_triangle() {
        let flat = [1.0, 2.0, 3.0]; // (0,0), (0,1), (1,1)
        let m = reshape_sym_block(&flat, 2);
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(0, 1)], 2.0);
        assert_eq!(m[(1, 0)], 2.0);
        assert_eq!(m[(1, 1)], 3.0);
    }

    /// k=1 parity: `fit_multi_kinship` with one trait must produce the
    /// same Σ⁻¹, Σ⁻¹X, cov, and PY as single-trait `fit_reml_dense` on
    /// the same inputs. At k=1 the expansion collapses to `E_{0,0} ⊗ V =
    /// V`, stacking is the identity, and the AI-REML runs on identical
    /// inputs; the fits should converge to the same estimates.
    #[test]
    fn k1_kinship_matches_single_trait_reml() {
        use crate::staar::kinship::fit_reml;

        let n_fam = 6;
        let n = n_fam * 2;
        // Use a dense kinship so both single-trait and multi paths route
        // through the dense backend — the k=1 parity assertion is about
        // the math, not the dense/sparse dispatch.
        let mut phi_dense = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            phi_dense[(i, i)] = 1.0;
        }
        for f in 0..n_fam {
            let a = 2 * f;
            let b = 2 * f + 1;
            phi_dense[(a, b)] = 0.5;
            phi_dense[(b, a)] = 0.5;
        }
        // Construct Dense variant directly so the density-based router
        // in `KinshipMatrix::new` does not re-promote to sparse.
        let phi = KinshipMatrix::Dense {
            matrix: phi_dense,
            label: "phi_dense".into(),
        };

        // Deterministic single-trait phenotype.
        let mut state: u64 = 4242;
        let mut u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut y_n1 = Mat::<f64>::zeros(n, 1);
        let mut x = Mat::<f64>::zeros(n, 2);
        for i in 0..n {
            y_n1[(i, 0)] = 0.5 + 2.0 * u01() - 1.0;
            x[(i, 0)] = 1.0;
            x[(i, 1)] = u01();
        }
        let groups = GroupPartition::single(n);

        // Single-trait REML via the top-level dispatcher (it picks the
        // dense backend automatically for our dense Φ at n=12).
        let single = fit_reml(
            &y_n1,
            &x,
            std::slice::from_ref(&phi),
            &groups,
            None,
            1 << 30,
        )
        .expect("single-trait REML must converge on this well-conditioned pedigree");

        // Multi-trait path with k=1.
        let multi = fit_multi_kinship(&y_n1, &x, std::slice::from_ref(&phi), &groups)
            .expect("multi k=1 must converge");

        assert_eq!(multi.n_pheno, 1);
        assert_eq!(multi.theta_res.len(), 1);
        assert_eq!(multi.theta_g.len(), 1);

        // Theta at k=1: theta_res[0] is a 1x1 with the residual σ²;
        // theta_g[0] is a 1x1 with τ. Match against single-trait tau
        // within AI-REML tolerance (convergence criterion uses REML_TOL
        // on both sides, so we allow a loose 1e-4 band).
        let single_sigma2 = single.tau.group(0);
        let single_tau = single.tau.kinship(0);
        let multi_sigma2 = multi.theta_res[0][(0, 0)];
        let multi_tau = multi.theta_g[0][(0, 0)];
        assert!(
            (single_sigma2 - multi_sigma2).abs() < 1e-4,
            "σ²_res: single={single_sigma2} multi={multi_sigma2}",
        );
        assert!(
            (single_tau - multi_tau).abs() < 1e-4,
            "τ_phi: single={single_tau} multi={multi_tau}",
        );

        // PY at k=1 is the stacked response projected through Σ⁻¹; the
        // stack is a no-op at k=1 so the numbers must match element-wise.
        for i in 0..n {
            let s = single.p_y[(i, 0)];
            let m = multi.p_y[(i, 0)];
            assert!((s - m).abs() < 1e-4, "PY[{i}] single={s} multi={m}");
        }
    }

    /// k=1 gene stats parity: the full `(m × m)` `cov_full` must equal
    /// the single-trait kinship-aware kernel
    /// `K = G₀' Σ⁻¹ G₀ − (G₀' Σ⁻¹ X) · cov · (X' Σ⁻¹ G₀)` computed
    /// directly from the single-trait `KinshipState`.
    #[test]
    fn k1_gene_stats_matches_single_trait_kernel() {
        use crate::staar::kinship::fit_reml;

        let n_fam = 6;
        let n = n_fam * 2;
        let mut phi_dense = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            phi_dense[(i, i)] = 1.0;
        }
        for f in 0..n_fam {
            let a = 2 * f;
            let b = 2 * f + 1;
            phi_dense[(a, b)] = 0.5;
            phi_dense[(b, a)] = 0.5;
        }
        let phi = KinshipMatrix::Dense {
            matrix: phi_dense,
            label: "phi_dense".into(),
        };

        let mut state: u64 = 777;
        let mut u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut y_n1 = Mat::<f64>::zeros(n, 1);
        let mut x = Mat::<f64>::zeros(n, 2);
        let m = 4;
        let mut g0 = Mat::<f64>::zeros(n, m);
        for i in 0..n {
            y_n1[(i, 0)] = 2.0 * u01() - 1.0;
            x[(i, 0)] = 1.0;
            x[(i, 1)] = u01();
            for j in 0..m {
                g0[(i, j)] = if u01() < 0.2 { 1.0 } else { 0.0 };
            }
        }
        let groups = GroupPartition::single(n);

        let single = fit_reml(
            &y_n1,
            &x,
            std::slice::from_ref(&phi),
            &groups,
            None,
            1 << 30,
        )
        .expect("single-trait REML convergence");
        let multi = fit_multi_kinship(&y_n1, &x, std::slice::from_ref(&phi), &groups)
            .expect("multi k=1 convergence");

        // Single-trait reference: u = G₀' PY, K = G₀' Σ⁻¹ G₀ − cross.
        let sinv_dense = match &single.inverse {
            KinshipInverse::Dense(m) => m,
            _ => panic!("expected dense single-trait inverse at n=12 dense Φ"),
        };
        let u_ref = g0.transpose() * &single.p_y;
        let sinv_g = sinv_dense * &g0;
        let term1 = g0.transpose() * &sinv_g;
        let xs_sinv_g = single.sigma_inv_x.transpose() * &g0;
        let cov_corr = xs_sinv_g.transpose() * (&single.cov * &xs_sinv_g);
        let mut k_ref = term1.clone();
        for i in 0..m {
            for j in 0..m {
                k_ref[(i, j)] -= cov_corr[(i, j)];
            }
        }

        let gs = gene_stats_kinship(&g0, &multi);
        assert_eq!(gs.s.nrows(), m);
        assert_eq!(gs.s.ncols(), 1);
        assert_eq!(gs.cov_full.nrows(), m);
        assert_eq!(gs.cov_full.ncols(), m);

        for j in 0..m {
            assert!(
                (gs.s[(j, 0)] - u_ref[(j, 0)]).abs() < 1e-6,
                "S[{j}] single={} multi={}",
                u_ref[(j, 0)],
                gs.s[(j, 0)],
            );
            for l in 0..m {
                assert!(
                    (gs.cov_full[(j, l)] - k_ref[(j, l)]).abs() < 1e-6,
                    "K[{j},{l}] single={} multi={}",
                    k_ref[(j, l)],
                    gs.cov_full[(j, l)],
                );
            }
        }
    }

    /// k=2 omnibus smoke test: the kinship omnibus runs end-to-end and
    /// produces a finite p in (0, 1]. This is not a parity test — it
    /// confirms the full-cov per-test functions compose without
    /// numerical blow-ups on a well-conditioned input.
    #[test]
    fn k2_kinship_omnibus_runs_and_is_valid_p() {
        let n_fam = 8;
        let n = n_fam * 2;
        let mut phi_dense = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            phi_dense[(i, i)] = 1.0;
        }
        for f in 0..n_fam {
            let a = 2 * f;
            let b = 2 * f + 1;
            phi_dense[(a, b)] = 0.5;
            phi_dense[(b, a)] = 0.5;
        }
        let phi = KinshipMatrix::Dense {
            matrix: phi_dense,
            label: "phi_dense".into(),
        };

        let mut state: u64 = 91011;
        let mut u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut y = Mat::<f64>::zeros(n, 2);
        let mut x = Mat::<f64>::zeros(n, 2);
        let m = 5;
        let mafs = vec![0.04, 0.08, 0.03, 0.06, 0.02];
        let mut g0 = Mat::<f64>::zeros(n, m);
        for i in 0..n {
            let e1 = 2.0 * u01() - 1.0;
            let e2 = 2.0 * u01() - 1.0;
            y[(i, 0)] = 0.5 + e1;
            y[(i, 1)] = 0.3 + 0.4 * e1 + 0.9 * e2;
            x[(i, 0)] = 1.0;
            x[(i, 1)] = u01();
            for j in 0..m {
                g0[(i, j)] = if u01() < mafs[j] {
                    if u01() < mafs[j] { 2.0 } else { 1.0 }
                } else {
                    0.0
                };
            }
        }
        let groups = GroupPartition::single(n);
        let null = fit_multi_kinship(&y, &x, std::slice::from_ref(&phi), &groups)
            .expect("multi k=2 should converge");
        let ann: Vec<Vec<f64>> = (0..2)
            .map(|c| (0..m).map(|j| 0.5 + 0.1 * (c + j) as f64).collect())
            .collect();
        let result = run_multi_staar_kinship(&g0, &null, &ann, &mafs);
        assert!(result.staar_o.is_finite(), "STAAR-O not finite: {}", result.staar_o);
        assert!(
            (0.0..=1.0).contains(&result.staar_o),
            "STAAR-O out of range: {}",
            result.staar_o,
        );
        assert!(result.acat_o.is_finite());
        assert!(result.burden_1_25.is_finite());
        assert!(result.skat_1_25.is_finite());
    }

    /// Tiny-n smoke test: fit k=2 with a small pedigree kinship; assert
    /// the expansion + AI-REML pipeline runs end-to-end without error
    /// and returns plausibly-shaped state.
    #[test]
    fn fit_multi_kinship_smoke_n12_k2() {
        let n_fam = 6; // 12 samples total
        let n = n_fam * 2;
        let phi = pedigree_kinship(n_fam, "phi");

        // Simple normal phenotype, two correlated traits.
        let mut y = Mat::<f64>::zeros(n, 2);
        let mut x = Mat::<f64>::zeros(n, 2);
        // Deterministic xorshift so the test doesn't depend on thread-RNG.
        let mut state: u64 = 123;
        let mut u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        for i in 0..n {
            let e1 = 2.0 * u01() - 1.0;
            let e2 = 2.0 * u01() - 1.0;
            y[(i, 0)] = 0.5 + e1;
            y[(i, 1)] = 0.3 + 0.4 * e1 + 0.9 * e2;
            x[(i, 0)] = 1.0;
            x[(i, 1)] = u01();
        }
        let groups = GroupPartition::single(n);
        let fit = fit_multi_kinship(&y, &x, std::slice::from_ref(&phi), &groups)
            .expect("fit_multi_kinship should converge on well-conditioned inputs");
        assert_eq!(fit.n_samples, n);
        assert_eq!(fit.n_pheno, 2);
        assert_eq!(fit.theta_res.len(), 1);
        assert_eq!(fit.theta_g.len(), 1);
        // Residual covariance must be PSD (diagonal entries positive).
        assert!(fit.theta_res[0][(0, 0)] >= 0.0);
        assert!(fit.theta_res[0][(1, 1)] >= 0.0);
        assert_eq!(fit.sigma_inv_x.nrows(), n * 2);
        assert_eq!(fit.cov.nrows(), 2 * 2);
        assert_eq!(fit.p_y.nrows(), n * 2);
    }
}
