//! Shared AI-REML algorithm and convergence loop for both backends.
//!
//! Ports `glmmkin.ai` from upstream GMMAT
//! (`hanchenphd/GMMAT@473b342 R/glmmkin.R:281-422`). The dense and sparse
//! backends differ in *how* they apply Σ⁻¹ to a vector and *how* they
//! compute the trace term `tr(P K_l)`. Everything else — the score vector
//! formula, the average information matrix, the step-halving and boundary
//! refit loops, the convergence criterion — is identical and lives here.
//!
//! Two abstractions tie the backends in:
//!
//! * [`SolverBuilder`] — per-iteration "freeze Σ at current τ" hook.
//!   Each iteration the convergence loop calls `builder.build(&tau)`,
//!   which returns a [`SigmaSolver`] holding either a materialized dense
//!   inverse or a sparse Cholesky factor.
//!
//! * [`SigmaSolver`] — the in-iteration handle the AI step uses to apply
//!   Σ⁻¹ to vectors, compute `diag(Σ⁻¹)`, and compute `tr(P K_l)`. Each
//!   variant routes to a backend-specific free function for `trace_p_k`
//!   so dense and sparse keep their original floating-point order
//!   (statistical invariance).
//!
//! At the bottom: [`run_reml`] is the entry point each backend calls. It
//! does step-halving and boundary refit on top of the AI iteration.
//!
//! ## Variable rename map (upstream R → here)
//!
//! See `UPSTREAM.md` for the full table. The most common ones:
//!
//! ```text
//! upstream R   →  here
//! Y            →  y
//! tau / theta  →  tau (kinship slots first, group slots last)
//! Sigma_i      →  inverse of Σ (in `SigmaSolver`)
//! Sigma_iX     →  sigma_inv_x
//! cov          →  cov
//! PY           →  p_y
//! diagP        →  diag_p
//! score        →  score
//! AI           →  ai
//! Dtau         →  d_tau
//! fixtau       →  fixtau
//! ```

use faer::prelude::Solve;
use faer::Mat;

use crate::error::FavorError;
use crate::staar::kinship::dense;
use crate::staar::kinship::sparse;
use crate::staar::kinship::sparse::takahashi::TakahashiNumeric;
use crate::staar::kinship::sparse::SparseFactor;
use crate::staar::kinship::types::{
    GroupPartition, KinshipInverse, KinshipMatrix, KinshipState, VarianceComponents,
};

/// REML convergence tolerance. Upstream `glmmkin.R:122` `tol = 1e-5`.
pub const REML_TOL: f64 = 1e-5;
/// Inner AI iteration cap. Upstream `glmmkin.R:122` `maxiter = 500`.
pub const REML_MAX_ITER: usize = 500;
/// Outer boundary-refit cap. Upstream's refit `while` loop in
/// `glmmkin.R:194-210` is unbounded; this is our safety net against
/// pathological oscillation. On hit we return an Analysis error.
pub const REML_MAX_OUTER_REFITS: usize = 10;
/// Step-halving attempts before giving up. Upstream's halving loop in
/// `glmmkin.R:362-366` is unbounded; this is our safety net.
pub const REML_STEP_HALVE_MAX: usize = 50;
/// Boundary factor for declaring τ at zero. Upstream uses the literal
/// `1.01 * tol` throughout `glmmkin.R` (e.g. lines 187, 191, 203, 207)
/// to flag a variance component as collapsed onto the boundary.
pub const BOUNDARY_FACTOR: f64 = 1.01;

/// Force `weights` into an owned vector, defaulting to all-ones.
pub fn weights_or_ones(n: usize, weights: Option<&[f64]>) -> Vec<f64> {
    match weights {
        Some(w) => {
            assert_eq!(w.len(), n);
            w.to_vec()
        }
        None => vec![1.0; n],
    }
}

/// Solve an SPD system via faer's column-pivoted QR. Used for the small
/// k × k systems (cov, AI free block) where pivoting is cheap and the
/// numerical robustness is worth the overhead vs Cholesky.
pub fn solve_spd(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    a.col_piv_qr().solve(b)
}

/// Invert an SPD matrix by solving against the identity.
pub fn invert_spd(a: &Mat<f64>) -> Mat<f64> {
    let n = a.nrows();
    debug_assert_eq!(n, a.ncols());
    solve_spd(a, &Mat::<f64>::identity(n, n))
}

/// Frobenius inner product `<A, B> = Σ_{ij} A[i,j] · B[i,j]`. Used by the
/// dense path to compute `tr(A^T B)` in the trace correction term without
/// forming the product matrix.
pub fn frob_inner(a: &Mat<f64>, b: &Mat<f64>) -> f64 {
    let (rows, cols) = (a.nrows(), a.ncols());
    debug_assert_eq!(rows, b.nrows());
    debug_assert_eq!(cols, b.ncols());
    let mut s = 0.0;
    for r in 0..rows {
        for c in 0..cols {
            s += a[(r, c)] * b[(r, c)];
        }
    }
    s
}

/// Manual matrix-vector product, kept as the canonical "dense Σ⁻¹ · v"
/// arithmetic. Bit-identical to the original `dense.rs::matvec`.
pub fn matvec(m: &Mat<f64>, v: &Mat<f64>) -> Mat<f64> {
    debug_assert_eq!(m.ncols(), v.nrows());
    debug_assert_eq!(v.ncols(), 1);
    let n = m.nrows();
    let mut out = Mat::<f64>::zeros(n, 1);
    for r in 0..n {
        let mut s = 0.0;
        for c in 0..m.ncols() {
            s += m[(r, c)] * v[(c, 0)];
        }
        out[(r, 0)] = s;
    }
    out
}

/// Polymorphic `K_l · v` for either dense or sparse kinship storage.
/// Each variant uses the bit-identical arithmetic from the original
/// `dense.rs::matvec` and `sparse.rs::sparse_matvec`.
pub fn matvec_kinship(k: &KinshipMatrix, v: &Mat<f64>) -> Mat<f64> {
    debug_assert_eq!(v.ncols(), 1);
    let n = v.nrows();
    let mut out = Mat::<f64>::zeros(n, 1);
    match k {
        KinshipMatrix::Dense { matrix, .. } => {
            debug_assert_eq!(matrix.nrows(), n);
            for r in 0..n {
                let mut s = 0.0;
                for c in 0..n {
                    s += matrix[(r, c)] * v[(c, 0)];
                }
                out[(r, 0)] = s;
            }
        }
        KinshipMatrix::Sparse { matrix, .. } => {
            let cp = matrix.symbolic().col_ptr();
            let ri = matrix.symbolic().row_idx();
            let val = matrix.val();
            for j in 0..n {
                let v_j = v[(j, 0)];
                let s = cp[j] as usize;
                let e = cp[j + 1] as usize;
                for kk in s..e {
                    let row = ri[kk] as usize;
                    out[(row, 0)] += val[kk] * v_j;
                }
            }
        }
    }
    out
}

/// In-iteration handle for applying Σ⁻¹ and reading the entries the
/// score / AI matrix need. Built fresh each AI iteration by the
/// [`SolverBuilder`].
///
/// Three variants:
///
/// * `Dense` — owns a materialized n × n inverse. tr and diag are exact.
///   Cost: O(n²) memory.
/// * `SparseHutchinson` — owns a high-level sparse Cholesky factor of Σ.
///   Σ⁻¹ · v is a triangular solve; tr(Σ⁻¹ K_l) and diag(Σ⁻¹) are
///   stochastic Hutchinson estimates. Kept as a fallback (and for the
///   parity test that confirms the stochastic path still produces
///   valid-within-tolerance estimates).
/// * `SparseTakahashi` — owns the simplicial Cholesky factor plus the
///   selected inverse Z = Σ⁻¹ at the union sparsity pattern, computed by
///   the Takahashi recursion. tr and diag are *exact* — 1:1 with
///   upstream `R/glmmkin.R::R_fitglmm_ai`'s `sum(Sigma_i * kins[[i]])`
///   formula, while using strictly less memory than upstream (we never
///   materialize the dense inverse). Default sparse path.
pub enum SigmaSolver {
    Dense(Mat<f64>),
    SparseHutchinson(SparseHutchinsonState),
    SparseTakahashi(SparseTakahashiState),
}

/// Per-iteration state for the Hutchinson-based sparse path. The factor
/// applies Σ⁻¹ via triangular solves; the probe count and `n` feed the
/// stochastic trace estimators.
pub struct SparseHutchinsonState {
    pub factor: faer::sparse::linalg::solvers::Llt<u32, f64>,
    pub n: usize,
    pub n_probes: usize,
}

/// Per-iteration state for the Takahashi-based sparse path. Owns the
/// simplicial Cholesky factor (for solves via `SimplicialLltRef`) and
/// the selected inverse `Z` at the L sparsity pattern (for exact trace
/// and diagonal reads).
pub struct SparseTakahashiState {
    pub numeric: TakahashiNumeric,
}

impl SigmaSolver {
    /// Σ⁻¹ · v_col where `v_col` is `n × 1`. Returns a fresh column vector.
    /// Dense path uses the manual matvec loop (bit-identical to the
    /// original `dense.rs::matvec`); sparse paths use a sparse triangular
    /// solve via the cached Cholesky factor.
    pub fn solve_into(&self, v: &Mat<f64>) -> Mat<f64> {
        debug_assert_eq!(v.ncols(), 1);
        match self {
            Self::Dense(sigma_inv) => matvec(sigma_inv, v),
            Self::SparseHutchinson(state) => {
                use faer::linalg::solvers::Solve;
                let mut out = v.clone();
                state.factor.solve_in_place(&mut out);
                out
            }
            Self::SparseTakahashi(state) => {
                let mut out = v.clone();
                state.numeric.solve_in_place(out.as_mut());
                out
            }
        }
    }

    /// Σ⁻¹ · M for an `n × k` matrix M. Returns a fresh matrix. Dense uses
    /// faer's matrix-matrix product (matches `dense.rs`'s `&sigma_inv * x`);
    /// sparse paths clone M and run the Cholesky solve in place.
    pub fn solve_columns(&self, m: &Mat<f64>) -> Mat<f64> {
        match self {
            Self::Dense(sigma_inv) => sigma_inv * m,
            Self::SparseHutchinson(state) => {
                use faer::linalg::solvers::Solve;
                let mut out = m.clone();
                state.factor.solve_in_place(&mut out);
                out
            }
            Self::SparseTakahashi(state) => {
                let mut out = m.clone();
                state.numeric.solve_in_place(out.as_mut());
                out
            }
        }
    }

    /// Full `diag(Σ⁻¹)` as a length-n vector.
    ///
    /// * Dense — exact, reads the diagonal of the materialized inverse.
    /// * SparseHutchinson — stochastic Rademacher estimator.
    /// * SparseTakahashi — exact, reads the diagonal of the selected
    ///   inverse computed by the Takahashi recursion. Bit-identical
    ///   replacement for upstream's `diag(Sigma_i)`.
    pub fn diag_sigma_inv(&self) -> Vec<f64> {
        match self {
            Self::Dense(sigma_inv) => {
                let n = sigma_inv.nrows();
                (0..n).map(|i| sigma_inv[(i, i)]).collect()
            }
            Self::SparseHutchinson(state) => sparse::hutchinson::diag_inverse_estimate(
                &state.factor,
                state.n,
                state.n_probes,
                sparse::hutchinson::HUTCHINSON_SEED,
            ),
            Self::SparseTakahashi(state) => state.numeric.selected.diag(),
        }
    }

    /// `tr(P K_l) = tr(Σ⁻¹ K_l) − tr(Σ⁻¹ X cov X' Σ⁻¹ K_l)`.
    ///
    /// Each backend computes this its own way to keep the floating-point
    /// summation order bit-identical to the corresponding original code
    /// path:
    ///
    /// * Dense — `tr(Σ⁻¹ K_l) − frob(Σ⁻¹X, K_l Σ⁻¹X cov)`, both terms
    ///   exact via the materialized inverse.
    /// * SparseHutchinson — `hutchinson_trace − Σ_{a,b} cov[a,b]·(B'K_l B)[b,a]`,
    ///   first term stochastic, second term exact.
    /// * SparseTakahashi — `Z.trace_with(K_l) − Σ_{a,b} cov[a,b]·(B'K_l B)[b,a]`,
    ///   both terms exact. The first term reads the selected inverse Z
    ///   directly at K_l's pattern, which matches upstream R's
    ///   `sum(Sigma_i * kins[[i]])` Frobenius product.
    pub fn trace_p_k(
        &self,
        kinship: &KinshipMatrix,
        l: usize,
        sigma_inv_x: &Mat<f64>,
        cov: &Mat<f64>,
    ) -> f64 {
        match self {
            Self::Dense(sigma_inv) => {
                dense::trace_p_k_dense(sigma_inv, kinship, sigma_inv_x, cov)
            }
            Self::SparseHutchinson(state) => {
                sparse::trace_p_k_sparse_hutchinson(state, kinship, l, sigma_inv_x, cov)
            }
            Self::SparseTakahashi(state) => {
                sparse::trace_p_k_sparse_takahashi(state, kinship, sigma_inv_x, cov)
            }
        }
    }

    /// Hand back the inverse representation for embedding in `KinshipState`.
    /// For the Takahashi variant the caller is expected to use
    /// [`SolverBuilder::finalize_inverse`] instead, which can refactor via
    /// the high-level Llt path needed by the score test.
    pub fn into_inverse(self) -> KinshipInverse {
        match self {
            Self::Dense(sigma_inv) => KinshipInverse::Dense(sigma_inv),
            Self::SparseHutchinson(state) => {
                KinshipInverse::Sparse(SparseFactor { llt: state.factor })
            }
            Self::SparseTakahashi(_) => panic!(
                "SigmaSolver::SparseTakahashi must go through SolverBuilder::finalize_inverse, \
                 not into_inverse — the score path needs a high-level Llt that the simplicial \
                 representation doesn't carry"
            ),
        }
    }
}

/// Builds a fresh [`SigmaSolver`] from the current variance components.
/// Implemented once per backend (`DenseBuilder` in `dense.rs`,
/// `HutchinsonBuilder` and `TakahashiBuilder` in `sparse/mod.rs`). The
/// shared convergence loop calls `build(&tau)` once per AI iteration.
///
/// `finalize_inverse` is the hook that converts the final per-iteration
/// solver into the [`KinshipInverse`] representation that gets carried in
/// `KinshipState` and consumed by the score test. The default
/// implementation just calls `solver.into_inverse()`. The Takahashi
/// builder overrides it because the simplicial Cholesky factor inside the
/// Takahashi solver isn't a high-level `Llt`, and the score path needs
/// the high-level form — so the override refactors Σ once at the final τ
/// to produce a fresh `Llt` for the score path.
pub trait SolverBuilder {
    fn build(&self, tau: &VarianceComponents) -> Result<SigmaSolver, FavorError>;

    fn finalize_inverse(
        &self,
        solver: SigmaSolver,
        tau: &VarianceComponents,
    ) -> Result<KinshipInverse, FavorError> {
        let _ = tau;
        Ok(solver.into_inverse())
    }
}

/// Per-iteration scratch returned by [`ai_step`].
pub struct AiStep {
    pub d_tau: Vec<f64>,
    pub alpha: Mat<f64>,
    pub p_y: Mat<f64>,
    pub sigma_inv_x: Mat<f64>,
    pub cov: Mat<f64>,
    pub solver: SigmaSolver,
}

/// One AI-REML iteration step. Backend-agnostic: every Σ⁻¹ application
/// goes through the [`SigmaSolver`] passed in, and the trace term
/// `tr(P K_l)` is dispatched per backend so each path keeps its original
/// floating-point arithmetic.
///
/// upstream: `R/glmmkin.R::R_fitglmm_ai_dense:712-769` for the dense form,
/// `R/glmmkin.R::R_fitglmm_ai:662-710` for the sparse form. Both compute
/// the same score and AI matrix:
///
/// ```text
/// score[l]   = (PY)' K_l (PY) − tr(P K_l)
/// score[L+g] = Σ_{i∈g} (PY[i])²/W[i] − Σ_{i∈g} diag_p[i]/W[i]
/// AI[i,j]    = (PY)' V_i P V_j (PY)
/// ```
///
/// where `V_l = K_l` for kinship components and `V_g = diag(1/W) over g`
/// for group components.
#[allow(clippy::too_many_arguments)]
pub fn ai_step(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    tau: &VarianceComponents,
    fixtau: &[bool],
    solver: SigmaSolver,
) -> Result<AiStep, FavorError> {
    let n = y.nrows();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    let k = x.ncols();
    debug_assert_eq!(fixtau.len(), n_comp);
    debug_assert_eq!(tau.n_total(), n_comp);
    let _ = weights;

    // Σ⁻¹ X via column-wise solve.
    let sigma_inv_x = solver.solve_columns(x);

    // (X' Σ⁻¹ X)⁻¹.
    let xt_sinv_x = x.transpose() * &sigma_inv_x;
    let cov = invert_spd(&xt_sinv_x);

    // α̂ = cov · X' · Σ⁻¹ · y, η = X α̂, residual = y − η, PY = Σ⁻¹ · residual.
    let sinv_y = solver.solve_into(y);
    let xt_sinv_y = x.transpose() * &sinv_y;
    let alpha = &cov * &xt_sinv_y;
    let eta = x * &alpha;
    let mut residual = Mat::<f64>::zeros(n, 1);
    for i in 0..n {
        residual[(i, 0)] = y[(i, 0)] - eta[(i, 0)];
    }
    let p_y = solver.solve_into(&residual);

    // diag(P) = diag(Σ⁻¹) − Σ_{a,b} sigma_inv_x[i,a] · cov[a,b] · sigma_inv_x[i,b].
    // The first term is exact for dense, stochastic (Hutchinson) for sparse;
    // the second term is always exact and computed inline.
    let diag_sigma_inv = solver.diag_sigma_inv();
    let mut diag_p = vec![0.0_f64; n];
    for i in 0..n {
        let mut hi = 0.0;
        for a in 0..k {
            for b in 0..k {
                hi += sigma_inv_x[(i, a)] * cov[(a, b)] * sigma_inv_x[(i, b)];
            }
        }
        diag_p[i] = diag_sigma_inv[i] - hi;
    }

    // K_l · PY for each kinship component, plus the score and the trace
    // correction.
    let mut score = vec![0.0_f64; n_comp];
    let mut k_py: Vec<Mat<f64>> = Vec::with_capacity(l);
    for li in 0..l {
        let v = matvec_kinship(&kinships[li], &p_y);
        let tr_p_k = solver.trace_p_k(&kinships[li], li, &sigma_inv_x, &cov);
        let mut py_k_py = 0.0;
        for i in 0..n {
            py_k_py += p_y[(i, 0)] * v[(i, 0)];
        }
        score[li] = py_k_py - tr_p_k;
        k_py.push(v);
    }

    // Group score: U_g = Σ_{i∈g} (PY[i])²/W[i] − Σ_{i∈g} diag_p[i]/W[i].
    for gi in 0..g {
        let mut s_data = 0.0;
        let mut s_trace = 0.0;
        for &row in groups.group(gi) {
            let i = row as usize;
            let w_i = weights[i];
            s_data += p_y[(i, 0)] * p_y[(i, 0)] / w_i;
            s_trace += diag_p[i] / w_i;
        }
        score[l + gi] = s_data - s_trace;
    }

    // AI matrix: AI[i,j] = (PY)' V_i P V_j (PY). Project K_l · PY and
    // V_g · PY through P (which routes solve_into to the backend), then
    // build the symmetric AI block.
    let mut ai = Mat::<f64>::zeros(n_comp, n_comp);

    let p_k_py: Vec<Mat<f64>> = (0..l)
        .map(|li| apply_projection(&solver, &sigma_inv_x, &cov, x, &k_py[li]))
        .collect();

    let mut v_py: Vec<Mat<f64>> = Vec::with_capacity(g);
    let mut p_v_py: Vec<Mat<f64>> = Vec::with_capacity(g);
    for gi in 0..g {
        let mut vpy = Mat::<f64>::zeros(n, 1);
        for &row in groups.group(gi) {
            let i = row as usize;
            vpy[(i, 0)] = p_y[(i, 0)] / weights[i];
        }
        let p_vpy = apply_projection(&solver, &sigma_inv_x, &cov, x, &vpy);
        v_py.push(vpy);
        p_v_py.push(p_vpy);
    }

    // Kinship × kinship block.
    for li in 0..l {
        for lj in 0..l {
            let mut s = 0.0;
            for i in 0..n {
                s += k_py[li][(i, 0)] * p_k_py[lj][(i, 0)];
            }
            ai[(li, lj)] = s;
        }
    }
    // Kinship × group block (and its transpose by symmetry).
    for li in 0..l {
        for gj in 0..g {
            let mut s = 0.0;
            for i in 0..n {
                s += k_py[li][(i, 0)] * p_v_py[gj][(i, 0)];
            }
            ai[(li, l + gj)] = s;
            ai[(l + gj, li)] = s;
        }
    }
    // Group × group block. Sums only over the (gi) group's row indices —
    // V_g is zero outside its group.
    for gi in 0..g {
        for gj in 0..g {
            let mut s = 0.0;
            for &row in groups.group(gi) {
                let i = row as usize;
                s += v_py[gi][(i, 0)] * p_v_py[gj][(i, 0)];
            }
            ai[(l + gi, l + gj)] = s;
        }
    }

    // Solve AI · Δτ = score over the free components only.
    let d_tau = solve_dtau(&ai, &score, fixtau);

    Ok(AiStep {
        d_tau,
        alpha,
        p_y,
        sigma_inv_x,
        cov,
        solver,
    })
}

/// Project a vector through `P = Σ⁻¹ − Σ⁻¹ X (X' Σ⁻¹ X)⁻¹ X' Σ⁻¹`. Calls
/// [`SigmaSolver::solve_into`] for the Σ⁻¹ application; everything else
/// is a small dense matmul on `(n × k)` and `(k × k)` matrices.
pub fn apply_projection(
    solver: &SigmaSolver,
    sigma_inv_x: &Mat<f64>,
    cov: &Mat<f64>,
    x: &Mat<f64>,
    v: &Mat<f64>,
) -> Mat<f64> {
    let sinv_v = solver.solve_into(v);
    let xt_sinv_v = x.transpose() * &sinv_v;
    let cov_xt_sinv_v = cov * &xt_sinv_v;
    let proj = sigma_inv_x * &cov_xt_sinv_v;
    let n = v.nrows();
    let mut out = Mat::<f64>::zeros(n, 1);
    for i in 0..n {
        out[(i, 0)] = sinv_v[(i, 0)] - proj[(i, 0)];
    }
    out
}

/// Solve `AI · Δτ = score` over the free components only (entries flagged
/// as fixed by the boundary refit get a zero step). Pulls the free block
/// out of `ai`, solves the small system, and scatters the answer back.
fn solve_dtau(ai: &Mat<f64>, score: &[f64], fixtau: &[bool]) -> Vec<f64> {
    let n_comp = score.len();
    debug_assert_eq!(ai.nrows(), n_comp);
    debug_assert_eq!(ai.ncols(), n_comp);
    debug_assert_eq!(fixtau.len(), n_comp);

    let n_free = fixtau.iter().filter(|&&f| !f).count();
    let mut d_tau = vec![0.0_f64; n_comp];
    if n_free == 0 {
        return d_tau;
    }
    let idx_map: Vec<usize> = (0..n_comp).filter(|&i| !fixtau[i]).collect();
    let mut ai_free = Mat::<f64>::zeros(n_free, n_free);
    let mut score_free = Mat::<f64>::zeros(n_free, 1);
    for (a, &fi) in idx_map.iter().enumerate() {
        score_free[(a, 0)] = score[fi];
        for (b, &fj) in idx_map.iter().enumerate() {
            ai_free[(a, b)] = ai[(fi, fj)];
        }
    }
    let dt_free = solve_spd(&ai_free, &score_free);
    for (a, &fi) in idx_map.iter().enumerate() {
        d_tau[fi] = dt_free[(a, 0)];
    }
    d_tau
}

/// Inner AI-REML iteration loop. Runs `ai_step` until the convergence
/// criterion is met or `REML_MAX_ITER` is exhausted. Step-halving is
/// applied per upstream `glmmkin.R:362-366` whenever the Newton step
/// would push a free τ negative.
#[allow(clippy::too_many_arguments)]
fn converge<B: SolverBuilder>(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    init_tau: VarianceComponents,
    fixtau: &[bool],
    builder: &B,
) -> Result<KinshipState, FavorError> {
    let n = y.nrows();
    let k = x.ncols();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    debug_assert_eq!(init_tau.n_total(), n_comp);
    debug_assert_eq!(fixtau.len(), n_comp);
    let _ = n;

    let mut tau = init_tau;
    let mut alpha_prev = Mat::<f64>::zeros(k, 1);
    let mut iter_used = 0;
    let mut last: Option<AiStep> = None;

    for iter in 0..REML_MAX_ITER {
        iter_used = iter + 1;
        let tau0 = tau.clone();
        let solver = builder.build(&tau)?;
        let step_out =
            ai_step(y, x, kinships, groups, weights, &tau, fixtau, solver)?;

        // Newton step + step-halving on negative τ. Matches upstream
        // `glmmkin.R:357-367`. Free components update; fixed stay put.
        let mut step = step_out.d_tau.clone();
        let mut try_count = 0;
        loop {
            for i in 0..n_comp {
                if !fixtau[i] {
                    tau.as_slice_mut()[i] = tau0.as_slice()[i] + step[i];
                }
            }
            for i in 0..n_comp {
                let v = tau.as_slice()[i];
                if v < REML_TOL && tau0.as_slice()[i] < REML_TOL {
                    tau.as_slice_mut()[i] = 0.0;
                }
            }
            if tau.as_slice().iter().all(|&t| t >= 0.0) {
                break;
            }
            for s in step.iter_mut() {
                *s /= 2.0;
            }
            try_count += 1;
            if try_count > REML_STEP_HALVE_MAX {
                return Err(FavorError::Analysis(format!(
                    "AI-REML step-halving failed to find a feasible τ in {REML_STEP_HALVE_MAX} attempts"
                )));
            }
        }
        for v in tau.as_slice_mut().iter_mut() {
            if *v < REML_TOL {
                *v = 0.0;
            }
        }

        // Convergence: 2 · max(rel_α, rel_τ) < REML_TOL. Upstream
        // `glmmkin.R:405`.
        let mut rel_alpha = 0.0_f64;
        for i in 0..k {
            let denom = step_out.alpha[(i, 0)].abs() + alpha_prev[(i, 0)].abs() + REML_TOL;
            let r = (step_out.alpha[(i, 0)] - alpha_prev[(i, 0)]).abs() / denom;
            if r > rel_alpha {
                rel_alpha = r;
            }
        }
        let mut rel_tau = 0.0_f64;
        for i in 0..n_comp {
            let denom = tau.as_slice()[i].abs() + tau0.as_slice()[i].abs() + REML_TOL;
            let r = (tau.as_slice()[i] - tau0.as_slice()[i]).abs() / denom;
            if r > rel_tau {
                rel_tau = r;
            }
        }

        alpha_prev = step_out.alpha.clone();
        last = Some(step_out);

        if 2.0 * rel_alpha.max(rel_tau) < REML_TOL {
            break;
        }

        // Divergence guard. Upstream `glmmkin.R:406-410` warns and
        // breaks; we error out so the caller knows the fit is no good.
        if tau.as_slice().iter().any(|&t| t.abs() > REML_TOL.powi(-2)) {
            return Err(FavorError::Analysis(format!(
                "AI-REML diverged: |τ_max| > {}",
                REML_TOL.powi(-2)
            )));
        }
    }

    let last = last
        .ok_or_else(|| FavorError::Analysis("AI-REML produced no iterations".into()))?;

    let total_var: f64 = tau.as_slice().iter().sum();
    let h2: Vec<f64> = if total_var > 0.0 {
        (0..l).map(|i| tau.kinship(i) / total_var).collect()
    } else {
        vec![0.0; l]
    };

    let inverse = builder.finalize_inverse(last.solver, &tau)?;

    Ok(KinshipState {
        tau,
        inverse,
        sigma_inv_x: last.sigma_inv_x,
        cov: last.cov,
        p_y: last.p_y,
        h2,
        n_iter: iter_used,
        outer_refits: 0,
    })
}

/// Run AI-REML to convergence with the boundary-refit outer loop.
/// Equivalent to upstream `glmmkin.fit + glmmkin.ai` taken together
/// (`glmmkin.R:122-279` + `glmmkin.R:281-422`). On each refit, free
/// components that have collapsed below `BOUNDARY_FACTOR · REML_TOL` get
/// pinned at zero and the inner converge loop reruns. Hits an error if
/// the refit oscillates beyond `REML_MAX_OUTER_REFITS`.
pub fn run_reml<B: SolverBuilder>(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    init_tau: VarianceComponents,
    builder: &B,
) -> Result<KinshipState, FavorError> {
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;

    let mut fixtau = vec![false; n_comp];
    let mut state = converge(
        y,
        x,
        kinships,
        groups,
        weights,
        init_tau,
        &fixtau,
        builder,
    )?;

    for refit_iter in 0..REML_MAX_OUTER_REFITS {
        let mut changed = false;
        for (i, fixed) in fixtau.iter_mut().enumerate().take(n_comp) {
            if !*fixed && state.tau.as_slice()[i] < BOUNDARY_FACTOR * REML_TOL {
                *fixed = true;
                changed = true;
            }
        }
        if !changed {
            state.outer_refits = refit_iter;
            return Ok(state);
        }
        let mut tau_refit = state.tau.clone();
        for (i, &fixed) in fixtau.iter().enumerate().take(n_comp) {
            if fixed {
                tau_refit.as_slice_mut()[i] = 0.0;
            }
        }
        state = converge(
            y,
            x,
            kinships,
            groups,
            weights,
            tau_refit,
            &fixtau,
            builder,
        )?;
        state.outer_refits = refit_iter + 1;
    }

    Err(FavorError::Analysis(format!(
        "AI-REML boundary refit did not stabilize within {REML_MAX_OUTER_REFITS} iterations"
    )))
}
