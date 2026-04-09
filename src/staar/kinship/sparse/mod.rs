//! Sparse AI-REML backend.
//!
//! When all loaded kinships are sparse and dense Σ⁻¹ would blow the memory
//! budget, this module assembles the sparse Σ via [`assembler::Assembler`],
//! factors it via faer's sparse Cholesky, and hands the factor to the shared
//! [`crate::staar::kinship::reml`] loop.
//!
//! upstream:
//!   - per-iteration AI step → `R/glmmkin.R::R_fitglmm_ai:662-710` (the
//!     sparse-aware variant). Upstream materializes Σ⁻¹ as a dense matrix
//!     because R's `chol2inv` does so even on sparse Cholesky output, then
//!     uses `sum(Sigma_i * kins[[i]])` to exploit kinship sparsity in the
//!     Frobenius product. We do not materialize Σ⁻¹: the default
//!     [`TakahashiBuilder`] runs the Erisman-Tinney recursion to fill in
//!     the entries of Σ⁻¹ at the union sparsity pattern, which is exactly
//!     what the Frobenius sum needs. [`HutchinsonBuilder`] is kept as a
//!     stochastic fallback.
//!   - convergence + boundary refit → shared in `reml.rs::run_reml`.

pub mod assembler;
pub mod hutchinson;
pub mod takahashi;

pub use assembler::Assembler;
pub use hutchinson::{
    sparse_matvec, trace_k_estimate, DEFAULT_HUTCHINSON_PROBES, HUTCHINSON_SEED,
};

use faer::sparse::linalg::solvers::Llt;
use faer::Mat;

use crate::error::CohortError;
use crate::staar::kinship::reml::{
    run_reml, weights_or_ones, SigmaSolver, SolverBuilder, SparseHutchinsonState,
    SparseTakahashiState,
};
use crate::staar::kinship::types::{
    GroupPartition, KinshipInverse, KinshipMatrix, KinshipState, VarianceComponents,
};

/// Which sparse trace estimator to use. `Takahashi` is exact (1:1 with
/// upstream `R/glmmkin.R::R_fitglmm_ai`'s `sum(Sigma_i * kins[[i]])`)
/// and is the default. `Hutchinson` is the stochastic fallback, kept
/// available for the regression test that confirms it still produces
/// valid-within-tolerance estimates.
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq)]
#[allow(dead_code)]
pub enum SparseSolverKind {
    #[default]
    Takahashi,
    Hutchinson,
}

/// Owned wrapper around the faer sparse Cholesky factor of Σ. Cloned into
/// `KinshipState::inverse` (the `Sparse` arm) by [`run_reml`] when the
/// fit takes the sparse path. Once held, every `Σ⁻¹·v` operation in
/// scoring becomes `factor.solve_in_place(v)`.
#[derive(Clone)]
pub struct SparseFactor {
    pub llt: Llt<u32, f64>,
}

/// Builder for the Hutchinson-based sparse path. Each AI iteration
/// rebuilds Σ via the cached assembler and wraps the resulting
/// high-level Llt in `SigmaSolver::SparseHutchinson`. Trace and diagonal
/// of Σ⁻¹ come from stochastic estimators in `hutchinson.rs`.
pub struct HutchinsonBuilder<'a> {
    pub groups: &'a GroupPartition,
    pub weights: &'a [f64],
    pub row_to_group: &'a [usize],
    pub assembler: &'a Assembler,
    pub n_probes: usize,
}

impl<'a> SolverBuilder for HutchinsonBuilder<'a> {
    fn build(&self, tau: &VarianceComponents) -> Result<SigmaSolver, CohortError> {
        let factor = self.assembler.factor(
            tau,
            self.groups,
            self.weights,
            self.row_to_group,
        )?;
        Ok(SigmaSolver::SparseHutchinson(SparseHutchinsonState {
            factor,
            n: self.assembler.n(),
            n_probes: self.n_probes,
        }))
    }
}

/// Builder for the Takahashi-based sparse path. Default sparse builder.
///
/// Each AI iteration:
///   1. assemble fresh Σ values at current τ via the assembler
///   2. refactor numerically via the simplicial Cholesky API (low-level)
///   3. run the Takahashi recursion to fill in the selected inverse
///   4. wrap everything in `SigmaSolver::SparseTakahashi`
///
/// `finalize_inverse` does one extra factorization at the final τ via
/// the high-level Llt API, because the score path
/// (`carrier/sparse_score.rs`) consumes a `KinshipInverse::Sparse`
/// holding a high-level `Llt`, not the simplicial form Takahashi uses
/// internally. The extra cost is one Cholesky factorization per
/// `fit_reml` call (not per iteration), which is negligible.
pub struct TakahashiBuilder<'a> {
    pub groups: &'a GroupPartition,
    pub weights: &'a [f64],
    pub row_to_group: &'a [usize],
    pub assembler: &'a Assembler,
    pub symbolic: &'a takahashi::SymbolicTakahashi,
}

impl<'a> SolverBuilder for TakahashiBuilder<'a> {
    fn build(&self, tau: &VarianceComponents) -> Result<SigmaSolver, CohortError> {
        // Assemble the upper-triangle-only view of Σ at the current τ
        // and hand it to the simplicial Cholesky API. See the comment
        // in `Assembler::new` for why we cannot use the full pattern
        // here.
        let pattern = self.assembler.pattern_owned_upper();
        let values = self.assembler.assemble_values_upper(
            tau,
            self.groups,
            self.weights,
            self.row_to_group,
        );
        let sigma = faer::sparse::SparseColMat::<u32, f64>::new(pattern, values);

        let numeric = takahashi::TakahashiNumeric::factor(self.symbolic, sigma.as_ref())?;
        Ok(SigmaSolver::SparseTakahashi(SparseTakahashiState { numeric }))
    }

    fn finalize_inverse(
        &self,
        _solver: SigmaSolver,
        tau: &VarianceComponents,
    ) -> Result<KinshipInverse, CohortError> {
        // Drop the simplicial state and refactor once via the high-level
        // path so the score test gets a `Llt<u32, f64>` it can call
        // `solve_in_place` on directly.
        let factor = self.assembler.factor(
            tau,
            self.groups,
            self.weights,
            self.row_to_group,
        )?;
        Ok(KinshipInverse::Sparse(SparseFactor { llt: factor }))
    }
}

/// Per-backend trace term for `tr(P K_l)` on the Hutchinson sparse path.
///
/// Computed as `tr(Σ⁻¹ K_l) − tr(Σ⁻¹X cov X' Σ⁻¹ K_l)`, with the first
/// term coming from a stochastic Hutchinson estimator and the second
/// term computed exactly via `B' K_l B` where `B = Σ⁻¹X`. The Takahashi
/// path replaces the first term with an exact computation; this
/// function is kept as the fallback (and as the regression test against
/// Takahashi).
///
/// upstream: `R/glmmkin.R::R_fitglmm_ai:698` —
/// `score[i] <- sum(Y * PAPY) - (sum(Sigma_i*kins[[i]]) - sum(Sigma_iX * crossprod(kins[[i]], Sigma_iXcov)))`.
/// Upstream computes `sum(Sigma_i * kins[[i]])` exactly via the dense
/// inverse it materializes; this function approximates it.
pub fn trace_p_k_sparse_hutchinson(
    state: &SparseHutchinsonState,
    kinship: &KinshipMatrix,
    l: usize,
    sigma_inv_x: &Mat<f64>,
    cov: &Mat<f64>,
) -> f64 {
    let k_sp = kinship
        .as_sparse()
        .expect("sparse trace_p_k called with non-sparse kinship");
    let n = k_sp.nrows();
    let k = sigma_inv_x.ncols();

    // tr(Σ⁻¹ K_l) — Hutchinson stochastic estimate. Mixing constant on
    // top of the global seed decorrelates probes across kinship indices.
    let tr_sinv_k = trace_k_estimate(
        &state.factor,
        k_sp,
        state.n_probes,
        HUTCHINSON_SEED ^ ((l as u64).wrapping_mul(0xC2B2AE3D27D4EB4F)),
    );

    // tr(Σ⁻¹X cov X' Σ⁻¹ K_l). Define B = Σ⁻¹X (n × k). Then
    // tr(B cov B' K_l) = tr(cov B' K_l B). Compute K_l · B column by
    // column via sparse matvec; form B' K_l B (k × k) via dense matmul;
    // then take the trace inner product with cov.
    let mut k_b = Mat::<f64>::zeros(n, k);
    for c in 0..k {
        let mut col_mat = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            col_mat[(i, 0)] = sigma_inv_x[(i, c)];
        }
        let kb_col = sparse_matvec(k_sp, &col_mat);
        for i in 0..n {
            k_b[(i, c)] = kb_col[(i, 0)];
        }
    }
    let bt_k_b = sigma_inv_x.transpose() * &k_b;
    let mut tr_correction = 0.0;
    for a in 0..k {
        for b in 0..k {
            tr_correction += cov[(a, b)] * bt_k_b[(b, a)];
        }
    }

    tr_sinv_k - tr_correction
}

/// Per-backend trace term for `tr(P K_l)` on the Takahashi sparse path.
///
/// Both terms are exact:
///
/// * `tr(Σ⁻¹ K_l)` — read from the selected inverse `Z` at K_l's
///   nonzero positions. Bit-identical to upstream R's
///   `sum(Sigma_i * kins[[i]])` Frobenius product.
/// * `tr(Σ⁻¹X cov X'Σ⁻¹ K_l)` — computed exactly via `B' K_l B` where
///   `B = Σ⁻¹X`, same as the Hutchinson path's correction term.
///
/// upstream: `R/glmmkin.R::R_fitglmm_ai:698`. With this function the
/// sparse path is finally 1:1 with upstream while still using strictly
/// less memory (we never materialize the dense inverse).
pub fn trace_p_k_sparse_takahashi(
    state: &SparseTakahashiState,
    kinship: &KinshipMatrix,
    sigma_inv_x: &Mat<f64>,
    cov: &Mat<f64>,
) -> f64 {
    let k_sp = kinship
        .as_sparse()
        .expect("sparse trace_p_k called with non-sparse kinship");
    let n = k_sp.nrows();
    let k = sigma_inv_x.ncols();

    // Exact tr(Σ⁻¹ K_l) via the selected inverse, summed at K_l's
    // nonzero positions.
    let tr_sinv_k = state.numeric.selected.trace_with(k_sp);

    // tr(Σ⁻¹X cov X' Σ⁻¹ K_l) — same correction as the Hutchinson path,
    // bit-identical floating-point order so the only difference between
    // the two sparse paths is the first term.
    let mut k_b = Mat::<f64>::zeros(n, k);
    for c in 0..k {
        let mut col_mat = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            col_mat[(i, 0)] = sigma_inv_x[(i, c)];
        }
        let kb_col = sparse_matvec(k_sp, &col_mat);
        for i in 0..n {
            k_b[(i, c)] = kb_col[(i, 0)];
        }
    }
    let bt_k_b = sigma_inv_x.transpose() * &k_b;
    let mut tr_correction = 0.0;
    for a in 0..k {
        for b in 0..k {
            tr_correction += cov[(a, b)] * bt_k_b[(b, a)];
        }
    }

    tr_sinv_k - tr_correction
}

/// Sparse fit_reml entry point. Picks between the Takahashi (exact) and
/// Hutchinson (stochastic) backends based on `kind` and dispatches to
/// the shared [`run_reml`] loop. Caller is `fit_reml` in `mod.rs` after
/// it has decided to take the sparse path; the public dispatcher passes
/// `SparseSolverKind::Takahashi` (the default).
///
/// `hutchinson_probes` is the M parameter for the stochastic Hutchinson
/// estimator. Ignored when `kind == Takahashi`. `None` for Hutchinson
/// uses [`DEFAULT_HUTCHINSON_PROBES`].
#[allow(clippy::too_many_arguments)]
pub fn fit_reml_sparse(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: Option<&[f64]>,
    init_tau: VarianceComponents,
    kind: SparseSolverKind,
    hutchinson_probes: Option<usize>,
    budget_bytes: u64,
) -> Result<KinshipState, CohortError> {
    let n = y.nrows();
    let w_vec = weights_or_ones(n, weights);

    // Row → group lookup, used by the assembler to layer the residual
    // term onto the diagonal entries.
    let mut row_to_group = vec![0_usize; n];
    for gi in 0..groups.n_groups() {
        for &row in groups.group(gi) {
            row_to_group[row as usize] = gi;
        }
    }

    let assembler = Assembler::new(n, kinships)?;

    let mut state = match kind {
        SparseSolverKind::Takahashi => {
            let pattern = assembler.pattern_owned_upper();
            let symbolic = takahashi::SymbolicTakahashi::new(pattern.as_ref())?;
            let builder = TakahashiBuilder {
                groups,
                weights: &w_vec,
                row_to_group: &row_to_group,
                assembler: &assembler,
                symbolic: &symbolic,
            };
            run_reml(y, x, kinships, groups, &w_vec, init_tau, &builder)?
        }
        SparseSolverKind::Hutchinson => {
            let builder = HutchinsonBuilder {
                groups,
                weights: &w_vec,
                row_to_group: &row_to_group,
                assembler: &assembler,
                n_probes: hutchinson_probes.unwrap_or(DEFAULT_HUTCHINSON_PROBES),
            };
            run_reml(y, x, kinships, groups, &w_vec, init_tau, &builder)?
        }
    };
    state.budget_bytes = budget_bytes;
    Ok(state)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::staar::kinship::dense::fit_reml_dense;
    use crate::staar::kinship::types::KinshipInverse;
    use faer::sparse::Triplet;

    fn pedigree_kinship(n_fam: usize, sibs: usize) -> KinshipMatrix {
        let n = n_fam * sibs;
        let mut entries: Vec<Triplet<u32, u32, f64>> = Vec::new();
        for f in 0..n_fam {
            let base = f * sibs;
            for i in 0..sibs {
                for j in 0..sibs {
                    let v = if i == j { 1.0 } else { 0.5 };
                    entries.push(Triplet::new((base + i) as u32, (base + j) as u32, v));
                }
            }
        }
        KinshipMatrix::from_triplets(n, entries, "ped".into()).unwrap()
    }

    fn dense_pedigree_kinship(n_fam: usize, sibs: usize) -> KinshipMatrix {
        let n = n_fam * sibs;
        let mut k = Mat::<f64>::zeros(n, n);
        for f in 0..n_fam {
            let base = f * sibs;
            for i in 0..sibs {
                for j in 0..sibs {
                    k[(base + i, base + j)] = if i == j { 1.0 } else { 0.5 };
                }
            }
        }
        KinshipMatrix::Dense {
            matrix: k,
            label: "ped-dense".into(),
        }
    }

    fn synthetic_pheno(n: usize) -> Mat<f64> {
        let mut y = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            let fam = i / 5;
            let g = ((fam as f64 * 0.91).sin()) * 0.6;
            let e = (i as f64 * 1.7).sin() * 0.4;
            y[(i, 0)] = g + e;
        }
        y
    }

    fn intercept(n: usize) -> Mat<f64> {
        let mut x = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }
        x
    }

    #[test]
    fn sparse_assembler_pattern_includes_diagonal() {
        let kin = pedigree_kinship(20, 5);
        let n = kin.n();
        let assembler = Assembler::new(n, std::slice::from_ref(&kin)).unwrap();
        for j in 0..n {
            let s = assembler.col_ptr[j] as usize;
            let e = assembler.col_ptr[j + 1] as usize;
            let mut has_diag = false;
            for k in s..e {
                if assembler.row_idx[k] as usize == j {
                    has_diag = true;
                    break;
                }
            }
            assert!(has_diag, "missing diagonal at column {j}");
        }
    }

    #[test]
    fn sparse_factor_succeeds_at_warm_start() {
        let kin = pedigree_kinship(20, 5);
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let assembler = Assembler::new(n, std::slice::from_ref(&kin)).unwrap();
        let mut tau = VarianceComponents::zeros(1, 1);
        tau.set_kinship(0, 0.5);
        tau.set_group(0, 0.5);
        let weights = vec![1.0; n];
        let row_to_group = vec![0_usize; n];
        let factor = assembler.factor(&tau, &groups, &weights, &row_to_group);
        assert!(factor.is_ok(), "Σ at warm-start τ should be SPD");
    }

    #[test]
    fn hutchinson_trace_reproducible() {
        let kin = pedigree_kinship(40, 5);
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let weights = vec![1.0; n];
        let row_to_group = vec![0_usize; n];
        let assembler = Assembler::new(n, std::slice::from_ref(&kin)).unwrap();
        let mut tau = VarianceComponents::zeros(1, 1);
        tau.set_kinship(0, 0.4);
        tau.set_group(0, 0.6);
        let factor = assembler
            .factor(&tau, &groups, &weights, &row_to_group)
            .unwrap();
        let k_sp = kin.as_sparse().unwrap();
        let t1 = trace_k_estimate(&factor, k_sp, 30, 12345);
        let t2 = trace_k_estimate(&factor, k_sp, 30, 12345);
        assert_eq!(t1, t2, "same seed must produce identical estimate");
    }

    #[test]
    fn sparse_fit_reml_pedigree_matches_dense_within_tolerance() {
        // 20 families × 5 sibs = 100 samples. Sparse path with M=200
        // probes should land within ~5% of the dense REML τ values.
        let kin_sparse = pedigree_kinship(20, 5);
        let kin_dense = dense_pedigree_kinship(20, 5);
        let n = kin_sparse.n();
        let groups = GroupPartition::single(n);
        let y = synthetic_pheno(n);
        let x = intercept(n);

        let weights = vec![1.0_f64; n];

        let y_mean = (0..n).map(|i| y[(i, 0)]).sum::<f64>() / n as f64;
        let y_var = (0..n)
            .map(|i| (y[(i, 0)] - y_mean).powi(2))
            .sum::<f64>()
            / (n as f64 - 1.0).max(1.0);
        let mut warm_init = VarianceComponents::zeros(1, 1);
        let warm = y_var / 2.0;
        warm_init.set_kinship(0, warm / kin_sparse.mean_diagonal());
        warm_init.set_group(0, warm);

        let budget = crate::staar::kinship::DEFAULT_KINSHIP_MEM_BYTES;
        let dense_state = fit_reml_dense(
            &y,
            &x,
            std::slice::from_ref(&kin_dense),
            &groups,
            &weights,
            warm_init.clone(),
            budget,
        )
        .expect("dense fit");

        let sparse_state = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin_sparse),
            &groups,
            None,
            warm_init,
            SparseSolverKind::Hutchinson,
            Some(200),
            budget,
        )
        .expect("sparse fit");

        let dense_tau_k = dense_state.tau.kinship(0);
        let sparse_tau_k = sparse_state.tau.kinship(0);
        let dense_tau_g = dense_state.tau.group(0);
        let sparse_tau_g = sparse_state.tau.group(0);

        let scale = (dense_tau_k.abs() + dense_tau_g.abs()).max(1e-3);
        let err_k = (dense_tau_k - sparse_tau_k).abs() / scale;
        let err_g = (dense_tau_g - sparse_tau_g).abs() / scale;

        assert!(
            err_k < 0.15,
            "τ_kinship: dense={dense_tau_k:.4} sparse={sparse_tau_k:.4} rel_err={err_k:.3}"
        );
        assert!(
            err_g < 0.15,
            "τ_group:   dense={dense_tau_g:.4} sparse={sparse_tau_g:.4} rel_err={err_g:.3}"
        );

        assert!(matches!(sparse_state.inverse, KinshipInverse::Sparse(_)));
    }

    #[test]
    fn sparse_fit_reml_takahashi_matches_dense() {
        // Same fixture as the Hutchinson tolerance test, but with the
        // exact Takahashi backend. Should match the dense path to ~1e-8
        // (limited only by floating-point reordering between the dense
        // and simplicial Cholesky paths).
        let kin_sparse = pedigree_kinship(20, 5);
        let kin_dense = dense_pedigree_kinship(20, 5);
        let n = kin_sparse.n();
        let groups = GroupPartition::single(n);
        let y = synthetic_pheno(n);
        let x = intercept(n);

        let weights = vec![1.0_f64; n];

        let y_mean = (0..n).map(|i| y[(i, 0)]).sum::<f64>() / n as f64;
        let y_var = (0..n)
            .map(|i| (y[(i, 0)] - y_mean).powi(2))
            .sum::<f64>()
            / (n as f64 - 1.0).max(1.0);
        let mut warm_init = VarianceComponents::zeros(1, 1);
        let warm = y_var / 2.0;
        warm_init.set_kinship(0, warm / kin_sparse.mean_diagonal());
        warm_init.set_group(0, warm);

        let budget = crate::staar::kinship::DEFAULT_KINSHIP_MEM_BYTES;
        let dense_state = fit_reml_dense(
            &y,
            &x,
            std::slice::from_ref(&kin_dense),
            &groups,
            &weights,
            warm_init.clone(),
            budget,
        )
        .expect("dense fit");

        let sparse_state = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin_sparse),
            &groups,
            None,
            warm_init,
            SparseSolverKind::Takahashi,
            None,
            budget,
        )
        .expect("takahashi fit");

        let dense_tau_k = dense_state.tau.kinship(0);
        let sparse_tau_k = sparse_state.tau.kinship(0);
        let dense_tau_g = dense_state.tau.group(0);
        let sparse_tau_g = sparse_state.tau.group(0);

        let scale = (dense_tau_k.abs() + dense_tau_g.abs()).max(1e-3);
        let err_k = (dense_tau_k - sparse_tau_k).abs() / scale;
        let err_g = (dense_tau_g - sparse_tau_g).abs() / scale;

        // Tight tolerance — Takahashi is exact, so the only difference
        // is floating-point reordering between dense and simplicial
        // Cholesky paths. 1e-8 leaves ample headroom.
        assert!(
            err_k < 1e-8,
            "τ_kinship: dense={dense_tau_k:.10} takahashi={sparse_tau_k:.10} rel_err={err_k:.3e}"
        );
        assert!(
            err_g < 1e-8,
            "τ_group:   dense={dense_tau_g:.10} takahashi={sparse_tau_g:.10} rel_err={err_g:.3e}"
        );

        // Result still carries a high-level Llt for the score path.
        assert!(matches!(sparse_state.inverse, KinshipInverse::Sparse(_)));
    }

    #[test]
    fn sparse_fit_reml_reproducible() {
        let kin = pedigree_kinship(15, 4);
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let y = synthetic_pheno(n);
        let x = intercept(n);

        let mut warm_init = VarianceComponents::zeros(1, 1);
        warm_init.set_kinship(0, 0.5);
        warm_init.set_group(0, 0.5);

        let budget = crate::staar::kinship::DEFAULT_KINSHIP_MEM_BYTES;
        let s1 = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin),
            &groups,
            None,
            warm_init.clone(),
            SparseSolverKind::Hutchinson,
            Some(64),
            budget,
        )
        .expect("first sparse fit");
        let s2 = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin),
            &groups,
            None,
            warm_init,
            SparseSolverKind::Hutchinson,
            Some(64),
            budget,
        )
        .expect("second sparse fit");

        assert_eq!(s1.tau.kinship(0), s2.tau.kinship(0));
        assert_eq!(s1.tau.group(0), s2.tau.group(0));
        assert_eq!(s1.n_iter, s2.n_iter);
    }
}
