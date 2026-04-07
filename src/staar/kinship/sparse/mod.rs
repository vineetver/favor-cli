//! Sparse AI-REML backend.
//!
//! When all loaded kinships are sparse and dense Σ⁻¹ would blow the memory
//! budget, this module assembles the sparse Σ via [`assembler::Assembler`],
//! factors it via faer's sparse Cholesky, and hands the factor to the shared
//! [`crate::staar::kinship::reml`] loop wrapped in
//! [`SigmaSolver::Sparse`](crate::staar::kinship::reml::SigmaSolver::Sparse).
//!
//! upstream:
//!   - per-iteration AI step → `R/glmmkin.R::R_fitglmm_ai:662-710` (the
//!     sparse-aware variant). Upstream still materializes Σ⁻¹ as a dense
//!     matrix because R's `chol2inv` does so even on sparse Cholesky output;
//!     it then uses `sum(Sigma_i * kins[[i]])` to exploit kinship sparsity
//!     in the Frobenius product. We do not materialize Σ⁻¹: the trace term
//!     comes from a stochastic Hutchinson estimator until the Takahashi
//!     selected inversion (#26 #27) lands and restores 1:1 correspondence.
//!   - convergence + boundary refit → shared in `reml.rs::run_reml`.
//!
//! See `UPSTREAM.md` "sparse path is not bit-identical to upstream" for
//! the full story.

pub mod assembler;
pub mod hutchinson;

pub use assembler::Assembler;
pub use hutchinson::{
    sparse_matvec, trace_k_estimate, DEFAULT_HUTCHINSON_PROBES, HUTCHINSON_SEED,
};

use faer::sparse::linalg::solvers::Llt;
use faer::Mat;

use crate::error::FavorError;
use crate::staar::kinship::reml::{
    run_reml, weights_or_ones, SigmaSolver, SolverBuilder, SparseSolverState,
};
use crate::staar::kinship::types::{
    GroupPartition, KinshipMatrix, KinshipState, VarianceComponents,
};

/// Owned wrapper around the faer sparse Cholesky factor of Σ. Cloned into
/// `KinshipState::inverse` (the `Sparse` arm) by [`run_reml`] when the
/// fit takes the sparse path. Once held, every `Σ⁻¹·v` operation in
/// scoring becomes `factor.solve_in_place(v)`.
#[derive(Clone)]
pub struct SparseFactor {
    pub llt: Llt<u32, f64>,
}

/// Builds a sparse [`SigmaSolver`] each AI iteration via the cached
/// assembler. Holds references to all the things the assembler needs
/// (kinships, groups, weights, the row→group lookup) plus the probe
/// count for the trace estimators.
pub struct SparseBuilder<'a> {
    pub groups: &'a GroupPartition,
    pub weights: &'a [f64],
    pub row_to_group: &'a [usize],
    pub assembler: &'a Assembler,
    pub n_probes: usize,
}

impl<'a> SolverBuilder for SparseBuilder<'a> {
    fn build(&self, tau: &VarianceComponents) -> Result<SigmaSolver, FavorError> {
        let factor = self.assembler.factor(
            tau,
            self.groups,
            self.weights,
            self.row_to_group,
        )?;
        Ok(SigmaSolver::Sparse(SparseSolverState {
            factor,
            n: self.assembler.n(),
            n_probes: self.n_probes,
        }))
    }
}

/// Per-backend trace term for `tr(P K_l)` on the sparse path.
///
/// Computed as `tr(Σ⁻¹ K_l) − tr(Σ⁻¹X cov X' Σ⁻¹ K_l)`, with the first
/// term coming from a stochastic Hutchinson estimator and the second
/// term computed exactly via `B' K_l B` where `B = Σ⁻¹X`. Replaced by an
/// exact computation once the Takahashi recursion lands (#26 #27).
///
/// upstream: `R/glmmkin.R::R_fitglmm_ai:698` —
/// `score[i] <- sum(Y * PAPY) - (sum(Sigma_i*kins[[i]]) - sum(Sigma_iX * crossprod(kins[[i]], Sigma_iXcov)))`.
/// Upstream computes `sum(Sigma_i * kins[[i]])` exactly via the dense
/// inverse it materializes; we approximate that term until Takahashi.
pub fn trace_p_k_sparse(
    state: &SparseSolverState,
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

/// Sparse fit_reml entry point. Builds a [`SparseBuilder`] and hands it
/// to the shared [`run_reml`] loop. Caller is `fit_reml` in `mod.rs`
/// after it has decided to take the sparse path.
pub fn fit_reml_sparse(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: Option<&[f64]>,
    init_tau: VarianceComponents,
    n_probes: usize,
) -> Result<KinshipState, FavorError> {
    let n = y.nrows();
    let w_vec = weights_or_ones(n, weights);

    // Row → group lookup, used by the assembler to layer the residual
    // term onto the diagonal entries. One pass over the partition.
    let mut row_to_group = vec![0_usize; n];
    for gi in 0..groups.n_groups() {
        for &row in groups.group(gi) {
            row_to_group[row as usize] = gi;
        }
    }

    let assembler = Assembler::new(n, kinships)?;
    let builder = SparseBuilder {
        groups,
        weights: &w_vec,
        row_to_group: &row_to_group,
        assembler: &assembler,
        n_probes,
    };

    run_reml(y, x, kinships, groups, &w_vec, init_tau, &builder)
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

        let dense_state = fit_reml_dense(
            &y,
            &x,
            std::slice::from_ref(&kin_dense),
            &groups,
            &weights,
            warm_init.clone(),
        )
        .expect("dense fit");

        let sparse_state = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin_sparse),
            &groups,
            None,
            warm_init,
            200,
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
    fn sparse_fit_reml_reproducible() {
        let kin = pedigree_kinship(15, 4);
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let y = synthetic_pheno(n);
        let x = intercept(n);

        let mut warm_init = VarianceComponents::zeros(1, 1);
        warm_init.set_kinship(0, 0.5);
        warm_init.set_group(0, 0.5);

        let s1 = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin),
            &groups,
            None,
            warm_init.clone(),
            64,
        )
        .expect("first sparse fit");
        let s2 = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin),
            &groups,
            None,
            warm_init,
            64,
        )
        .expect("second sparse fit");

        assert_eq!(s1.tau.kinship(0), s2.tau.kinship(0));
        assert_eq!(s1.tau.group(0), s2.tau.group(0));
        assert_eq!(s1.n_iter, s2.n_iter);
    }
}
