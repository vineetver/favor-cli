//! Dense AI-REML backend.
//!
//! Materializes Σ as an `n × n` `Mat<f64>` each AI iteration, inverts it,
//! and hands the inverse to the shared [`crate::staar::kinship::reml`] loop
//! via [`SigmaSolver::Dense`]. The dense path is bit-identical to upstream
//! GMMAT for `family = gaussian, method.optim = "AI"`.
//!
//! upstream:
//!   - per-iteration AI step → `R/glmmkin.R::R_fitglmm_ai_dense:712-769`
//!     (and the C++ equivalent `src/fitglmm.cpp::fitglmm_ai:541+`,
//!     used by upstream when `n < 2^15.5`).
//!   - convergence + boundary refit → shared in `reml.rs::run_reml`,
//!     which ports `glmmkin.fit + glmmkin.ai` (`R/glmmkin.R:122-422`).

use faer::prelude::Solve;
use faer::Mat;

use crate::error::FavorError;
use crate::staar::kinship::reml::{
    frob_inner, run_reml, SigmaSolver, SolverBuilder,
};
use crate::staar::kinship::types::{
    GroupPartition, KinshipMatrix, KinshipState, VarianceComponents,
};

/// Σ assembly for the dense backend:
///
/// ```text
/// Σ[i,i] += τ_g(i) / W[i]              (residual block, group g(i))
/// Σ[i,j] += Σ_l τ_l · K_l[i,j]         (kinship blocks)
/// ```
///
/// All kinships must be the dense variant — sparse kinships flow through
/// the sparse path. Matches `R_fitglmm_ai_dense:716-721` lines that build
/// `Sigma <- diag(diagSigma); Sigma <- Sigma + tau[i+ng] * kins[[i]]`.
pub fn assemble_sigma(
    n: usize,
    tau: &VarianceComponents,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
) -> Mat<f64> {
    debug_assert_eq!(tau.n_kinship(), kinships.len());
    debug_assert_eq!(tau.n_group(), groups.n_groups());

    let mut sigma = Mat::<f64>::zeros(n, n);
    for (li, kin) in kinships.iter().enumerate() {
        let tau_l = tau.kinship(li);
        let k = kin
            .as_dense()
            .expect("dense path called with non-dense kinship — fit_reml dispatch bug");
        for r in 0..n {
            for c in 0..n {
                sigma[(r, c)] += tau_l * k[(r, c)];
            }
        }
    }
    for gi in 0..groups.n_groups() {
        let tau_g = tau.group(gi);
        for &row in groups.group(gi) {
            let i = row as usize;
            sigma[(i, i)] += tau_g / weights[i];
        }
    }
    sigma
}

/// Build a dense Σ⁻¹ wrapped in a [`SigmaSolver`]. Each AI iteration calls
/// this once. The shared loop hands the solver to `ai_step`.
pub struct DenseBuilder<'a> {
    pub n: usize,
    pub kinships: &'a [KinshipMatrix],
    pub groups: &'a GroupPartition,
    pub weights: &'a [f64],
}

impl<'a> SolverBuilder for DenseBuilder<'a> {
    fn build(&self, tau: &VarianceComponents) -> Result<SigmaSolver, FavorError> {
        let sigma = assemble_sigma(self.n, tau, self.kinships, self.groups, self.weights);
        let sigma_inv = invert_dense(&sigma);
        Ok(SigmaSolver::Dense(sigma_inv))
    }
}

/// Per-backend trace term for `tr(P K_l)`. Computed as
/// `tr(Σ⁻¹ K_l) − tr(Σ⁻¹X cov X' Σ⁻¹ K_l)` with both terms exact.
///
/// upstream: `R/glmmkin.R::R_fitglmm_ai_dense:755-757` —
/// `score[i] <- sum(Y * PAPY) - sum(P * kins[[idxtau[i]-ng]])`. Upstream
/// folds `tr(P K_l) = sum(P * K_l)` because it has materialized `P`. We
/// don't materialize `P`; we substitute the algebraic identity
/// `tr(P K_l) = tr(Σ⁻¹ K_l) − tr(Σ⁻¹X cov X'Σ⁻¹ K_l)` and compute the
/// two terms separately. Identical math; the bit pattern matches what
/// the previous standalone `dense.rs` produced.
pub fn trace_p_k_dense(
    sigma_inv: &Mat<f64>,
    kinship: &KinshipMatrix,
    sigma_inv_x: &Mat<f64>,
    cov: &Mat<f64>,
) -> f64 {
    let k_mat = kinship
        .as_dense()
        .expect("dense trace_p_k called with non-dense kinship");
    let n = sigma_inv.nrows();

    // tr(Σ⁻¹ K_l) — form the dense product, sum the diagonal.
    let sinv_k = sigma_inv * k_mat;
    let tr_sinv_k = (0..n).map(|i| sinv_k[(i, i)]).sum::<f64>();

    // tr(Σ⁻¹X cov X' Σ⁻¹ K_l) = frob(Σ⁻¹X, K_l Σ⁻¹X cov).
    let k_sinv_x = k_mat * sigma_inv_x;
    let k_sinv_x_cov = &k_sinv_x * cov;
    let tr_correction = frob_inner(sigma_inv_x, &k_sinv_x_cov);

    tr_sinv_k - tr_correction
}

/// Solve `A · x = b` for SPD `A` via column-pivoted QR.
fn solve_spd(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    a.col_piv_qr().solve(b)
}

/// Invert an SPD matrix by solving against the identity. Used by the
/// dense path to materialize Σ⁻¹ each AI iteration.
fn invert_dense(a: &Mat<f64>) -> Mat<f64> {
    let n = a.nrows();
    debug_assert_eq!(n, a.ncols());
    let eye = Mat::<f64>::identity(n, n);
    solve_spd(a, &eye)
}

/// Dense fit_reml entry point. Builds a [`DenseBuilder`] and hands it to
/// the shared [`run_reml`] loop. Caller is `fit_reml` in `mod.rs` after it
/// has decided to take the dense path; sparse kinships are pre-promoted
/// to dense before reaching here.
pub fn fit_reml_dense(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    init_tau: VarianceComponents,
) -> Result<KinshipState, FavorError> {
    let n = y.nrows();
    let builder = DenseBuilder {
        n,
        kinships,
        groups,
        weights,
    };
    run_reml(y, x, kinships, groups, weights, init_tau, &builder)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assemble_sigma_single_kinship_single_group() {
        let n = 4;
        let mut k = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            k[(i, i)] = 1.0;
        }
        k[(0, 1)] = 0.5;
        k[(1, 0)] = 0.5;
        let kin = KinshipMatrix::new(k, "K".into()).unwrap();
        let groups = GroupPartition::single(n);
        let weights = vec![1.0; n];
        let mut tau = VarianceComponents::zeros(1, 1);
        tau.set_kinship(0, 2.0);
        tau.set_group(0, 0.5);

        let sigma = assemble_sigma(n, &tau, std::slice::from_ref(&kin), &groups, &weights);

        assert!((sigma[(0, 0)] - 2.5).abs() < 1e-12);
        assert!((sigma[(0, 1)] - 1.0).abs() < 1e-12);
        assert!((sigma[(1, 0)] - 1.0).abs() < 1e-12);
        assert!((sigma[(2, 2)] - 2.5).abs() < 1e-12);
        assert!((sigma[(2, 3)]).abs() < 1e-12);
    }
}
