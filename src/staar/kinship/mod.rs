//! Kinship-aware AI-REML and PQL for the STAAR null model.
//!
//! Ports `GMMAT::glmmkin` from upstream R
//! (`hanchenphd/GMMAT@473b342 R/glmmkin.R`):
//!  - `glmmkin.fit` (family + boundary refit) → `reml::run_reml`
//!  - `glmmkin.ai`  (AI iteration loop)        → `reml::run_reml`
//!  - `R_fitglmm_ai_dense` (dense AI step)     → `reml::ai_step` + `dense::DenseBuilder`
//!  - `R_fitglmm_ai`       (sparse AI step)    → `reml::ai_step` + `sparse::SparseBuilder`
//!  - PQL outer loop                            → `pql::fit_pql_glmm`
//!
//! Cross-checked against GENESIS `R/runAIREMLgaussian.R` (the AI-REML
//! formulas) and `R/iterateAIREMLworkingY.R` (the PQL outer loop
//! structure we follow). See `UPSTREAM.md` for the pinned SHAs and the
//! function map.
//!
//! ## Module layout
//!
//! - [`types`] — `KinshipMatrix` (Dense | Sparse), `GroupPartition`,
//!   `VarianceComponents`, `KinshipState`, `KinshipInverse`.
//! - [`load`] — TSV → `KinshipMatrix` and phenotype-column → `GroupPartition`.
//! - [`budget`] — Dense memory cap + dense/sparse path selector.
//! - [`reml`] — Shared AI-REML algorithm: `SigmaSolver`, `SolverBuilder`,
//!   `ai_step`, `run_reml`. Read this file to understand the math; both
//!   backends use it.
//! - [`dense`] — Dense backend: `DenseBuilder`, `assemble_sigma`,
//!   `trace_p_k_dense`, `fit_reml_dense` entry point.
//! - [`sparse`] — Sparse backend: `TakahashiBuilder` (default),
//!   `HutchinsonBuilder` (fallback), `fit_reml_sparse` entry point.
//!   Submodules `assembler` (cached pattern + symbolic Cholesky),
//!   `takahashi` (exact selected inversion), `hutchinson` (stochastic
//!   estimators kept as fallback).
//! - [`pql`] — Binary GLMM PQL outer loop wrapping `fit_reml`.
//!
//! ## Path selection
//!
//! `fit_reml` picks dense vs sparse at entry:
//! - All kinships dense → dense path always.
//! - At least one kinship sparse AND dense working set fits → dense
//!   (cheap and exact, sparse kinships are pre-promoted).
//! - At least one kinship sparse AND dense working set exceeds budget →
//!   sparse path with the Takahashi backend (1:1 with upstream R).
//!
//! Statistical outputs (h², α, p-values) are identical across paths to
//! within numerical noise. The dense path is bit-identical to GMMAT.

pub mod budget;
pub mod dense;
pub mod load;
pub mod pql;
pub mod reml;
pub mod sparse;
pub mod types;

pub use budget::{check_memory_budget, dense_path_fits};
pub use load::{load_groups, load_kinship};
pub use pql::fit_pql_glmm;
pub use types::{
    GroupPartition, KinshipInverse, KinshipMatrix, KinshipState, VarianceComponents,
};

use faer::Mat;

use crate::error::FavorError;

/// Fit AI-REML for the variance components in
/// `[τ_kinship_1..τ_L, τ_group_1..τ_G]`. Dispatches between the dense
/// kernel (small n, bit-identical to GMMAT) and the sparse kernel
/// (large n + sparse kinships, exact via the Takahashi recursion).
///
/// Returns `FavorError::Resource` if neither path can run within the
/// configured memory budget.
pub fn fit_reml(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: Option<&[f64]>,
) -> Result<KinshipState, FavorError> {
    let n = y.nrows();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    if n_comp == 0 {
        return Err(FavorError::Input(
            "fit_reml requires at least one variance component (kinship or group)".into(),
        ));
    }
    if groups.n_samples() != n {
        return Err(FavorError::Input(format!(
            "group partition covers {} samples but y has {n} rows",
            groups.n_samples()
        )));
    }
    for k in kinships {
        if k.n() != n {
            return Err(FavorError::Input(format!(
                "kinship matrix '{}' has shape ({}, {}) but n_samples = {n}",
                k.label(),
                k.n(),
                k.n()
            )));
        }
    }

    let any_sparse = kinships.iter().any(|k| k.is_sparse());
    let all_sparse = !kinships.is_empty() && kinships.iter().all(|k| k.is_sparse());
    let dense_fits = dense_path_fits(n);

    if !dense_fits && !all_sparse {
        // Dense doesn't fit AND we have at least one inherently dense
        // kinship. Bail out clearly.
        check_memory_budget(n)?;
        unreachable!("check_memory_budget should have errored above");
    }

    let w_vec = reml::weights_or_ones(n, weights);

    // Warm start: var(Y)/(L+G) per slot, kinship slots rescaled by
    // mean(diag(K_l)). Matches `glmmkin.R:309-316`.
    let y_mean = (0..n).map(|i| y[(i, 0)]).sum::<f64>() / n as f64;
    let y_var = (0..n)
        .map(|i| (y[(i, 0)] - y_mean).powi(2))
        .sum::<f64>()
        / (n as f64 - 1.0).max(1.0);
    let mut init_tau = VarianceComponents::zeros(l, g);
    let warm = y_var / n_comp as f64;
    for (li, kin) in kinships.iter().enumerate() {
        let mean_diag = kin.mean_diagonal();
        let scale = if mean_diag > 0.0 { mean_diag } else { 1.0 };
        init_tau.set_kinship(li, warm / scale);
    }
    for gi in 0..g {
        init_tau.set_group(gi, warm);
    }

    // Path selection. Sparse only when dense won't fit. Otherwise dense —
    // it's exact and slightly faster on small studies. Sparse path uses
    // Takahashi (exact, 1:1 with upstream R) by default; Hutchinson
    // remains available as a fallback for tests.
    let take_sparse = !dense_fits && all_sparse;
    if take_sparse {
        return sparse::fit_reml_sparse(
            y,
            x,
            kinships,
            groups,
            weights,
            init_tau,
            sparse::SparseSolverKind::default(),
            None,
        );
    }

    if any_sparse && dense_fits {
        // Dense fits but kinships are sparse-stored — promote them once.
        // Cheap when dense fits.
        let promoted: Vec<KinshipMatrix> = kinships
            .iter()
            .map(|k| match k.as_sparse() {
                Some(sp) => {
                    let n = sp.nrows();
                    let mut dense = Mat::<f64>::zeros(n, n);
                    let col_ptr = sp.symbolic().col_ptr();
                    let row_idx = sp.symbolic().row_idx();
                    let val = sp.val();
                    for j in 0..n {
                        let start = col_ptr[j] as usize;
                        let end = col_ptr[j + 1] as usize;
                        for k in start..end {
                            let i = row_idx[k] as usize;
                            dense[(i, j)] = val[k];
                        }
                    }
                    KinshipMatrix::Dense {
                        matrix: dense,
                        label: k.label().to_string(),
                    }
                }
                None => k.clone(),
            })
            .collect();
        check_memory_budget(n)?;
        return dense::fit_reml_dense(y, x, &promoted, groups, &w_vec, init_tau);
    }

    check_memory_budget(n)?;
    dense::fit_reml_dense(y, x, kinships, groups, &w_vec, init_tau)
}

#[cfg(test)]
mod tests {
    //! Integration tests across the public dispatcher. Per-component tests
    //! live in their own files (`types.rs`, `dense.rs`, `budget.rs`,
    //! `sparse/mod.rs`).

    use super::*;
    use crate::output::{create, Output, OutputMode};
    use faer::Mat;

    fn make_x_intercept(n: usize) -> Mat<f64> {
        let mut x = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }
        x
    }

    fn block_family_kinship(n_fam: usize, family_size: usize, label: &str) -> KinshipMatrix {
        let n = n_fam * family_size;
        let mut k = Mat::<f64>::zeros(n, n);
        for f in 0..n_fam {
            let base = f * family_size;
            for i in 0..family_size {
                for j in 0..family_size {
                    k[(base + i, base + j)] = if i == j { 1.0 } else { 0.5 };
                }
            }
        }
        KinshipMatrix::new(k, label.into()).unwrap()
    }

    fn null_output() -> Box<dyn Output> {
        create(&OutputMode::Machine)
    }

    #[test]
    fn fit_reml_groups_only_matches_ols_sigma2() {
        use crate::staar::model::fit_glm;
        let n = 150;
        let mut y = Mat::<f64>::zeros(n, 1);
        let mut x = Mat::<f64>::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            let xi = (i as f64 / n as f64) * 4.0 - 2.0;
            x[(i, 1)] = xi;
            y[(i, 0)] = 0.7 * xi + 0.3 * (i as f64 * 0.41).sin();
        }
        let glm = fit_glm(&y, &x);
        let groups = GroupPartition::single(n);
        let state = fit_reml(&y, &x, &[], &groups, None).expect("fit_reml groups-only");
        assert!(state.n_iter < 500);
        let rel_err = (state.tau.group(0) - glm.sigma2).abs() / glm.sigma2;
        assert!(
            rel_err < 1e-3,
            "groups-only τ_g0 = {:.6} vs OLS σ² = {:.6} (rel err {:.3e})",
            state.tau.group(0),
            glm.sigma2,
            rel_err
        );
    }

    #[test]
    fn fit_reml_three_groups_recovers_per_group_variance() {
        let n = 300;
        let mut y = Mat::<f64>::zeros(n, 1);
        let x = make_x_intercept(n);
        let scales = [0.5_f64, 1.5, 3.0];
        let mut assignments = vec![0_usize; n];
        for i in 0..n {
            let g = i % 3;
            assignments[i] = g;
            let r = ((i as f64 * 1.13).sin() + (i as f64 * 0.41).cos()) * scales[g].sqrt();
            y[(i, 0)] = r;
        }
        let group_labels: Vec<String> = vec!["A".into(), "B".into(), "C".into()];
        let groups = GroupPartition::from_assignments(&assignments, &group_labels).unwrap();
        let state = fit_reml(&y, &x, &[], &groups, None).expect("fit_reml 3 groups");
        assert!(state.n_iter < 500);
        assert!(state.tau.group(0) < state.tau.group(1));
        assert!(state.tau.group(1) < state.tau.group(2));
    }

    #[test]
    fn fit_reml_kinship_block_family_converges() {
        let kin = block_family_kinship(20, 10, "fam");
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let mut y = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            let fam = i / 10;
            let g = ((fam as f64 * 0.91).sin()) * 0.6;
            let e = (i as f64 * 1.7).sin() * 0.4;
            y[(i, 0)] = g + e;
        }
        let x = make_x_intercept(n);
        let state =
            fit_reml(&y, &x, std::slice::from_ref(&kin), &groups, None).expect("fit_reml");
        assert!(state.n_iter < 500);
        let total: f64 = state.tau.as_slice().iter().sum();
        assert!(total > 0.0 && total.is_finite());
        assert!(state.tau.kinship(0) >= 0.0);
        assert!(state.tau.group(0) >= 0.0);
    }

    #[test]
    fn fit_reml_rejects_group_partition_size_mismatch() {
        // GroupPartition covers more samples than y has rows.
        let n = 10;
        let y = Mat::<f64>::zeros(n, 1);
        let x = make_x_intercept(n);
        let groups = GroupPartition::single(n + 5);
        match fit_reml(&y, &x, &[], &groups, None) {
            Err(FavorError::Input(msg)) => {
                assert!(msg.contains("group partition covers"), "msg = {msg}");
            }
            Err(other) => panic!("expected Input, got {other:?}"),
            Ok(_) => panic!("expected Input error, got Ok"),
        }
    }

    #[test]
    fn fit_reml_rejects_kinship_size_mismatch() {
        let n = 10;
        let y = Mat::<f64>::zeros(n, 1);
        let x = make_x_intercept(n);
        let groups = GroupPartition::single(n);
        // Build a kinship over n+5 samples — must be rejected.
        let kin = block_family_kinship(3, 5, "wrong-size");
        assert_ne!(kin.n(), n);
        match fit_reml(&y, &x, std::slice::from_ref(&kin), &groups, None) {
            Err(FavorError::Input(msg)) => {
                assert!(msg.contains("n_samples"), "msg = {msg}");
                assert!(msg.contains("wrong-size"), "msg = {msg}");
            }
            Err(other) => panic!("expected Input, got {other:?}"),
            Ok(_) => panic!("expected Input error, got Ok"),
        }
    }

    #[test]
    fn fit_pql_glmm_runs_and_returns_state() {
        let kin = block_family_kinship(20, 10, "fam");
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let mut y = Mat::<f64>::zeros(n, 1);
        let mut x = Mat::<f64>::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            let xi = (i as f64 / n as f64) * 2.0 - 1.0;
            x[(i, 1)] = xi;
            let p = 1.0 / (1.0 + (-(0.5 * xi)).exp());
            let u = 0.5 + 0.5 * ((i as f64 * 1.234567).sin() * 0.999);
            y[(i, 0)] = if u < p { 1.0 } else { 0.0 };
        }
        let out = null_output();
        let state = fit_pql_glmm(&y, &x, std::slice::from_ref(&kin), &groups, out.as_ref())
            .expect("fit_pql_glmm");
        assert_eq!(state.tau.n_total(), 2);
        assert!(state.tau.as_slice().iter().all(|&t| t.is_finite()));
        assert!(state.p_y.nrows() == n);
    }
}
