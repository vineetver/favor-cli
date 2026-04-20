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
#[cfg(test)]
pub use budget::DEFAULT_KINSHIP_MEM_BYTES;
pub use load::{load_groups, load_kinship, load_random_slope};
pub use pql::fit_pql_glmm;
pub use types::{
    CovarianceIdx, GroupPartition, KinshipInverse, KinshipMatrix, KinshipState,
    VarianceComponents,
};

use faer::Mat;

use crate::error::CohortError;

/// Fit AI-REML for the variance components in
/// `[τ_kinship_1..τ_L, τ_group_1..τ_G]`. Dispatches between the dense
/// kernel (small n, bit-identical to GMMAT) and the sparse kernel
/// (large n + sparse kinships, exact via the Takahashi recursion).
///
/// `budget_bytes` is the dense AI-REML cap, sourced from
/// `Resources::kinship_budget_bytes`. Returns `CohortError::Resource` if
/// neither path can run within it.
pub fn fit_reml(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: Option<&[f64]>,
    budget_bytes: u64,
) -> Result<KinshipState, CohortError> {
    let n = y.nrows();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    if n_comp == 0 {
        return Err(CohortError::Input(
            "fit_reml requires at least one variance component (kinship or group)".into(),
        ));
    }
    if groups.n_samples() != n {
        return Err(CohortError::Input(format!(
            "group partition covers {} samples but y has {n} rows",
            groups.n_samples()
        )));
    }
    for k in kinships {
        if k.n() != n {
            return Err(CohortError::Input(format!(
                "kinship matrix '{}' has shape ({}, {}) but n_samples = {n}",
                k.label(),
                k.n(),
                k.n()
            )));
        }
    }

    let any_sparse = kinships.iter().any(|k| k.is_sparse());
    let all_sparse = !kinships.is_empty() && kinships.iter().all(|k| k.is_sparse());
    let dense_fits = dense_path_fits(n, budget_bytes);

    if !dense_fits && !all_sparse {
        // Dense doesn't fit AND we have at least one inherently dense
        // kinship. Bail out clearly.
        check_memory_budget(n, budget_bytes)?;
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
            budget_bytes,
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
        check_memory_budget(n, budget_bytes)?;
        return dense::fit_reml_dense(y, x, &promoted, groups, &w_vec, init_tau, budget_bytes);
    }

    check_memory_budget(n, budget_bytes)?;
    dense::fit_reml_dense(y, x, kinships, groups, &w_vec, init_tau, budget_bytes)
}

/// Expand an input kinship list for random-slope longitudinal LMM.
///
/// Ports `GMMAT::glmmkin.R:174-183`. For each input `Φ_l` and a
/// per-sample time vector `t`, emit three entries:
///
/// * intercept variance `Φ_l` — left in place (same matrix as input).
/// * intercept-slope covariance `Φ_l · diag(t) + diag(t) · Φ_l`.
/// * slope variance `Φ_l .* (t · tᵀ)` — element-wise product of `Φ_l`
///   with the outer product of `t` with itself.
///
/// Output order is `[int_1..int_L, cov_1..cov_L, slope_1..slope_L]`, so
/// a caller building the `VarianceComponents` layout gets the natural
/// (kinship-major, purpose-minor) ordering. This matches GMMAT's
/// `kins`/`kins[[q+i]]`/`kins[[2*q+i]]` indexing at upstream
/// `glmmkin.R:176-178`.
///
/// Both input variants (Dense / Sparse) are handled; output is always
/// a dense matrix since the column-scaled Kronecker structure destroys
/// the input sparsity pattern (the scaled matrix has `2·nnz(Φ)` non-
/// zeros at most, but the structure changes per entry).
pub fn expand_for_random_slope(
    kinships: &[KinshipMatrix],
    time_var: &[f64],
) -> Result<Vec<KinshipMatrix>, CohortError> {
    let l = kinships.len();
    if l == 0 {
        return Err(CohortError::Input(
            "random.slope requires at least one kinship matrix — GMMAT `glmmkin.R:64` \
             drops the argument when kinship is NULL and samples are unrelated"
                .into(),
        ));
    }
    let n = time_var.len();
    for k in kinships {
        if k.n() != n {
            return Err(CohortError::Input(format!(
                "random.slope time vector has length {} but kinship '{}' has n={}",
                n,
                k.label(),
                k.n(),
            )));
        }
    }

    let to_dense = |k: &KinshipMatrix| -> Mat<f64> {
        match k {
            KinshipMatrix::Dense { matrix, .. } => matrix.clone(),
            KinshipMatrix::Sparse { matrix, .. } => {
                let mut d = Mat::<f64>::zeros(n, n);
                let cp = matrix.symbolic().col_ptr();
                let ri = matrix.symbolic().row_idx();
                let val = matrix.val();
                for c in 0..n {
                    let s = cp[c] as usize;
                    let e = cp[c + 1] as usize;
                    for kk in s..e {
                        let r = ri[kk] as usize;
                        d[(r, c)] = val[kk];
                    }
                }
                d
            }
        }
    };

    let mut intercepts = Vec::with_capacity(l);
    let mut covariances = Vec::with_capacity(l);
    let mut slopes = Vec::with_capacity(l);

    for (li, kin) in kinships.iter().enumerate() {
        let k_dense = to_dense(kin);
        let label = kin.label().to_string();
        intercepts.push(KinshipMatrix::Dense {
            matrix: k_dense.clone(),
            label: format!("{label}_int_{li}"),
        });

        // cov = Φ · diag(t) + diag(t) · Φ
        //      (i, j)  = Φ[i, j] · (t_i + t_j)
        let mut cov = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                cov[(i, j)] = k_dense[(i, j)] * (time_var[i] + time_var[j]);
            }
        }
        covariances.push(KinshipMatrix::Dense {
            matrix: cov,
            label: format!("{label}_cov_{li}"),
        });

        // slope = Φ .* (t · tᵀ)
        //        (i, j) = Φ[i, j] · t_i · t_j
        let mut slope = Mat::<f64>::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                slope[(i, j)] = k_dense[(i, j)] * time_var[i] * time_var[j];
            }
        }
        slopes.push(KinshipMatrix::Dense {
            matrix: slope,
            label: format!("{label}_slope_{li}"),
        });
    }

    let mut out = Vec::with_capacity(3 * l);
    out.extend(intercepts);
    out.extend(covariances);
    out.extend(slopes);
    Ok(out)
}

/// Fit AI-REML with random slopes on the time variable. See
/// [`expand_for_random_slope`] for the variance-component structure.
///
/// Input layout (caller-supplied):
/// * `kinships` — `L` original kinship matrices `Φ_1..Φ_L`.
/// * `groups` — residual group partition (pass `GroupPartition::single`
///   if no per-group heteroscedasticity is wanted).
/// * `time_var` — per-sample time values. GMMAT's semantic is "time
///   since baseline at the measurement"; any numeric column works.
///
/// Output layout (on the returned `KinshipState.tau`):
/// * indices `0..L` — `τ` for the intercept variance of each input `Φ`.
/// * indices `L..2L` — `τ` for the intercept-slope covariance of each `Φ`.
/// * indices `2L..3L` — `τ` for the slope variance of each `Φ`.
/// * indices `3L..3L+G` — per-group residual variances (standard groups).
///
/// The heritability field `state.h2` holds per-input-`Φ` *intercept*
/// heritability (`τ[int_l] / Σ τ`) — it is not meaningful for the cov
/// or slope components, which are not standalone variance fractions.
/// Cross-check: at `time_var = 0`, this function reduces to `fit_reml`
/// with the input kinship list (the cov and slope entries collapse to
/// zero matrices and AI-REML pins their `τ` at zero via boundary refit).
pub fn fit_reml_random_slope(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    time_var: &[f64],
    budget_bytes: u64,
) -> Result<KinshipState, CohortError> {
    if kinships.is_empty() {
        return Err(CohortError::Input(
            "random.slope requires at least one kinship matrix".into(),
        ));
    }
    let n = y.nrows();
    if time_var.len() != n {
        return Err(CohortError::Input(format!(
            "random.slope time vector length {} does not match phenotype n={n}",
            time_var.len(),
        )));
    }
    if groups.n_samples() != n {
        return Err(CohortError::Input(format!(
            "group partition covers {} samples but y has {n} rows",
            groups.n_samples()
        )));
    }

    let expanded = expand_for_random_slope(kinships, time_var)?;
    let l = kinships.len();
    let g = groups.n_groups();
    let triples: Vec<CovarianceIdx> = (0..l)
        .map(|li| CovarianceIdx {
            cov_idx: l + li,
            var_int_idx: li,
            var_slope_idx: 2 * l + li,
        })
        .collect();

    // Warm start: variance entries at var(Y)/(3L+G), covariance entries
    // at 0 (per `glmmkin.R:446-447` which sets diagonal of the
    // Kronecker-flattened (k×k) to var/q and off-diagonal to 0).
    let y_mean: f64 = (0..n).map(|i| y[(i, 0)]).sum::<f64>() / n as f64;
    let y_var: f64 = (0..n).map(|i| (y[(i, 0)] - y_mean).powi(2)).sum::<f64>()
        / (n as f64 - 1.0).max(1.0);
    let n_comp = 3 * l + g;
    let warm = (y_var / n_comp as f64).max(1e-6);
    let mut init_tau = VarianceComponents::zeros(3 * l, g);
    for (li, kin) in kinships.iter().enumerate() {
        let mean_diag = kin.mean_diagonal().max(1e-6);
        init_tau.set_kinship(li, warm / mean_diag); // intercept
        init_tau.set_kinship(l + li, 0.0); // covariance
        init_tau.set_kinship(2 * l + li, warm / mean_diag); // slope
    }
    for gi in 0..g {
        init_tau.set_group(gi, warm);
    }

    let weights = vec![1.0; n];
    let state = reml::run_reml_constrained(
        y,
        x,
        &expanded,
        groups,
        &weights,
        init_tau,
        &triples,
        budget_bytes,
    )?;
    Ok(state)
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
        let state = fit_reml(&y, &x, &[], &groups, None, DEFAULT_KINSHIP_MEM_BYTES).expect("fit_reml groups-only");
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
        let state = fit_reml(&y, &x, &[], &groups, None, DEFAULT_KINSHIP_MEM_BYTES).expect("fit_reml 3 groups");
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
            fit_reml(&y, &x, std::slice::from_ref(&kin), &groups, None, DEFAULT_KINSHIP_MEM_BYTES).expect("fit_reml");
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
        match fit_reml(&y, &x, &[], &groups, None, DEFAULT_KINSHIP_MEM_BYTES) {
            Err(CohortError::Input(msg)) => {
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
        match fit_reml(&y, &x, std::slice::from_ref(&kin), &groups, None, DEFAULT_KINSHIP_MEM_BYTES) {
            Err(CohortError::Input(msg)) => {
                assert!(msg.contains("n_samples"), "msg = {msg}");
                assert!(msg.contains("wrong-size"), "msg = {msg}");
            }
            Err(other) => panic!("expected Input, got {other:?}"),
            Ok(_) => panic!("expected Input error, got Ok"),
        }
    }

    #[test]
    fn expand_for_random_slope_shapes_and_symmetries() {
        let n = 4;
        let kin = {
            let mut m = Mat::<f64>::zeros(n, n);
            for i in 0..n {
                m[(i, i)] = 1.0;
            }
            m[(0, 1)] = 0.5;
            m[(1, 0)] = 0.5;
            KinshipMatrix::Dense {
                matrix: m,
                label: "k".into(),
            }
        };
        let t = vec![1.0, 2.0, 3.0, 4.0];
        let expanded = expand_for_random_slope(std::slice::from_ref(&kin), &t).unwrap();
        // 3L = 3 entries (intercept, cov, slope) for one input kinship.
        assert_eq!(expanded.len(), 3);

        // Intercept is the input kinship untouched.
        let int_mat = expanded[0].as_dense().unwrap();
        assert!((int_mat[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((int_mat[(0, 1)] - 0.5).abs() < 1e-12);

        // Covariance at (0,1) is K[0,1] · (t_0 + t_1) = 0.5 · 3 = 1.5.
        let cov_mat = expanded[1].as_dense().unwrap();
        assert!((cov_mat[(0, 1)] - 1.5).abs() < 1e-12);
        // Symmetric: cov[1,0] = cov[0,1].
        assert!((cov_mat[(1, 0)] - 1.5).abs() < 1e-12);

        // Slope at (0,1) is K[0,1] · t_0 · t_1 = 0.5 · 2 = 1.0.
        let slope_mat = expanded[2].as_dense().unwrap();
        assert!((slope_mat[(0, 1)] - 1.0).abs() < 1e-12);
        // Diagonal slope[i,i] = K[i,i] · t_i² = 1 · t_i².
        for i in 0..n {
            assert!((slope_mat[(i, i)] - t[i].powi(2)).abs() < 1e-12);
        }
    }

    /// Random-slope convergence smoke: repeated-measures data with a
    /// within-subject time trend should fit to convergence and return
    /// a PSD random-effects matrix. No upstream fixture yet — this is
    /// a feasibility + structural check, not a numerical parity check.
    #[test]
    fn fit_reml_random_slope_converges_on_longitudinal_data() {
        // 10 subjects × 3 repeated measurements each = n=30.
        let n_subj = 10;
        let reps = 3;
        let n = n_subj * reps;
        // Subject-level kinship: identity at the row level is the
        // "repeated-measures" kinship that GMMAT adds when duplicated
        // ids are detected (see `glmmkin.R:55-62`). For our smoke test
        // we use a block-diagonal kinship where within-subject blocks
        // are 1 and between-subject off-diagonal is 0 — i.e., repeated
        // measurements share the same random intercept.
        let mut k = Mat::<f64>::zeros(n, n);
        for s in 0..n_subj {
            for r1 in 0..reps {
                for r2 in 0..reps {
                    k[(s * reps + r1, s * reps + r2)] = 1.0;
                }
            }
        }
        let kin = KinshipMatrix::Dense {
            matrix: k,
            label: "subject".into(),
        };

        // Time values cycle 0, 1, 2 per subject.
        let time_var: Vec<f64> = (0..n).map(|i| (i % reps) as f64).collect();

        // Deterministic phenotype with per-subject random intercept and
        // slope plus residual noise.
        let mut state: u64 = 9999;
        let mut u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut subj_int = vec![0.0; n_subj];
        let mut subj_slope = vec![0.0; n_subj];
        for s in 0..n_subj {
            subj_int[s] = 2.0 * u01() - 1.0;
            subj_slope[s] = (2.0 * u01() - 1.0) * 0.3;
        }
        let mut y = Mat::<f64>::zeros(n, 1);
        let mut x = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            let s = i / reps;
            let t = time_var[i];
            let noise = (2.0 * u01() - 1.0) * 0.5;
            y[(i, 0)] = 1.0 + subj_int[s] + subj_slope[s] * t + noise;
            x[(i, 0)] = 1.0;
        }

        let groups = GroupPartition::single(n);
        let state = fit_reml_random_slope(
            &y,
            &x,
            std::slice::from_ref(&kin),
            &groups,
            &time_var,
            1 << 30,
        )
        .expect("random.slope should converge on this longitudinal fixture");

        // τ layout: [var_int, var_cov, var_slope, var_group].
        assert_eq!(state.tau.n_total(), 4);
        let var_int = state.tau.as_slice()[0];
        let var_cov = state.tau.as_slice()[1];
        let var_slope = state.tau.as_slice()[2];
        let var_res = state.tau.as_slice()[3];

        assert!(var_int.is_finite() && var_int >= 0.0, "var_int: {var_int}");
        assert!(
            var_slope.is_finite() && var_slope >= 0.0,
            "var_slope: {var_slope}",
        );
        assert!(var_res.is_finite() && var_res >= 0.0, "var_res: {var_res}");
        assert!(var_cov.is_finite(), "var_cov: {var_cov}");
        // PSD invariant: |cov| ≤ √(var_int · var_slope).
        let psd_limit = (var_int * var_slope).sqrt();
        assert!(
            var_cov.abs() <= psd_limit * (1.0 + 1e-6),
            "PSD violated: |cov|={} > √(var_int · var_slope)={}",
            var_cov.abs(),
            psd_limit,
        );
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
        let state = fit_pql_glmm(&y, &x, std::slice::from_ref(&kin), &groups, out.as_ref(), DEFAULT_KINSHIP_MEM_BYTES)
            .expect("fit_pql_glmm");
        assert_eq!(state.tau.n_total(), 2);
        assert!(state.tau.as_slice().iter().all(|&t| t.is_finite()));
        assert!(state.p_y.nrows() == n);
    }
}
