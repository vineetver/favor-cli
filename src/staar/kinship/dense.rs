//! Dense AI-REML kernel — the in-memory path used for small-to-medium n.
//!
//! Materializes Σ, Σ⁻¹, Σ⁻¹X as `Mat<f64>` and runs the standard AI loop
//! straight from `R/glmmkin.R::R_fitglmm_ai`. Identical math to upstream
//! GMMAT for `family = gaussian, method.optim = "AI"`. Bounded by the
//! dense memory budget — see `crate::staar::kinship::budget`.
//!
//! Sparse path lives in `sparse.rs`; selection between the two happens in
//! `crate::staar::kinship::fit_reml`.

use faer::prelude::*;
use faer::Mat;

use crate::error::FavorError;
use crate::staar::kinship::types::{
    GroupPartition, KinshipInverse, KinshipMatrix, KinshipState, VarianceComponents,
};

/// REML convergence tolerance, matches `glmmkin` default.
pub(super) const REML_TOL: f64 = 1e-5;
/// Inner AI iteration cap, matches `glmmkin` `maxiter` default.
pub(super) const REML_MAX_ITER: usize = 500;
/// Outer boundary-refit cap. Upstream has no explicit cap; 10 is a hard
/// ceiling against pathological oscillation.
pub(super) const REML_MAX_OUTER_REFITS: usize = 10;
/// Step-halving attempts before giving up on a non-negative τ feasible step.
pub(super) const REML_STEP_HALVE_MAX: usize = 50;
/// Boundary factor for fixing τ at 0 in the outer refit loop.
pub(super) const BOUNDARY_FACTOR: f64 = 1.01;

// ─── Linear algebra helpers ─────────────────────────────────────────────────

pub(super) fn weights_or_one(n: usize, weights: Option<&[f64]>) -> Vec<f64> {
    match weights {
        Some(w) => {
            assert_eq!(w.len(), n);
            w.to_vec()
        }
        None => vec![1.0; n],
    }
}

fn frob_inner(a: &Mat<f64>, b: &Mat<f64>) -> f64 {
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

fn matvec(m: &Mat<f64>, v: &Mat<f64>) -> Mat<f64> {
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

fn solve_spd(a: &Mat<f64>, b: &Mat<f64>) -> Mat<f64> {
    a.col_piv_qr().solve(b)
}

fn invert_dense(a: &Mat<f64>) -> Mat<f64> {
    let n = a.nrows();
    debug_assert_eq!(n, a.ncols());
    let eye = Mat::<f64>::identity(n, n);
    solve_spd(a, &eye)
}

// ─── Σ assembly ─────────────────────────────────────────────────────────────

/// Build dense Σ from current τ:
///
/// ```text
/// diagP[i] = τ_{L + g(i)} / W[i]      for i in group g
/// Σ        = diag(diagP) + Σ_l τ_l · K_l
/// ```
///
/// All kinships must be the dense variant; sparse kinships route through
/// `sparse.rs`.
pub(super) fn assemble_sigma(
    n: usize,
    tau: &VarianceComponents,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
) -> Mat<f64> {
    debug_assert_eq!(tau.n_kinship(), kinships.len());
    debug_assert_eq!(tau.n_group(), groups.n_groups());

    let mut sigma = Mat::<f64>::zeros(n, n);
    for li in 0..kinships.len() {
        let tau_l = tau.kinship(li);
        let k = kinships[li]
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

// ─── AI step ────────────────────────────────────────────────────────────────

/// Per-iteration AI-REML scratch returned by [`ai_step`].
struct AiStep {
    d_tau: Vec<f64>,
    alpha: Mat<f64>,
    p_y: Mat<f64>,
    sigma_inv: Mat<f64>,
    sigma_inv_x: Mat<f64>,
    cov: Mat<f64>,
}

/// One AI-REML iteration step (`R/glmmkin.R::R_fitglmm_ai`).
fn ai_step(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    tau: &VarianceComponents,
    fixtau: &[bool],
) -> AiStep {
    let n = y.nrows();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    debug_assert_eq!(fixtau.len(), n_comp);

    let sigma = assemble_sigma(n, tau, kinships, groups, weights);
    let sigma_inv = invert_dense(&sigma);
    let sigma_inv_x = &sigma_inv * x;
    let xt_sinv_x = x.transpose() * &sigma_inv_x;
    let cov = invert_dense(&xt_sinv_x);

    let sinv_y = matvec(&sigma_inv, y);
    let xt_sinv_y = x.transpose() * &sinv_y;
    let alpha = &cov * &xt_sinv_y;
    let eta = x * &alpha;
    let mut residual = Mat::<f64>::zeros(n, 1);
    for i in 0..n {
        residual[(i, 0)] = y[(i, 0)] - eta[(i, 0)];
    }
    let p_y = matvec(&sigma_inv, &residual);

    // diag(P) where P = Σ⁻¹ − Σ⁻¹X (X'Σ⁻¹X)⁻¹ X'Σ⁻¹.
    let k_cov = sigma_inv_x.ncols();
    let mut diag_p = vec![0.0; n];
    for i in 0..n {
        let mut hi = 0.0;
        for a in 0..k_cov {
            for b in 0..k_cov {
                hi += sigma_inv_x[(i, a)] * cov[(a, b)] * sigma_inv_x[(i, b)];
            }
        }
        diag_p[i] = sigma_inv[(i, i)] - hi;
    }

    // Score:
    //   kinship-l: U_l = (PY)' K_l (PY) − tr(P K_l)
    //   group-g:   U_g = Σ_{i∈g} (PY[i])²/W[i] − Σ_{i∈g} diag_p[i]/W[i]
    let mut score = vec![0.0; n_comp];
    let mut k_py: Vec<Mat<f64>> = Vec::with_capacity(l);
    for li in 0..l {
        let k_mat = kinships[li].as_dense().unwrap();
        let v = matvec(k_mat, &p_y);
        // tr(P K_l) = tr(Σ⁻¹ K_l) − tr(Σ⁻¹X cov X'Σ⁻¹ K_l).
        let sinv_k = &sigma_inv * k_mat;
        let tr_sinv_k = (0..n).map(|i| sinv_k[(i, i)]).sum::<f64>();
        let k_sinv_x = k_mat * &sigma_inv_x;
        let k_sinv_x_cov = &k_sinv_x * &cov;
        let tr_correction = frob_inner(&sigma_inv_x, &k_sinv_x_cov);
        let tr_p_k = tr_sinv_k - tr_correction;

        let mut py_k_py = 0.0;
        for i in 0..n {
            py_k_py += p_y[(i, 0)] * v[(i, 0)];
        }
        score[li] = py_k_py - tr_p_k;
        k_py.push(v);
    }

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

    // AI matrix: AI[i,j] = (PY)' V_i P V_j (PY).
    let mut ai = Mat::<f64>::zeros(n_comp, n_comp);

    let project_through_p = |v: &Mat<f64>| -> Mat<f64> {
        let sinv_v = matvec(&sigma_inv, v);
        let xt_sinv_v = x.transpose() * &sinv_v;
        let cov_xt_sinv_v = &cov * &xt_sinv_v;
        let proj = &sigma_inv_x * &cov_xt_sinv_v;
        let mut out = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            out[(i, 0)] = sinv_v[(i, 0)] - proj[(i, 0)];
        }
        out
    };

    let mut p_k_py: Vec<Mat<f64>> = Vec::with_capacity(l);
    for li in 0..l {
        p_k_py.push(project_through_p(&k_py[li]));
    }

    let mut v_py: Vec<Mat<f64>> = Vec::with_capacity(g);
    let mut p_v_py: Vec<Mat<f64>> = Vec::with_capacity(g);
    for gi in 0..g {
        let mut vpy = Mat::<f64>::zeros(n, 1);
        for &row in groups.group(gi) {
            let i = row as usize;
            vpy[(i, 0)] = p_y[(i, 0)] / weights[i];
        }
        let p_vpy = project_through_p(&vpy);
        v_py.push(vpy);
        p_v_py.push(p_vpy);
    }

    // (l, l'): (K_l PY)' (P K_l' PY)
    for li in 0..l {
        for lj in 0..l {
            let mut s = 0.0;
            for i in 0..n {
                s += k_py[li][(i, 0)] * p_k_py[lj][(i, 0)];
            }
            ai[(li, lj)] = s;
        }
    }

    // (l, g): (K_l PY)' (P V_g PY) — symmetric.
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

    // (g, g'): (V_g PY)' (P V_g' PY).
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

    // Solve AI · Δτ = score for free components only.
    let mut d_tau = vec![0.0; n_comp];
    let n_free = fixtau.iter().filter(|&&f| !f).count();
    if n_free > 0 {
        let mut ai_free = Mat::<f64>::zeros(n_free, n_free);
        let mut score_free = Mat::<f64>::zeros(n_free, 1);
        let idx_map: Vec<usize> = (0..n_comp).filter(|&i| !fixtau[i]).collect();
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
    }

    AiStep {
        d_tau,
        alpha,
        p_y,
        sigma_inv,
        sigma_inv_x,
        cov,
    }
}

// ─── Inner / outer fit loops ────────────────────────────────────────────────

fn fit_reml_inner(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    init_tau: VarianceComponents,
    fixtau: &[bool],
) -> Result<KinshipState, FavorError> {
    let n = y.nrows();
    let k = x.ncols();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    debug_assert_eq!(init_tau.n_total(), n_comp);
    debug_assert_eq!(fixtau.len(), n_comp);

    let mut tau = init_tau;
    let mut alpha_prev = Mat::<f64>::zeros(k, 1);
    let mut iter_used = 0;

    let mut last_p_y = Mat::<f64>::zeros(n, 1);
    let mut last_sigma_inv = Mat::<f64>::zeros(n, n);
    let mut last_sigma_inv_x = Mat::<f64>::zeros(n, k);
    let mut last_cov = Mat::<f64>::zeros(k, k);

    for iter in 0..REML_MAX_ITER {
        iter_used = iter + 1;
        let tau0 = tau.clone();
        let step_out = ai_step(y, x, kinships, groups, weights, &tau, fixtau);

        // Newton step + step-halving on negative τ.
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

        // Convergence: 2·max(rel_α, rel_τ) < tol.
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

        last_p_y = step_out.p_y;
        last_sigma_inv = step_out.sigma_inv;
        last_sigma_inv_x = step_out.sigma_inv_x;
        last_cov = step_out.cov;
        alpha_prev = step_out.alpha;

        if 2.0 * rel_alpha.max(rel_tau) < REML_TOL {
            break;
        }

        // Divergence guard.
        if tau.as_slice().iter().any(|&t| t.abs() > REML_TOL.powi(-2)) {
            return Err(FavorError::Analysis(format!(
                "AI-REML diverged: |τ_max| > {}",
                REML_TOL.powi(-2)
            )));
        }
    }

    let total_var: f64 = tau.as_slice().iter().sum();
    let h2: Vec<f64> = if total_var > 0.0 {
        (0..l).map(|i| tau.kinship(i) / total_var).collect()
    } else {
        vec![0.0; l]
    };

    Ok(KinshipState {
        tau,
        inverse: KinshipInverse::Dense(last_sigma_inv),
        sigma_inv_x: last_sigma_inv_x,
        cov: last_cov,
        p_y: last_p_y,
        h2,
        n_iter: iter_used,
        outer_refits: 0,
    })
}

/// Outer boundary-refit loop. Caller is `crate::staar::kinship::fit_reml`
/// after it has decided to take the dense path.
pub(super) fn fit_reml_dense(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    init_tau: VarianceComponents,
) -> Result<KinshipState, FavorError> {
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;

    let mut fixtau = vec![false; n_comp];
    let mut state = fit_reml_inner(y, x, kinships, groups, weights, init_tau, &fixtau)?;

    for refit_iter in 0..REML_MAX_OUTER_REFITS {
        let mut changed = false;
        for i in 0..n_comp {
            if !fixtau[i] && state.tau.as_slice()[i] < BOUNDARY_FACTOR * REML_TOL {
                fixtau[i] = true;
                changed = true;
            }
        }
        if !changed {
            state.outer_refits = refit_iter;
            return Ok(state);
        }
        let mut tau_refit = state.tau.clone();
        for i in 0..n_comp {
            if fixtau[i] {
                tau_refit.as_slice_mut()[i] = 0.0;
            }
        }
        state = fit_reml_inner(y, x, kinships, groups, weights, tau_refit, &fixtau)?;
        state.outer_refits = refit_iter + 1;
    }

    Err(FavorError::Analysis(format!(
        "AI-REML boundary refit did not stabilize within {REML_MAX_OUTER_REFITS} iterations"
    )))
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
