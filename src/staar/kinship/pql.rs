//! Binary GLMM via PQL (`R/glmmkin.R:326` outer loop).
//!
//! Wraps `fit_reml` on the working pseudo-response `Y_w = η + (y − μ)/μ'`
//! with IRLS weights `W = μ(1 − μ)`. Each PQL iteration calls back into
//! the dense or sparse `fit_reml` dispatcher with the current weights;
//! the inner kernel is unchanged. Convergence is on max |Δη| between
//! iterations. Failure to converge is a warning, not an error — the
//! caller still gets the most recent state and is told to treat it with
//! caution.

use faer::Mat;

use crate::error::FavorError;
use crate::output::Output;
use crate::staar::kinship::budget::check_memory_budget;
use crate::staar::kinship::fit_reml;
use crate::staar::kinship::types::{GroupPartition, KinshipMatrix, KinshipState};

const PQL_MAX_ITER: usize = 100;
const PQL_TOL: f64 = 1e-5;

pub fn fit_pql_glmm(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    out: &dyn Output,
) -> Result<KinshipState, FavorError> {
    use crate::staar::model::fit_logistic;

    let n = y.nrows();
    // Skip the dense memory check on the sparse path — sparse fit_reml has
    // its own budget logic.
    if kinships.iter().all(|k| !k.is_sparse()) {
        check_memory_budget(n)?;
    }

    let init = fit_logistic(y, x, 25);
    let Some(init_mu) = init.fitted_values.as_ref() else {
        return Err(FavorError::Analysis(
            "fit_logistic returned without fitted values — cannot start PQL".into(),
        ));
    };
    let mut mu = Mat::<f64>::zeros(n, 1);
    let mut eta = Mat::<f64>::zeros(n, 1);
    for i in 0..n {
        let p = init_mu[i].clamp(1e-6, 1.0 - 1e-6);
        mu[(i, 0)] = p;
        eta[(i, 0)] = (p / (1.0 - p)).ln();
    }

    // Row → group lookup for BLUP recovery and per-group residual scaling.
    let mut row_to_group = vec![0_usize; n];
    for gi in 0..groups.n_groups() {
        for &row in groups.group(gi) {
            row_to_group[row as usize] = gi;
        }
    }

    let mut last_state: Option<KinshipState> = None;
    let mut converged = false;
    let mut last_diff = f64::INFINITY;
    let mut pql_iters = 0;
    for iter in 0..PQL_MAX_ITER {
        pql_iters = iter + 1;
        let mut y_work = Mat::<f64>::zeros(n, 1);
        let mut w_irls = vec![0.0; n];
        for i in 0..n {
            // GMMAT clamps μ at 1e-5/(1-1e-5) (`R/glmmkin.R::glmmkin.fit`).
            let mu_i = mu[(i, 0)].clamp(1e-5, 1.0 - 1e-5);
            let w_i = mu_i * (1.0 - mu_i);
            w_irls[i] = w_i;
            y_work[(i, 0)] = eta[(i, 0)] + (y[(i, 0)] - mu_i) / w_i;
        }

        let state = fit_reml(&y_work, x, kinships, groups, Some(&w_irls))?;

        // BLUP: η_new = Xα + D PY where D = Σ - R, so η_new = Y_w - R PY.
        // R is diagonal with R[i,i] = τ_{g(i)} / W[i].
        let mut new_eta = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            let tau_g = state.tau.group(row_to_group[i]);
            new_eta[(i, 0)] = y_work[(i, 0)] - tau_g * state.p_y[(i, 0)] / w_irls[i];
        }

        let mut max_diff = 0.0_f64;
        for i in 0..n {
            let d = (new_eta[(i, 0)] - eta[(i, 0)]).abs();
            if d > max_diff {
                max_diff = d;
            }
        }
        last_diff = max_diff;
        eta = new_eta;
        for i in 0..n {
            let p = 1.0 / (1.0 + (-eta[(i, 0)].clamp(-500.0, 500.0)).exp());
            mu[(i, 0)] = p;
        }

        last_state = Some(state);
        if max_diff < PQL_TOL {
            converged = true;
            break;
        }
    }

    if !converged {
        out.warn(&format!(
            "PQL did not converge after {pql_iters} iterations (last max |Δη| = {last_diff:.3e}). \
             Score test results use the last PQL state; treat with caution."
        ));
    }

    last_state.ok_or_else(|| FavorError::Analysis("PQL outer loop did not run".into()))
}
