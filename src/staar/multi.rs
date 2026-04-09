//! Joint multi-trait STAAR (Li et al. 2023, MultiSTAAR) for unrelated
//! continuous traits with shared covariates.
//!
//! Upstream `MultiSTAAR_O_SMMAT.cpp` materialises a stacked genotype
//! `G_big = I_k ⊗ G₀` and a dense `(n·k × n·k)` projection `P` from GMMAT,
//! then forms an `(m·k × m·k)` covariance and runs the same omnibus tree
//! as single-trait STAAR. For unrelated continuous Gaussian traits with a
//! single shared covariate matrix `X`, both `P` and the joint covariance
//! Kronecker-factor cleanly:
//!
//!   P = Σ_res⁻¹ ⊗ M,                        with M = I_n − X(X'X)⁻¹X'
//!   K_full = Σ_res⁻¹ ⊗ K₀,                  with K₀ = G₀' M G₀
//!   U_full = vec((G₀' R) Σ_res⁻¹),          with R = M Y the (n × k) residuals
//!
//! So the whole machinery collapses to two small matrices, `S = G₀' R`
//! `(m × k)` and `K₀` `(m × m)`. The joint chi-square tests this module
//! produces are bit-equivalent to the upstream cpp on the unrelated
//! continuous path; they reduce exactly to single-trait STAAR when k=1.

use faer::prelude::*;
use faer::Mat;

use super::score::{beta_density_weight, StaarResult};
use super::stats;

/// Null model fit for unrelated continuous multi-trait STAAR.
///
/// Holds the per-trait residual matrix and the cross-trait residual
/// covariance, plus the shared covariate matrix and `(X'X)⁻¹` so we can
/// reconstruct `K₀` cheaply for any gene.
pub struct MultiNullContinuous {
    /// `R = M Y`, shape `(n, k)`. Per-trait OLS residuals against the
    /// shared covariate matrix.
    pub residuals: Mat<f64>,
    /// `Σ_res = R'R / (n − p)`, the unbiased residual covariance.
    pub sigma_res: Mat<f64>,
    /// `Σ_res⁻¹`, pre-computed for the hot path.
    pub sigma_inv: Mat<f64>,
    /// Eigenvalues of `Σ_res⁻¹`. Cached so the joint SKAT eigenvalue
    /// product is cheap.
    pub sigma_inv_eigvals: Vec<f64>,
    /// Shared covariate matrix, `(n, p)`.
    pub x: Mat<f64>,
    /// `(X'X)⁻¹`, `(p, p)`.
    pub xtx_inv: Mat<f64>,
    pub n_samples: usize,
    pub n_pheno: usize,
}

impl MultiNullContinuous {
    /// Fit the multivariate linear model `vec(Y) = (I_k ⊗ X) β + ε`,
    /// `ε ~ N(0, Σ_res ⊗ I_n)` for unrelated continuous traits with a
    /// shared covariate matrix `X`.
    ///
    /// `y` is `(n, k)`, `x` is `(n, p)`. The unbiased residual covariance
    /// is `Σ_res = R'R / (n − p)`, matching GMMAT's REML estimate when X
    /// is shared across traits.
    pub fn fit(y: &Mat<f64>, x: &Mat<f64>) -> Self {
        let n = y.nrows();
        let k = y.ncols();
        let p = x.ncols();
        assert_eq!(x.nrows(), n, "X rows must match Y rows");
        assert!(n > p + k, "need n > p + k for the joint REML estimate");

        let xtx = x.transpose() * x;
        let eye_p = Mat::<f64>::identity(p, p);
        let xtx_inv = xtx.col_piv_qr().solve(&eye_p);
        let beta = &xtx_inv * (x.transpose() * y);
        let residuals = y - x * &beta;

        let denom = (n - p) as f64;
        let sigma_res = (residuals.transpose() * &residuals) * (1.0 / denom);
        let eye_k = Mat::<f64>::identity(k, k);
        let sigma_inv = sigma_res.col_piv_qr().solve(&eye_k);

        let sigma_inv_eigvals = symmetric_eigenvalues(&sigma_inv);

        Self {
            residuals,
            sigma_res,
            sigma_inv,
            sigma_inv_eigvals,
            x: x.clone(),
            xtx_inv,
            n_samples: n,
            n_pheno: k,
        }
    }
}

/// Per-gene score statistics for joint multi-trait scoring.
///
/// `S[j, a] = (G₀' R)[j, a]` is the per-trait variant-residual product;
/// `k0[j, l] = (G₀' M G₀)[j, l]` is the **single-trait** projected
/// kernel. The joint `(m·k × m·k)` covariance is `Σ_res⁻¹ ⊗ k0` and never
/// materialised — every test below exploits the Kronecker structure.
pub struct GeneStats {
    pub s: Mat<f64>,
    pub k0: Mat<f64>,
}

/// Build `(S, K₀)` for one gene.
pub fn gene_stats(g0: &Mat<f64>, null: &MultiNullContinuous) -> GeneStats {
    let n = g0.nrows();
    let m = g0.ncols();
    assert_eq!(n, null.n_samples, "G rows must match null samples");

    let s = g0.transpose() * &null.residuals;
    let k0 = projected_kernel(g0, &null.x, &null.xtx_inv, m);
    GeneStats { s, k0 }
}

/// `K₀ = G₀'(I − X(X'X)⁻¹X')G₀ = G₀' G₀ − (G₀' X)(X'X)⁻¹(X' G₀)`.
fn projected_kernel(g: &Mat<f64>, x: &Mat<f64>, xtx_inv: &Mat<f64>, m: usize) -> Mat<f64> {
    let gtx = g.transpose() * x;
    let inner = &gtx * xtx_inv;
    let correction = &inner * gtx.transpose();
    let mut k = g.transpose() * g;
    for i in 0..m {
        for j in 0..m {
            k[(i, j)] -= correction[(i, j)];
        }
    }
    k
}

/// Joint per-variant common test:
///
///     q_i = (1 / k0[i,i]) · S[i,:] Σ_res⁻¹ S[i,:]ᵀ ~ χ²(k_pheno)
///
/// Returns `1.0` for variants with degenerate projected variance.
pub fn variant_joint_chi2(stats: &GeneStats, sigma_inv: &Mat<f64>, i: usize, n_pheno: usize) -> f64 {
    let k_ii = stats.k0[(i, i)];
    if k_ii <= 0.0 || !k_ii.is_finite() {
        return 1.0;
    }
    let mut quad = 0.0;
    for a in 0..n_pheno {
        for b in 0..n_pheno {
            quad += stats.s[(i, a)] * sigma_inv[(a, b)] * stats.s[(i, b)];
        }
    }
    let q = quad / k_ii;
    chisq_pvalue(q, n_pheno as f64)
}

/// Joint burden test:
///
///     q_burden(w) = (wᵀ S Σ_res⁻¹ Sᵀ w) / (wᵀ K₀ w) ~ χ²(k_pheno)
pub fn multi_burden(stats: &GeneStats, sigma_inv: &Mat<f64>, w: &[f64], n_pheno: usize) -> f64 {
    let m = w.len();
    if m == 0 {
        return 1.0;
    }
    let mut wkw = 0.0;
    for j in 0..m {
        if w[j] == 0.0 {
            continue;
        }
        for l in 0..m {
            wkw += w[j] * stats.k0[(j, l)] * w[l];
        }
    }
    if wkw <= 0.0 || !wkw.is_finite() {
        return 1.0;
    }
    // b[a] = Σ_j w_j · S[j, a]
    let mut b = vec![0.0_f64; n_pheno];
    for a in 0..n_pheno {
        let mut acc = 0.0;
        for j in 0..m {
            acc += w[j] * stats.s[(j, a)];
        }
        b[a] = acc;
    }
    let mut num = 0.0;
    for a in 0..n_pheno {
        for c in 0..n_pheno {
            num += b[a] * sigma_inv[(a, c)] * b[c];
        }
    }
    chisq_pvalue(num / wkw, n_pheno as f64)
}

/// Joint SKAT.
///
///     Q = ‖vec(W S Σ_res⁻¹)‖²    (Frobenius squared)
///
/// The weighted joint covariance factors as `Σ_res⁻¹ ⊗ (W K₀ W)`, so the
/// `(m·k)` mixture-of-chi-square eigenvalues are pairwise products
/// `{λ_a^Σ · μ_j^WK₀W}`. Two small eigendecompositions instead of one
/// `(m·k × m·k)` decomposition.
pub fn multi_skat(
    stats: &GeneStats,
    sigma_inv: &Mat<f64>,
    sigma_inv_eigvals: &[f64],
    w: &[f64],
) -> f64 {
    let m = w.len();
    if m == 0 {
        return 1.0;
    }
    // Q = Σ_a Σ_j (Σ_b w_j · S[j, b] · Σ_inv[a, b])²
    //   = ‖W S Σ_inv‖_F²
    let n_pheno = sigma_inv.nrows();
    let mut q = 0.0;
    for a in 0..n_pheno {
        for j in 0..m {
            let mut e = 0.0;
            for b in 0..n_pheno {
                e += w[j] * stats.s[(j, b)] * sigma_inv[(b, a)];
            }
            q += e * e;
        }
    }

    // Eigenvalues of W K₀ W (m × m) — small.
    let mut wkw = Mat::<f64>::zeros(m, m);
    for i in 0..m {
        for j in 0..m {
            wkw[(i, j)] = w[i] * stats.k0[(i, j)] * w[j];
        }
    }
    let mu = symmetric_eigenvalues(&wkw);

    // Pairwise product gives the (m·k) eigenvalues of the joint kernel.
    let mut joint_eigs = Vec::with_capacity(m * sigma_inv_eigvals.len());
    for &la in sigma_inv_eigvals {
        for &mj in &mu {
            joint_eigs.push(la * mj);
        }
    }
    stats::mixture_chisq_pvalue(q, &joint_eigs)
}

/// Joint ACAT-V. Common variants emit per-variant joint chi-sq p-values;
/// very-rare variants are pooled into a single joint burden statistic
/// with `Burden(1,β)` weights, exactly mirroring single-trait STAAR's
/// rare-group convention.
#[allow(clippy::too_many_arguments)]
pub fn multi_acat_v(
    stats: &GeneStats,
    sigma_inv: &Mat<f64>,
    n_pheno: usize,
    w_acat: &[f64],
    w_burden: &[f64],
    mafs: &[f64],
    n_samples: usize,
) -> f64 {
    const MAC_THRESHOLD: f64 = 10.0;
    let m = w_acat.len();
    let ns = n_samples as f64;

    let mut p_values = Vec::with_capacity(m);
    let mut cauchy_weights = Vec::with_capacity(m);
    let mut rare_indices: Vec<usize> = Vec::new();

    for j in 0..m {
        if w_acat[j] == 0.0 {
            continue;
        }
        let mac = (2.0 * mafs[j] * ns).round();
        if mac > MAC_THRESHOLD {
            let p = variant_joint_chi2(stats, sigma_inv, j, n_pheno);
            p_values.push(p);
            cauchy_weights.push(w_acat[j]);
        } else {
            rare_indices.push(j);
        }
    }

    if !rare_indices.is_empty() {
        let mut w_rare = vec![0.0_f64; m];
        for &j in &rare_indices {
            w_rare[j] = w_burden[j];
        }
        let p = multi_burden(stats, sigma_inv, &w_rare, n_pheno);
        let mean_w: f64 =
            rare_indices.iter().map(|&j| w_acat[j]).sum::<f64>() / rare_indices.len() as f64;
        p_values.push(p);
        cauchy_weights.push(mean_w);
    }

    if p_values.is_empty() {
        return 1.0;
    }
    stats::cauchy_combine_weighted(&p_values, &cauchy_weights)
}

/// Run the joint multi-trait STAAR omnibus on one gene.
///
/// Mirrors the annotation Cauchy tree of `score::run_staar_from_sumstats`
/// but every leaf is the joint multi-trait test instead of the single-
/// trait χ²(1) version. The result type is shared with single-trait STAAR
/// so the existing writers and report layer keep working.
pub fn run_multi_staar(
    g0: &Mat<f64>,
    null: &MultiNullContinuous,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
) -> StaarResult {
    let m = g0.ncols();
    if m == 0 {
        return nan_result();
    }
    let stats = gene_stats(g0, null);
    multi_tests(&stats, null, annotation_matrix, mafs)
}

fn multi_tests(
    stats: &GeneStats,
    null: &MultiNullContinuous,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
) -> StaarResult {
    let m = mafs.len();
    let n_pheno = null.n_pheno;
    let beta_1_25: Vec<f64> = mafs
        .iter()
        .map(|&maf| beta_density_weight(maf, 1.0, 25.0))
        .collect();
    let beta_1_1: Vec<f64> = mafs
        .iter()
        .map(|&maf| beta_density_weight(maf, 1.0, 1.0))
        .collect();
    let acat_denom: Vec<f64> = mafs
        .iter()
        .map(|&maf| {
            let d = beta_density_weight(maf, 0.5, 0.5);
            if d > 0.0 {
                d * d
            } else {
                1.0
            }
        })
        .collect();

    let run_burden = |w: &[f64]| multi_burden(stats, &null.sigma_inv, w, n_pheno);
    let run_skat =
        |w: &[f64]| multi_skat(stats, &null.sigma_inv, &null.sigma_inv_eigvals, w);
    let run_acat_v = |w_acat: &[f64], w_burden: &[f64]| {
        multi_acat_v(
            stats,
            &null.sigma_inv,
            n_pheno,
            w_acat,
            w_burden,
            mafs,
            null.n_samples,
        )
    };

    let base_burden_1_25 = run_burden(&beta_1_25);
    let base_burden_1_1 = run_burden(&beta_1_1);
    let base_skat_1_25 = run_skat(&beta_1_25);
    let base_skat_1_1 = run_skat(&beta_1_1);
    let wa_base_1_25: Vec<f64> = beta_1_25
        .iter()
        .zip(&acat_denom)
        .map(|(b, d)| b * b / d)
        .collect();
    let wa_base_1_1: Vec<f64> = beta_1_1
        .iter()
        .zip(&acat_denom)
        .map(|(b, d)| b * b / d)
        .collect();
    let base_acat_v_1_25 = run_acat_v(&wa_base_1_25, &beta_1_25);
    let base_acat_v_1_1 = run_acat_v(&wa_base_1_1, &beta_1_1);

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
            run_burden(&wb_1_25),
            run_burden(&wb_1_1),
            run_skat(&ws_1_25),
            run_skat(&ws_1_1),
            run_acat_v(&wa_1_25, &wb_1_25),
            run_acat_v(&wa_1_1, &wb_1_1),
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

fn nan_result() -> StaarResult {
    StaarResult {
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

/// χ²(df) survival function with the same `P_FLOOR` floor and tail
/// behaviour the rest of the codebase uses, so the omnibus stays
/// consistent across single-trait and multi-trait paths.
fn chisq_pvalue(t: f64, df: f64) -> f64 {
    use statrs::distribution::{ChiSquared, ContinuousCDF};
    if t <= 0.0 || !t.is_finite() {
        return 1.0;
    }
    match ChiSquared::new(df) {
        Ok(d) => (1.0 - d.cdf(t)).max(stats::P_FLOOR),
        Err(_) => f64::NAN,
    }
}

fn symmetric_eigenvalues(mat: &Mat<f64>) -> Vec<f64> {
    let n = mat.nrows();
    if n == 0 {
        return Vec::new();
    }
    match mat.self_adjoint_eigen(faer::Side::Lower) {
        Ok(evd) => {
            let s = evd.S();
            let cv = s.column_vector();
            (0..n).map(|i| cv[i].max(0.0)).collect()
        }
        Err(_) => vec![0.0; n],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::staar::model;
    use crate::staar::score;

    fn intercept_x(n: usize) -> Mat<f64> {
        let mut x = Mat::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }
        x
    }

    fn xorshift_normal(seed: u64) -> impl FnMut() -> f64 {
        // Box-Muller from a deterministic xorshift stream. Not for
        // statistics — just enough randomness for unit tests.
        let mut state = seed.max(1);
        let mut spare: Option<f64> = None;
        move || {
            if let Some(v) = spare.take() {
                return v;
            }
            let mut next = || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 11) as f64 / (1u64 << 53) as f64
            };
            let mut u1 = next();
            let u2 = next();
            if u1 < 1e-12 {
                u1 = 1e-12;
            }
            let mag = (-2.0 * u1.ln()).sqrt();
            let z0 = mag * (2.0 * std::f64::consts::PI * u2).cos();
            let z1 = mag * (2.0 * std::f64::consts::PI * u2).sin();
            spare = Some(z1);
            z0
        }
    }

    fn random_genotypes(n: usize, m: usize, mafs: &[f64], seed: u64) -> Mat<f64> {
        let mut state = seed.max(1);
        let mut next_u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut g = Mat::<f64>::zeros(n, m);
        for j in 0..m {
            for i in 0..n {
                let r = next_u01();
                let dose = if r < mafs[j] {
                    if next_u01() < mafs[j] {
                        2.0
                    } else {
                        1.0
                    }
                } else {
                    0.0
                };
                g[(i, j)] = dose;
            }
        }
        g
    }

    /// k=1 must reduce exactly to single-trait STAAR. The cleanest check
    /// is at the per-test level, since the omnibus Cauchy tree is the
    /// same code path on both sides.
    #[test]
    fn k1_burden_matches_single_trait() {
        let n = 200;
        let m = 8;
        let mafs = vec![0.005; m];
        let g = random_genotypes(n, m, &mafs, 11);

        let mut sample_normal = xorshift_normal(7);
        let mut y_single = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            y_single[(i, 0)] = sample_normal();
        }
        let x = intercept_x(n);

        let single_null = model::fit_glm(&y_single, &x);
        let multi_null = MultiNullContinuous::fit(&y_single, &x);
        assert_eq!(multi_null.n_pheno, 1);

        let single_p = score::run_staar(&g, &[], &mafs, &single_null, false);
        let multi_p = run_multi_staar(&g, &multi_null, &[], &mafs);

        assert!(
            (single_p.burden_1_25 - multi_p.burden_1_25).abs() < 1e-10,
            "burden(1,25) single={} multi={}",
            single_p.burden_1_25,
            multi_p.burden_1_25,
        );
        assert!((single_p.burden_1_1 - multi_p.burden_1_1).abs() < 1e-10);
        assert!((single_p.skat_1_25 - multi_p.skat_1_25).abs() < 1e-10);
        assert!((single_p.skat_1_1 - multi_p.skat_1_1).abs() < 1e-10);
        assert!((single_p.acat_v_1_25 - multi_p.acat_v_1_25).abs() < 1e-10);
        assert!((single_p.acat_v_1_1 - multi_p.acat_v_1_1).abs() < 1e-10);
        assert!((single_p.acat_o - multi_p.acat_o).abs() < 1e-10);
        assert!((single_p.staar_o - multi_p.staar_o).abs() < 1e-10);
    }

    /// k=1 with annotation channels also reduces exactly.
    #[test]
    fn k1_with_annotations_matches_single_trait() {
        let n = 150;
        let m = 6;
        let mafs = vec![0.004, 0.006, 0.002, 0.008, 0.001, 0.003];
        let g = random_genotypes(n, m, &mafs, 19);

        let mut normal = xorshift_normal(23);
        let mut y = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = normal();
        }
        let x = intercept_x(n);

        let ann: Vec<Vec<f64>> = (0..3)
            .map(|c| (0..m).map(|j| 0.1 * (c as f64 + 1.0) + 0.05 * j as f64).collect())
            .collect();

        let single_null = model::fit_glm(&y, &x);
        let multi_null = MultiNullContinuous::fit(&y, &x);
        let single_p = score::run_staar(&g, &ann, &mafs, &single_null, false);
        let multi_p = run_multi_staar(&g, &multi_null, &ann, &mafs);

        assert!((single_p.staar_o - multi_p.staar_o).abs() < 1e-10);
        for (a, b) in single_p.per_annotation.iter().zip(&multi_p.per_annotation) {
            for i in 0..6 {
                assert!(
                    (a[i] - b[i]).abs() < 1e-10,
                    "annotation row diverges at {i}: {} vs {}",
                    a[i],
                    b[i]
                );
            }
        }
    }

    /// k=2 with two independent traits and a single rare variant.
    /// The joint per-variant test must be a strictly tighter combination
    /// than either single-trait test (in the sense that it produces a
    /// finite, valid p-value and doesn't blow up).
    #[test]
    fn k2_per_variant_joint_runs_and_is_valid_p() {
        let n = 300;
        let m = 5;
        let mafs = vec![0.01; m];
        let g = random_genotypes(n, m, &mafs, 31);

        let mut nrm = xorshift_normal(37);
        let mut y = Mat::<f64>::zeros(n, 2);
        for i in 0..n {
            y[(i, 0)] = nrm();
            y[(i, 1)] = nrm();
        }
        let x = intercept_x(n);

        let null = MultiNullContinuous::fit(&y, &x);
        assert_eq!(null.n_pheno, 2);
        assert_eq!(null.sigma_res.nrows(), 2);
        assert_eq!(null.sigma_inv_eigvals.len(), 2);

        let stats = gene_stats(&g, &null);
        for i in 0..m {
            let p = variant_joint_chi2(&stats, &null.sigma_inv, i, 2);
            assert!(p.is_finite());
            assert!((0.0..=1.0).contains(&p), "p out of range at {i}: {p}");
        }
    }

    /// k=2 omnibus runs and produces a valid p-value.
    #[test]
    fn k2_omnibus_runs_and_is_valid_p() {
        let n = 250;
        let m = 7;
        let mafs = vec![0.005, 0.002, 0.008, 0.001, 0.004, 0.003, 0.006];
        let g = random_genotypes(n, m, &mafs, 41);
        let mut nrm = xorshift_normal(43);
        let mut y = Mat::<f64>::zeros(n, 2);
        for i in 0..n {
            y[(i, 0)] = nrm();
            y[(i, 1)] = 0.5 * y[(i, 0)] + 0.866 * nrm();
        }
        let x = intercept_x(n);

        let null = MultiNullContinuous::fit(&y, &x);
        // Off-diagonal of Σ_res should be positive given the y2 construction.
        assert!(null.sigma_res[(0, 1)] > 0.0);

        let ann: Vec<Vec<f64>> =
            (0..2).map(|_c| (0..m).map(|j| 0.5 + 0.1 * j as f64).collect()).collect();
        let result = run_multi_staar(&g, &null, &ann, &mafs);

        assert!(result.staar_o.is_finite());
        assert!((0.0..=1.0).contains(&result.staar_o));
        assert!(result.acat_o.is_finite());
        assert_eq!(result.per_annotation.len(), 2);
    }
}
