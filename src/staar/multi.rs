// `needless_range_loop` is silenced because the kernels are matrix
// arithmetic where index loops read more naturally than enumerate
// chains over `Mat` rows.
#![allow(clippy::needless_range_loop)]

//! Joint multi-trait STAAR (Li et al. 2023, MultiSTAAR) for unrelated
//! traits with shared covariates.
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
//!
//! **Binary extension.** GMMAT `glmmkin.R:43` rejects non-gaussian families
//! in the multi-pheno path, so MultiSTAAR has no upstream joint-binary null
//! fitter. Our `fit_binary` is a lift that combines two upstreams:
//! per-trait logistic IRLS (GMMAT single-trait binomial, mirrored in
//! `model::fit_logistic`) with MultiSTAAR's gaussian Kronecker pattern.
//! Concretely we fit each of k traits independently via single-trait
//! logistic IRLS, form R = Y − μ and take `Σ_res` as the cross-trait
//! residual **correlation** with unit diagonal (binomial dispersion is 1,
//! so using the raw covariance would shift k=1 off single-trait parity),
//! and use a shared working weight `W̄ = (1/k) Σ_j W_j` to build the
//! weighted projection `K₀ = G₀'(W̄ − W̄X(X'W̄X)⁻¹X'W̄)G₀`, i.e. the
//! single-trait Fisher-information kernel at the mean working weight. At
//! k=1 this reduces exactly to single-trait logistic: `Σ_res = [1]`,
//! `W̄ = W_1`, `K₀ = K_single_trait`. For k > 1 with unequal per-trait
//! weights this is a **shared-W approximation**: the joint covariance is
//! approximated by `Σ_res_corr ⊗ K̄` instead of the block-varying
//! `blockdiag(K_j)` a per-trait-weighted formulation would produce. This
//! preserves the Kronecker hot path and is the pragmatic first cut in the
//! absence of an upstream reference.

use faer::prelude::*;
use faer::Mat;

use super::score::{beta_density_weight, StaarResult};
use super::stats;

/// Null model fit for unrelated multi-trait STAAR, continuous or binary.
///
/// The gaussian path stores OLS residuals against the shared covariate
/// matrix plus `(X'X)⁻¹`. The binary path stores per-trait IRLS residuals
/// `Y − μ`, the mean working-weight vector, and `(X'W̄X)⁻¹`. Which path
/// produced the fit is recorded in `working_weights`: `None` for gaussian,
/// `Some(W̄)` for binary. Every hot-path function branches on this one
/// field; all other fields keep the same Kronecker semantics.
pub struct MultiNull {
    /// Gaussian: `R = M Y`, shape `(n, k)`. Binary: `R = Y − μ` per trait,
    /// raw score residuals (not working responses).
    pub residuals: Mat<f64>,
    /// Gaussian: `Σ_res = R'R / (n − p)`, the unbiased residual covariance.
    /// Binary: `Σ_res` is the cross-trait residual **correlation** with
    /// unit diagonal. The unit diagonal keeps k=1 binary aligned with
    /// single-trait logistic (binomial dispersion is 1), and preserves the
    /// Kronecker structure `Σ_res ⊗ K̄` for the joint tests.
    pub sigma_res: Mat<f64>,
    /// `Σ_res⁻¹`, pre-computed for the hot path.
    pub sigma_inv: Mat<f64>,
    /// Eigenvalues of `Σ_res⁻¹`. Cached so the joint SKAT eigenvalue
    /// product is cheap.
    pub sigma_inv_eigvals: Vec<f64>,
    /// Shared covariate matrix, `(n, p)`.
    pub x: Mat<f64>,
    /// Gaussian: `(X'X)⁻¹`. Binary: `(X'W̄X)⁻¹` with `W̄ = (1/k) Σ_j W_j`.
    /// Used by `projected_kernel` to build the single-trait projected
    /// kernel `K₀` from any gene's `G₀`.
    pub xtx_inv: Mat<f64>,
    /// Per-sample mean working weight `W̄_i = (1/k) Σ_j μ_{ij}(1 − μ_{ij})`
    /// for binary traits, populated by `fit_binary`. `None` on the
    /// gaussian path. When `Some`, `projected_kernel` uses the
    /// Fisher-information form `K₀ = G₀'(W̄ − W̄X(X'W̄X)⁻¹X'W̄)G₀`; when
    /// `None`, it uses the OLS hat-removing form `K₀ = G₀'(I − H)G₀`.
    pub working_weights: Option<Vec<f64>>,
    /// Converged multi-trait AI-REML state when kinship is in use
    /// (`--kinship` on the CLI). `None` on both the unrelated gaussian
    /// and unrelated binary paths. When `Some`, `run_multi_staar`
    /// dispatches to the kinship-aware score kernel in
    /// `crate::staar::multi_kinship`, which contracts against the full
    /// `(mk × mk)` joint covariance.
    pub kinship: Option<crate::staar::multi_kinship::MultiKinshipState>,
    pub n_samples: usize,
    pub n_pheno: usize,
}

impl MultiNull {
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
            working_weights: None,
            kinship: None,
            n_samples: n,
            n_pheno: k,
        }
    }

    /// Fit k binary traits jointly via per-trait logistic IRLS + empirical
    /// cross-trait residual correlation. See the module-level docstring
    /// for the derivation. At k=1 this is bit-equivalent to
    /// `model::fit_logistic`: `Σ_res = [1]`, `W̄ = W_1`, `xtx_inv = (X'W_1 X)⁻¹`.
    pub fn fit_binary(y: &Mat<f64>, x: &Mat<f64>, max_iter: usize) -> Self {
        let n = y.nrows();
        let k = y.ncols();
        let p = x.ncols();
        assert_eq!(x.nrows(), n, "X rows must match Y rows");
        assert!(n > p + k, "need n > p + k samples");

        let mut residuals = Mat::<f64>::zeros(n, k);
        let mut w_bar = vec![0.0_f64; n];
        let mut y_col = Mat::<f64>::zeros(n, 1);
        for j in 0..k {
            for i in 0..n {
                y_col[(i, 0)] = y[(i, j)];
            }
            let fit = super::model::fit_logistic(&y_col, x, max_iter);
            for i in 0..n {
                residuals[(i, j)] = fit.residuals[(i, 0)];
            }
            let w_j = fit
                .working_weights
                .as_ref()
                .expect("fit_logistic sets working_weights");
            for i in 0..n {
                w_bar[i] += w_j[i];
            }
        }
        let inv_k = 1.0 / k as f64;
        for i in 0..n {
            w_bar[i] *= inv_k;
        }

        // Cross-trait residual correlation with unit diagonal. Using the
        // empirical covariance R'R/(n-p) directly would give `Σ_res[j,j] =
        // mean(μ_j(1-μ_j))` for k=1 binary, breaking parity with
        // single-trait logistic (which assumes dispersion = 1). The unit
        // diagonal captures the fact that binomial dispersion is fixed; the
        // off-diagonal correlation captures the cross-trait coupling.
        let cov = (residuals.transpose() * &residuals) * (1.0 / (n - p) as f64);
        let mut sigma_res = Mat::<f64>::zeros(k, k);
        let mut d_inv = vec![0.0_f64; k];
        for j in 0..k {
            let v = cov[(j, j)].max(1e-30);
            d_inv[j] = 1.0 / v.sqrt();
        }
        for a in 0..k {
            for b in 0..k {
                sigma_res[(a, b)] = cov[(a, b)] * d_inv[a] * d_inv[b];
            }
        }
        // Numerical guard: force exact unit diagonal so sigma_inv at k=1
        // is exactly 1.0.
        for a in 0..k {
            sigma_res[(a, a)] = 1.0;
        }

        let eye_k = Mat::<f64>::identity(k, k);
        let sigma_inv = sigma_res.col_piv_qr().solve(&eye_k);
        let sigma_inv_eigvals = symmetric_eigenvalues(&sigma_inv);

        // `(X'W̄X)⁻¹` — the weighted normal-equation inverse at the mean
        // working weight. Identical to single-trait `fit_logistic` when k=1.
        let mut xtwx: Mat<f64> = Mat::zeros(p, p);
        for i in 0..n {
            let wi = w_bar[i].max(1e-30);
            for a in 0..p {
                for b in 0..p {
                    xtwx[(a, b)] += x[(i, a)] * wi * x[(i, b)];
                }
            }
        }
        let eye_p = Mat::<f64>::identity(p, p);
        let xtx_inv = xtwx.col_piv_qr().solve(&eye_p);

        Self {
            residuals,
            sigma_res,
            sigma_inv,
            sigma_inv_eigvals,
            x: x.clone(),
            xtx_inv,
            working_weights: Some(w_bar),
            kinship: None,
            n_samples: n,
            n_pheno: k,
        }
    }

    /// Fit the joint multi-trait GLMM with shared kinship matrices.
    /// Wraps `multi_kinship::fit_multi_kinship` and packages the state
    /// into a `MultiNull` so the pipeline's single dispatch point can
    /// consume any of the three variants uniformly. The no-kinship
    /// fields on `MultiNull` are populated with placeholder values in
    /// this variant (the score path routes through `kinship` and never
    /// reads them).
    pub fn fit_with_kinship(
        y: &Mat<f64>,
        x: &Mat<f64>,
        kinships: &[crate::staar::kinship::types::KinshipMatrix],
        groups: &crate::staar::kinship::types::GroupPartition,
    ) -> Result<Self, crate::error::CohortError> {
        let state = crate::staar::multi_kinship::fit_multi_kinship(y, x, kinships, groups)?;
        let n = y.nrows();
        let k = y.ncols();
        Ok(Self {
            // Placeholders for the kinship path. `run_multi_staar`
            // branches on `kinship.is_some()` and never reads these
            // fields; they exist to keep the struct layout uniform so
            // callers that only want `n_samples`/`n_pheno` still work.
            residuals: Mat::<f64>::zeros(n, k),
            sigma_res: Mat::<f64>::identity(k, k),
            sigma_inv: Mat::<f64>::identity(k, k),
            sigma_inv_eigvals: vec![1.0; k],
            x: x.clone(),
            xtx_inv: Mat::<f64>::zeros(x.ncols(), x.ncols()),
            working_weights: None,
            kinship: Some(state),
            n_samples: n,
            n_pheno: k,
        })
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

/// Per-gene scratch buffers threaded through the joint test kernels.
///
/// Allocated once at the top of `multi_tests` and reused across every
/// per-channel call to `multi_burden` / `multi_skat` / `multi_acat_v`.
/// Mirrors the `kernel_buf` pattern in `score::staar_tests`: allocate
/// once per gene, not once per test-per-channel.
pub struct MultiScratch {
    /// `b[a] = Σ_j w_j · S[j, a]`, length `n_pheno`.
    b: Vec<f64>,
    /// Weighted kernel `W K₀ W`, shape `m × m`, overwritten per call.
    wkw: Mat<f64>,
    /// Eigenvalues of `W K₀ W`, length `m`.
    mu_eigs: Vec<f64>,
    /// Kronecker pairwise products `{λ_a · μ_j}`, length `m · n_pheno`.
    joint_eigs: Vec<f64>,
    /// Burden-weight vector staged by the ACAT-V rare-group path,
    /// length `m`, cleared per call before repopulation.
    w_rare: Vec<f64>,
}

impl MultiScratch {
    pub fn with_capacity(m: usize, n_pheno: usize) -> Self {
        Self {
            b: vec![0.0; n_pheno],
            wkw: Mat::<f64>::zeros(m, m),
            mu_eigs: Vec::with_capacity(m),
            joint_eigs: Vec::with_capacity(m * n_pheno),
            w_rare: vec![0.0; m],
        }
    }
}

/// Build `(S, K₀)` for one gene.
pub fn gene_stats(g0: &Mat<f64>, null: &MultiNull) -> GeneStats {
    let n = g0.nrows();
    let m = g0.ncols();
    assert_eq!(n, null.n_samples, "G rows must match null samples");

    let s = g0.transpose() * &null.residuals;
    let k0 = match null.working_weights.as_ref() {
        None => projected_kernel_ols(g0, &null.x, &null.xtx_inv, m),
        Some(w) => projected_kernel_weighted(g0, &null.x, &null.xtx_inv, w, m),
    };
    GeneStats { s, k0 }
}

/// `K₀ = G₀'(I − X(X'X)⁻¹X')G₀ = G₀' G₀ − (G₀' X)(X'X)⁻¹(X' G₀)`.
fn projected_kernel_ols(g: &Mat<f64>, x: &Mat<f64>, xtx_inv: &Mat<f64>, m: usize) -> Mat<f64> {
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

/// `K₀ = G₀'(W̄ − W̄X(X'W̄X)⁻¹X'W̄)G₀` for the binary path, using a shared
/// per-sample mean working weight `W̄`. Matches single-trait
/// `NullModel::compute_kernel` at k=1.
fn projected_kernel_weighted(
    g: &Mat<f64>,
    x: &Mat<f64>,
    xtwx_inv: &Mat<f64>,
    w: &[f64],
    m: usize,
) -> Mat<f64> {
    let n = g.nrows();
    debug_assert_eq!(w.len(), n);
    let mut wg = Mat::<f64>::zeros(n, m);
    for j in 0..m {
        for i in 0..n {
            wg[(i, j)] = w[i] * g[(i, j)];
        }
    }
    let gtwg = g.transpose() * &wg;
    let xtwg = x.transpose() * &wg;
    let inner = xtwx_inv * &xtwg;
    let correction = wg.transpose() * (x * &inner);
    let mut k = gtwg;
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
///
/// `b_buf` is the caller-owned `n_pheno`-sized scratch used to stage
/// `b[a] = Σ_j w_j · S[j, a]` without allocating per call.
pub fn multi_burden(
    stats: &GeneStats,
    sigma_inv: &Mat<f64>,
    w: &[f64],
    n_pheno: usize,
    b_buf: &mut [f64],
) -> f64 {
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
    debug_assert_eq!(b_buf.len(), n_pheno);
    // b[a] = Σ_j w_j · S[j, a]
    for a in 0..n_pheno {
        let mut acc = 0.0;
        for j in 0..m {
            acc += w[j] * stats.s[(j, a)];
        }
        b_buf[a] = acc;
    }
    let mut num = 0.0;
    for a in 0..n_pheno {
        for c in 0..n_pheno {
            num += b_buf[a] * sigma_inv[(a, c)] * b_buf[c];
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
    wkw_buf: &mut Mat<f64>,
    mu_buf: &mut Vec<f64>,
    joint_eigs_buf: &mut Vec<f64>,
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

    // Eigenvalues of W K₀ W (m × m) — small. Reuse the scratch buffer;
    // every cell gets overwritten so there is no need to zero first.
    debug_assert_eq!(wkw_buf.nrows(), m);
    debug_assert_eq!(wkw_buf.ncols(), m);
    for i in 0..m {
        for j in 0..m {
            wkw_buf[(i, j)] = w[i] * stats.k0[(i, j)] * w[j];
        }
    }
    symmetric_eigenvalues_into(wkw_buf, mu_buf);

    // Pairwise product gives the (m·k) eigenvalues of the joint kernel.
    joint_eigs_buf.clear();
    for &la in sigma_inv_eigvals {
        for &mj in mu_buf.iter() {
            joint_eigs_buf.push(la * mj);
        }
    }
    stats::mixture_chisq_pvalue(q, joint_eigs_buf)
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
    w_rare_buf: &mut [f64],
    b_buf: &mut [f64],
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
        // Reuse the caller-owned scratch; a previous call may have left
        // stale values, so zero before repopulating.
        debug_assert_eq!(w_rare_buf.len(), m);
        for v in w_rare_buf.iter_mut() {
            *v = 0.0;
        }
        for &j in &rare_indices {
            w_rare_buf[j] = w_burden[j];
        }
        let p = multi_burden(stats, sigma_inv, w_rare_buf, n_pheno, b_buf);
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
    null: &MultiNull,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
) -> StaarResult {
    let m = g0.ncols();
    if m == 0 {
        return nan_result();
    }
    if let Some(ks) = null.kinship.as_ref() {
        return crate::staar::multi_kinship::run_multi_staar_kinship(
            g0,
            ks,
            annotation_matrix,
            mafs,
        );
    }
    let stats = gene_stats(g0, null);
    let mut scratch = MultiScratch::with_capacity(m, null.n_pheno);
    multi_tests(&stats, null, annotation_matrix, mafs, &mut scratch)
}

fn multi_tests(
    stats: &GeneStats,
    null: &MultiNull,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
    scratch: &mut MultiScratch,
) -> StaarResult {
    let m = mafs.len();
    let n_pheno = null.n_pheno;
    let sigma_inv = &null.sigma_inv;
    let sigma_inv_eigvals = &null.sigma_inv_eigvals;
    let n_samples = null.n_samples;

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

    // Each test call borrows only the scratch fields it actually needs —
    // disjoint-field borrows let us re-enter `scratch` on the next call
    // without aliasing.
    let base_burden_1_25 =
        multi_burden(stats, sigma_inv, &beta_1_25, n_pheno, &mut scratch.b);
    let base_burden_1_1 =
        multi_burden(stats, sigma_inv, &beta_1_1, n_pheno, &mut scratch.b);
    let base_skat_1_25 = multi_skat(
        stats,
        sigma_inv,
        sigma_inv_eigvals,
        &beta_1_25,
        &mut scratch.wkw,
        &mut scratch.mu_eigs,
        &mut scratch.joint_eigs,
    );
    let base_skat_1_1 = multi_skat(
        stats,
        sigma_inv,
        sigma_inv_eigvals,
        &beta_1_1,
        &mut scratch.wkw,
        &mut scratch.mu_eigs,
        &mut scratch.joint_eigs,
    );
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
    let base_acat_v_1_25 = multi_acat_v(
        stats,
        sigma_inv,
        n_pheno,
        &wa_base_1_25,
        &beta_1_25,
        mafs,
        n_samples,
        &mut scratch.w_rare,
        &mut scratch.b,
    );
    let base_acat_v_1_1 = multi_acat_v(
        stats,
        sigma_inv,
        n_pheno,
        &wa_base_1_1,
        &beta_1_1,
        mafs,
        n_samples,
        &mut scratch.w_rare,
        &mut scratch.b,
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
            multi_burden(stats, sigma_inv, &wb_1_25, n_pheno, &mut scratch.b),
            multi_burden(stats, sigma_inv, &wb_1_1, n_pheno, &mut scratch.b),
            multi_skat(
                stats,
                sigma_inv,
                sigma_inv_eigvals,
                &ws_1_25,
                &mut scratch.wkw,
                &mut scratch.mu_eigs,
                &mut scratch.joint_eigs,
            ),
            multi_skat(
                stats,
                sigma_inv,
                sigma_inv_eigvals,
                &ws_1_1,
                &mut scratch.wkw,
                &mut scratch.mu_eigs,
                &mut scratch.joint_eigs,
            ),
            multi_acat_v(
                stats,
                sigma_inv,
                n_pheno,
                &wa_1_25,
                &wb_1_25,
                mafs,
                n_samples,
                &mut scratch.w_rare,
                &mut scratch.b,
            ),
            multi_acat_v(
                stats,
                sigma_inv,
                n_pheno,
                &wa_1_1,
                &wb_1_1,
                mafs,
                n_samples,
                &mut scratch.w_rare,
                &mut scratch.b,
            ),
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
    let mut out = Vec::new();
    symmetric_eigenvalues_into(mat, &mut out);
    out
}

/// Overwrite `out` with the clamped eigenvalues of `mat`. Reuses the
/// caller-owned buffer so the hot path in `multi_skat` avoids a fresh
/// `Vec<f64>` per call.
fn symmetric_eigenvalues_into(mat: &Mat<f64>, out: &mut Vec<f64>) {
    out.clear();
    let n = mat.nrows();
    if n == 0 {
        return;
    }
    match mat.self_adjoint_eigen(faer::Side::Lower) {
        Ok(evd) => {
            let s = evd.S();
            let cv = s.column_vector();
            for i in 0..n {
                out.push(cv[i].max(0.0));
            }
        }
        Err(_) => {
            out.resize(n, 0.0);
        }
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
        let multi_null = MultiNull::fit(&y_single, &x);
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
        let multi_null = MultiNull::fit(&y, &x);
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

        let null = MultiNull::fit(&y, &x);
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

        let null = MultiNull::fit(&y, &x);
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

    /// Deterministic Bernoulli draws from xorshift so the binary tests
    /// stay reproducible without pulling in a heavy RNG dep.
    fn random_binary(n: usize, p: f64, seed: u64) -> Mat<f64> {
        let mut state = seed.max(1);
        let mut next_u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut y = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = if next_u01() < p { 1.0 } else { 0.0 };
        }
        y
    }

    /// k=1 binary must reduce to single-trait logistic. `MultiNull::fit_binary`
    /// with k=1 shares the same IRLS, residuals, and `(X'WX)⁻¹`, so the
    /// per-test p-values should match to numerical tolerance (floating
    /// point path differs slightly between direct and Kronecker-routed
    /// kernels, but well below any analytically meaningful threshold).
    #[test]
    fn k1_binary_matches_single_trait_logistic() {
        let n = 300;
        let m = 8;
        let mafs = vec![0.01; m];
        let g = random_genotypes(n, m, &mafs, 53);

        let y = random_binary(n, 0.3, 59);
        let x = intercept_x(n);

        let single_null = model::fit_logistic(&y, &x, 25);
        let multi_null = MultiNull::fit_binary(&y, &x, 25);
        assert_eq!(multi_null.n_pheno, 1);
        assert!(multi_null.working_weights.is_some());
        assert!((multi_null.sigma_res[(0, 0)] - 1.0).abs() < 1e-12);

        let single_p = score::run_staar(&g, &[], &mafs, &single_null, false);
        let multi_p = run_multi_staar(&g, &multi_null, &[], &mafs);

        let tol = 1e-8;
        assert!(
            (single_p.burden_1_25 - multi_p.burden_1_25).abs() < tol,
            "burden(1,25) single={} multi={}",
            single_p.burden_1_25, multi_p.burden_1_25,
        );
        assert!((single_p.burden_1_1 - multi_p.burden_1_1).abs() < tol);
        assert!((single_p.skat_1_25 - multi_p.skat_1_25).abs() < tol);
        assert!((single_p.skat_1_1 - multi_p.skat_1_1).abs() < tol);
        assert!((single_p.acat_v_1_25 - multi_p.acat_v_1_25).abs() < tol);
        assert!((single_p.acat_v_1_1 - multi_p.acat_v_1_1).abs() < tol);
        assert!((single_p.acat_o - multi_p.acat_o).abs() < tol);
        assert!((single_p.staar_o - multi_p.staar_o).abs() < tol);
    }

    /// k=2 binary omnibus runs and emits a valid p-value. Exercises the
    /// shared-W̄ projected kernel and the residual-correlation coupling.
    #[test]
    fn k2_binary_omnibus_runs_and_is_valid_p() {
        let n = 400;
        let m = 7;
        let mafs = vec![0.005, 0.008, 0.003, 0.006, 0.002, 0.004, 0.007];
        let g = random_genotypes(n, m, &mafs, 67);

        // Correlated binaries: y2 = y1 xor coin(0.3) to induce dependence.
        let y1 = random_binary(n, 0.35, 71);
        let mut state: u64 = 73;
        let mut next_u01 = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut y = Mat::<f64>::zeros(n, 2);
        for i in 0..n {
            y[(i, 0)] = y1[(i, 0)];
            let flip = if next_u01() < 0.3 { 1.0 } else { 0.0 };
            let xor = if (y1[(i, 0)] == 1.0) ^ (flip == 1.0) {
                1.0
            } else {
                0.0
            };
            y[(i, 1)] = xor;
        }
        let x = intercept_x(n);

        let null = MultiNull::fit_binary(&y, &x, 25);
        assert_eq!(null.n_pheno, 2);
        assert!(null.working_weights.is_some());
        // Unit diagonal by construction; off-diagonal is the cross-trait
        // residual correlation. We only assert finiteness here since the
        // exact value is a function of the xorshift stream.
        assert!((null.sigma_res[(0, 0)] - 1.0).abs() < 1e-12);
        assert!((null.sigma_res[(1, 1)] - 1.0).abs() < 1e-12);
        assert!(null.sigma_res[(0, 1)].is_finite());

        let ann: Vec<Vec<f64>> =
            (0..2).map(|_c| (0..m).map(|j| 0.5 + 0.1 * j as f64).collect()).collect();
        let result = run_multi_staar(&g, &null, &ann, &mafs);

        assert!(result.staar_o.is_finite());
        assert!((0.0..=1.0).contains(&result.staar_o));
        assert!(result.acat_o.is_finite());
        assert_eq!(result.per_annotation.len(), 2);
    }
}
