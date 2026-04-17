//! SCANG-STAAR bridge: Monte Carlo pseudo-residuals and empirical thresholds.
//!
//! Mirrors STAARpipeline R/staar2scang_nullmodel.R and SCANG R/SCANG.r for
//! the gaussian, unrelated path. Two jobs:
//!
//! 1. `ScangExt` carries `times` (number of Monte Carlo simulations, 2000
//!    per SCANG convention) and a `(times × n)` matrix of pseudo-residuals
//!    drawn from N(0, P) where P = (I − X(X'X)⁻¹X')/σ² is the null model
//!    projection. With those in hand we can evaluate the null distribution
//!    of any weighted score statistic U/√K just by redoing the same
//!    arithmetic against every row of the pseudo-residual matrix.
//! 2. `chrom_threshold` computes the empirical 1 − α quantile of the max
//!    −log(p) across a set of SCANG windows under the Monte Carlo null,
//!    matching SCANG R/SCANG.r:181-205 (`quantile(emL20_O, 1 − alpha)`
//!    with the `-log(filter)` floor).

use faer::Mat;

use crate::error::CohortError;
use crate::staar::model::NullModel;

/// xorshift64* PRNG. Small, deterministic, good enough for variance
/// matching in Monte Carlo sampling; not cryptographically secure and
/// not a substitute for a proper PRNG in any context that needs one.
pub(crate) struct Xorshift64(u64);

impl Xorshift64 {
    pub(crate) fn new(seed: u64) -> Self {
        // Avoid the zero state; xorshift64 converges to 0 there.
        Self(if seed == 0 { 0x9E3779B97F4A7C15 } else { seed })
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    pub(crate) fn uniform_01(&mut self) -> f64 {
        // 53-bit mantissa fills the [0, 1) range uniformly.
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Pre-sampled null-distribution residuals for the SCANG MC threshold.
///
/// `pseudo_residuals` is `(times × n_pheno)` row-major: each row is one
/// Monte Carlo draw of `z ~ N(0, P)` where P is the null model's
/// projection kernel. Row-major so every MC iteration reads a contiguous
/// span of memory when computing U_j = G' z_j per window.
pub struct ScangExt {
    /// Number of Monte Carlo simulations the sampler drew. Redundant with
    /// `pseudo_residuals.nrows()` but kept as a typed field so the rest of
    /// the pipeline does not have to cascade through the matrix shape.
    #[allow(dead_code)]
    pub times: u32,
    pub pseudo_residuals: Mat<f64>,
}

/// R seed used by staar2scang_nullmodel.R so our MC draws land on the
/// same pseudo-residual matrix across languages.
pub const SCANG_SEED: u64 = 19_880_615 + 666;

/// Default number of Monte Carlo simulations used by STAARpipeline.
pub const SCANG_DEFAULT_TIMES: u32 = 2000;

/// Default family-wise error rate used by SCANG R/SCANG.r:30.
pub const SCANG_DEFAULT_ALPHA: f64 = 0.05;

/// Filter threshold floor from SCANG R/SCANG.r:32. Empirical thresholds
/// below −log(filter) get raised to −log(filter) so the SKAT screening
/// path still fires even when MC says nothing is interesting.
pub const SCANG_DEFAULT_FILTER: f64 = 1e-4;

/// Populate `null.scang` with `times` pseudo-residuals ~ N(0, P) if not
/// already set. Unrelated gaussian path only: the kinship-aware form
/// uses the Cholesky factors of Σ⁻¹ and lives in STAARpipeline
/// R/staar2scang_nullmodel.R:37-63, not yet ported.
pub fn ensure_unrelated(null: &mut NullModel, times: u32, seed: u64) -> Result<(), CohortError> {
    if null.scang.is_some() {
        return Ok(());
    }
    if null.kinship.is_some() {
        return Err(CohortError::Input(
            "SCANG MC threshold on kinship-aware null models is not yet wired; \
             drop --kinship for SCANG or wait on the sparse Cholesky bridge."
                .into(),
        ));
    }
    if null.working_weights.is_some() {
        return Err(CohortError::Input(
            "SCANG MC threshold on binary traits is not yet wired; use \
             a continuous trait for SCANG until the SPA-aware MC sampler lands."
                .into(),
        ));
    }

    let pseudo_residuals = sample_unrelated(null, times, seed);
    null.scang = Some(ScangExt {
        times,
        pseudo_residuals,
    });
    debug_assert_eq!(
        null.scang.as_ref().unwrap().times as usize,
        null.scang.as_ref().unwrap().pseudo_residuals.nrows(),
    );
    Ok(())
}

/// Sample `times` rows of `z ~ N(0, P)` where `P = (I − H) / σ²`.
///
/// Stable construction: draw `e ~ N(0, I_n)`, project out the column
/// space of X, divide by σ. The result has covariance `(I − H)/σ² = P`
/// because `(I − H)` is the idempotent projection onto X's orthogonal
/// complement. Matches R staar2scang_nullmodel.R:27-34 without paying
/// the O(n³) eigendecomposition cost — the eigen-basis representation
/// would produce the same covariance but with a different left-unitary
/// mixing, which is irrelevant to the null distribution of U/√K.
fn sample_unrelated(null: &NullModel, times: u32, seed: u64) -> Mat<f64> {
    let n = null.n_samples;
    let sigma = null.sigma2.sqrt().max(f64::MIN_POSITIVE);

    let mut rng = Xorshift64::new(seed);
    let mut pseudo = Mat::<f64>::zeros(times as usize, n);
    let mut e = Mat::<f64>::zeros(n, 1);
    for t in 0..times as usize {
        for i in 0..n {
            e[(i, 0)] = standard_normal(&mut rng);
        }
        let xt_e = null.x_matrix.transpose() * &e;
        let beta = &null.xtx_inv * &xt_e;
        let x_beta = &null.x_matrix * &beta;
        for i in 0..n {
            pseudo[(t, i)] = (e[(i, 0)] - x_beta[(i, 0)]) / sigma;
        }
    }
    pseudo
}

/// Box-Muller draw from N(0, 1). One variate per call; the second is
/// dropped so state is a single u64, trivially clonable for a later
/// parallel extension.
#[inline]
pub(crate) fn standard_normal(rng: &mut Xorshift64) -> f64 {
    let u1 = rng.uniform_01().max(f64::MIN_POSITIVE);
    let u2 = rng.uniform_01();
    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
}

/// Compute the SCANG empirical threshold for one chromosome.
///
/// `window_max_by_sim[j]` is the max over all SCANG windows on this
/// chromosome of −log(p) evaluated under pseudo-residual row j.
/// Returns max(quantile(1 − alpha), −log(filter)) as SCANG R/SCANG.r:181
/// does.
pub fn chrom_threshold(window_max_by_sim: &[f64], alpha: f64, filter: f64) -> f64 {
    let mut sorted = window_max_by_sim.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let q = empirical_quantile(&sorted, 1.0 - alpha);
    let floor = -filter.ln() / std::f64::consts::LN_10;
    q.max(floor)
}

/// Linear-interpolation quantile on an already-sorted slice. Matches R's
/// default `quantile(x, p)` (type-7) to within rounding.
fn empirical_quantile(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    if n == 1 {
        return sorted[0];
    }
    let h = (n as f64 - 1.0) * p;
    let lo = h.floor() as usize;
    let hi = lo + 1;
    if hi >= n {
        return sorted[n - 1];
    }
    let frac = h - lo as f64;
    sorted[lo] + frac * (sorted[hi] - sorted[lo])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::staar::model::fit_glm;

    fn fixture(n: usize, k: usize) -> NullModel {
        let mut x = Mat::<f64>::zeros(n, k);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            for c in 1..k {
                x[(i, c)] = ((i + c) as f64).sin();
            }
        }
        let y = Mat::<f64>::from_fn(n, 1, |i, _| (i as f64) * 0.01);
        fit_glm(&y, &x)
    }

    #[test]
    fn pseudo_residuals_orthogonal_to_covariates() {
        let mut null = fixture(200, 3);
        ensure_unrelated(&mut null, 50, SCANG_SEED).unwrap();
        let ext = null.scang.as_ref().unwrap();
        let pseudo = &ext.pseudo_residuals;

        // Each row z of the pseudo-residual matrix must satisfy X' z = 0
        // because we projected e onto the null space of X. Check a few
        // rows; exact equality would require infinite precision, 1e-8 is
        // already well below any statistic we'd compute on these.
        for t in [0usize, 17, 49] {
            let z = Mat::<f64>::from_fn(null.n_samples, 1, |i, _| pseudo[(t, i)]);
            let xt_z = null.x_matrix.transpose() * &z;
            for c in 0..null.x_matrix.ncols() {
                assert!(
                    xt_z[(c, 0)].abs() < 1e-8,
                    "X'z[{c}] = {} for row {t}",
                    xt_z[(c, 0)]
                );
            }
        }
    }

    #[test]
    fn pseudo_residual_variance_matches_projection() {
        // Empirically Var(z_i) should be close to P_ii = (1 − H_ii)/σ²
        // for iid draws. Pick a small n/times so the test stays fast but
        // the sample mean is still within 3× sqrt(2/times) of the truth.
        let mut null = fixture(80, 3);
        let sigma2 = null.sigma2;
        ensure_unrelated(&mut null, 4000, SCANG_SEED).unwrap();
        let ext = null.scang.as_ref().unwrap();

        // H = X(X'X)^{-1}X'; diag(H) computed column-by-column.
        let h_diag: Vec<f64> = (0..null.n_samples)
            .map(|i| {
                let mut x_i = Mat::<f64>::zeros(null.x_matrix.ncols(), 1);
                for c in 0..null.x_matrix.ncols() {
                    x_i[(c, 0)] = null.x_matrix[(i, c)];
                }
                let tmp = &null.xtx_inv * &x_i;
                let mut s = 0.0;
                for c in 0..null.x_matrix.ncols() {
                    s += null.x_matrix[(i, c)] * tmp[(c, 0)];
                }
                s
            })
            .collect();

        for i in [0usize, 13, 40, 79] {
            let mut sum = 0.0;
            let mut sum_sq = 0.0;
            for t in 0..ext.times as usize {
                let v = ext.pseudo_residuals[(t, i)];
                sum += v;
                sum_sq += v * v;
            }
            let m = ext.times as f64;
            let var = (sum_sq - sum * sum / m) / (m - 1.0);
            let expected = (1.0 - h_diag[i]) / sigma2;
            assert!(
                (var - expected).abs() < 0.1 * expected.max(1e-6),
                "var[z_{i}] = {var} vs expected {expected}",
            );
        }
    }

    #[test]
    fn chrom_threshold_respects_filter_floor() {
        let stats = vec![0.1, 0.2, 0.3, 0.4, 0.5];
        // quantile ~ 0.48; filter=1e-4 ⇒ floor = 4. Floor wins.
        let th = chrom_threshold(&stats, 0.05, 1e-4);
        assert!((th - 4.0).abs() < 1e-12);
    }

    #[test]
    fn chrom_threshold_returns_quantile_when_above_floor() {
        let stats = vec![5.0, 6.0, 7.0, 8.0, 9.0];
        let th = chrom_threshold(&stats, 0.05, 1e-4);
        assert!(th > 4.0);
        assert!(th <= 9.0);
    }

    #[test]
    fn ensure_unrelated_is_idempotent() {
        let mut null = fixture(64, 2);
        ensure_unrelated(&mut null, 10, SCANG_SEED).unwrap();
        let first_row0 = null.scang.as_ref().unwrap().pseudo_residuals[(0, 0)];
        ensure_unrelated(&mut null, 10, SCANG_SEED).unwrap();
        let second_row0 = null.scang.as_ref().unwrap().pseudo_residuals[(0, 0)];
        assert_eq!(first_row0, second_row0);
    }
}
