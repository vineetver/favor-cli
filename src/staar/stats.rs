use std::f64::consts::PI;

use statrs::distribution::{ChiSquared, Continuous, ContinuousCDF, Normal};

/// Cauchy combination test (CCT) — equal weights.
///
/// T = (1/K) * sum(tan((0.5 - p_i) * pi))
/// p = 0.5 - arctan(T) / pi
///
/// Handles edge cases: p near 0 or 1 via clamping. Skips NaN.
pub fn cauchy_combine(p_values: &[f64]) -> f64 {
    cauchy_combine_weighted(p_values, &[])
}

/// Smallest positive f64 (denormal). Floor for combined p-values so the
/// omnibus never underflows to exactly 0 (which becomes -log10(p) = +inf
/// in summary reports). Shared with the per-variant chi-squared survival
/// in `score.rs` so a single underflow doesn't poison Cauchy combination.
pub(super) const P_FLOOR: f64 = f64::MIN_POSITIVE * f64::EPSILON; // 5e-324

/// Threshold above which `atan(T)/π` saturates at 0.5 in f64. Past this we
/// switch to the upper-tail asymptotic `1/(π·T)` (STAAR/R/CCT.R).
const T_ASYMP: f64 = 1.0e15;

/// Weighted Cauchy combination.
///
/// T = sum(w_i * tan((0.5 - p_i) * pi)) / sum(w_i)
///
/// If weights is empty, uses equal weights. Weights need not sum to 1.
/// Used by ACAT-V where annotation weights modulate each variant's contribution.
pub fn cauchy_combine_weighted(p_values: &[f64], weights: &[f64]) -> f64 {
    let use_weights = !weights.is_empty();
    if use_weights {
        assert_eq!(p_values.len(), weights.len());
    }

    let mut t_sum = 0.0;
    let mut w_sum = 0.0;
    let mut count = 0;

    for (i, &p) in p_values.iter().enumerate() {
        if p.is_nan() {
            continue;
        }
        // Clamp into (0, 1). The lower bound is the smallest positive double
        // so the small-p branch still represents fully underflowed inputs.
        let p_clamped = p.clamp(P_FLOOR, 1.0 - 1e-15);
        let w = if use_weights {
            weights[i].max(0.0)
        } else {
            1.0
        };
        if w == 0.0 {
            continue;
        }
        // For very small p, tan((0.5-p)*pi) ≈ 1/(p*pi) avoids floating point
        // precision loss where 0.5-p rounds to 0.5. Matches R CCT_pval.cpp.
        if p_clamped < 1e-16 {
            t_sum += w / (p_clamped * PI);
        } else {
            t_sum += w * ((0.5 - p_clamped) * PI).tan();
        }
        w_sum += w;
        count += 1;
    }

    if count == 0 || w_sum == 0.0 {
        return f64::NAN;
    }

    let t = t_sum / w_sum;

    // |T| > T_ASYMP: atan(T)/π saturates at 0.5 in f64. Switch to the
    // upper-tail asymptotic 1/(π·T) (STAAR/R/CCT.R).
    let p = if t > T_ASYMP {
        1.0 / (PI * t)
    } else if t < -T_ASYMP {
        1.0 - 1.0 / (PI * -t)
    } else {
        0.5 - t.atan() / PI
    };

    p.clamp(P_FLOOR, 1.0)
}

/// SKAT p-value via Liu et al. (2009) moment-matching.
///
/// Q ~ Σ λ_j χ²_1. Approximated as a*χ²(l, ncp=δ) + b where
/// parameters (a, l, δ) are matched from the first 4 cumulants.
///
/// Matches SKAT R package: Get_Liu_Params_Mod + Get_Liu_PVal_MOD.
///
/// Cumulants (raw power sums, NOT multiplied by chi-sq moments):
///   c1 = Σ λ_j,  c2 = Σ λ_j²,  c3 = Σ λ_j³,  c4 = Σ λ_j⁴
pub fn mixture_chisq_pvalue(statistic: f64, eigenvalues: &[f64]) -> f64 {
    if eigenvalues.is_empty() || !statistic.is_finite() {
        return f64::NAN;
    }

    // Match SKAT R: Get_Lambda_Org — keep eigenvalues > mean(positive) / 100000
    let positive: Vec<f64> = eigenvalues.iter().copied().filter(|&l| l >= 0.0).collect();
    if positive.is_empty() {
        return 1.0;
    }
    let threshold = positive.iter().sum::<f64>() / positive.len() as f64 / 100000.0;
    let lambdas: Vec<f64> = eigenvalues
        .iter()
        .copied()
        .filter(|&l| l > threshold)
        .collect();
    if lambdas.is_empty() {
        return 1.0;
    }

    // Raw cumulants (power sums of eigenvalues)
    let c1: f64 = lambdas.iter().sum();
    let c2: f64 = lambdas.iter().map(|l| l * l).sum();
    let c3: f64 = lambdas.iter().map(|l| l * l * l).sum();
    let c4: f64 = lambdas.iter().map(|l| l.powi(4)).sum();

    if c2 <= 0.0 {
        return 1.0;
    }

    // Mean and standard deviation of Q under H0
    let mu_q = c1;
    let sigma_q = (2.0 * c2).sqrt();

    // Normalized cumulants (following SKAT R: Get_Liu_Params_Mod)
    let s1 = c3 / c2.powf(1.5);
    let s2 = c4 / (c2 * c2);

    // Fit parameters: Q ≈ a * (χ²(l, δ) - l - δ) * σ_Q + μ_Q
    // Equivalently: standardize Q, then transform to χ²(l, δ)
    let (l, delta, a) = if s1 * s1 > s2 {
        // Noncentral chi-squared approximation
        let a_val = 1.0 / (s1 - (s1 * s1 - s2).sqrt());
        let delta_val = (s1 * a_val.powi(3) - a_val * a_val).max(0.0);
        let l_val = (a_val * a_val - 2.0 * delta_val).max(0.5);
        (l_val, delta_val, a_val)
    } else {
        // Central chi-squared approximation (delta = 0)
        // This is the common case for most genetic data
        let a_val = 1.0 / s1;
        let l_val = 1.0 / (s1 * s1);
        (l_val.max(0.5), 0.0, a_val)
    };

    // Transform: Q_std = (Q - μ_Q)/σ_Q, then χ²_val = Q_std * σ_X + μ_X
    let mu_x = l + delta;
    let sigma_x = (2.0_f64).sqrt() * a;

    let q_std = (statistic - mu_q) / sigma_q;
    let chisq_val = q_std * sigma_x + mu_x;

    if chisq_val <= 0.0 {
        return 1.0;
    }

    // P-value from chi-squared distribution
    // When delta = 0: use central chi-squared (exact for this approximation)
    // When delta > 0: use central chi-squared with adjusted df as fallback
    //   (statrs doesn't have noncentral chi-squared; the central approx is
    //    accurate for small delta, which is the typical genetics case)
    let effective_df = if delta > 0.0 {
        // Noncentral → central approximation: χ²(l, δ) ≈ (l+δ)/(l+2δ) * χ²(l + 2δ)
        // This is the Patnaik (1949) two-moment approximation
        let df_adj = (l + delta).powi(2) / (l + 2.0 * delta);
        let scale = (l + 2.0 * delta) / (l + delta);
        // Rescale chisq_val
        let chisq_rescaled = chisq_val * scale;
        return match ChiSquared::new(df_adj) {
            Ok(dist) => (1.0 - dist.cdf(chisq_rescaled)).max(P_FLOOR),
            Err(_) => f64::NAN,
        };
    } else {
        l
    };

    match ChiSquared::new(effective_df) {
        Ok(dist) => (1.0 - dist.cdf(chisq_val)).max(P_FLOOR),
        Err(_) => f64::NAN,
    }
}

/// Saddlepoint approximation for binary trait score tests.
///
/// Provides calibrated p-values for imbalanced case-control designs where
/// the chi-squared approximation inflates Type I error. Falls back to the
/// normal approximation when p > 0.05 (root-finding is expensive).
///
/// Reference: Dey et al. (2017), "A Fast and Accurate Algorithm to Test
/// for Binary Phenotypes and Its Application to PheWAS"
const P_FILTER: f64 = 0.05;
const ROOT_TOL: f64 = 1e-9;
const ROOT_MAX_ITER: usize = 200;
const EXP_CLAMP: f64 = 500.0;

/// Two-sided SPA p-value for Burden or per-variant score tests.
///
/// `score`: the observed score statistic (w'G'(Y-μ) for Burden, G_j'(Y-μ) for single variant).
/// `mu`: fitted probabilities from logistic null model (length n).
/// `g`: (weighted) genotype vector (length n).
pub fn spa_pvalue(score: f64, mu: &[f64], g: &[f64]) -> f64 {
    debug_assert_eq!(mu.len(), g.len());
    if !score.is_finite() {
        return 1.0;
    }

    let norm = Normal::new(0.0, 1.0).expect("standard normal params are valid");

    // Variance of the centered score S = G'(Y-μ) under H0: Var(S) = Σ μ_i(1-μ_i)g_i²
    let var: f64 = mu
        .iter()
        .zip(g)
        .map(|(&m, &gi)| m * (1.0 - m) * gi * gi)
        .sum();
    if var <= 0.0 {
        return 1.0;
    }

    // Normal pre-filter: skip expensive root-finding for non-significant results.
    // CGF is centered (K'(0) = 0), so mean of S under H0 is 0.
    let z = score / var.sqrt();
    let p_normal = 2.0 * norm.cdf(-z.abs());
    if p_normal > P_FILTER {
        return p_normal;
    }

    // Two-sided: P(|S| >= |score|) = P(S >= |score|) + P(S <= -|score|)
    // Matches R: Saddle(|q|, lower=false) + Saddle(-|q|, lower=true)
    let abs_score = score.abs();
    let p_upper = tail_prob(abs_score, mu, g, &norm);
    let p_lower = 1.0 - tail_prob(-abs_score, mu, g, &norm);

    (p_upper + p_lower).clamp(0.0, 1.0)
}

/// Upper tail probability P(S >= q) via Lugannani-Rice formula.
fn tail_prob(q: f64, mu: &[f64], g: &[f64], norm: &Normal) -> f64 {
    let t_hat = match find_root(q, mu, g) {
        Some(t) if t.abs() > 1e-12 => t,
        _ => {
            // Fallback to normal approximation when root is near zero
            let var: f64 = mu
                .iter()
                .zip(g)
                .map(|(&m, &gi)| m * (1.0 - m) * gi * gi)
                .sum();
            if var <= 0.0 {
                return if q > 0.0 { 0.0 } else { 1.0 };
            }
            return norm.cdf(-q / var.sqrt());
        }
    };

    let k = cgf(t_hat, mu, g);
    let k2 = cgf_d2(t_hat, mu, g);
    if k2 <= 0.0 {
        return 0.5;
    }

    let w_sq = 2.0 * (t_hat * q - k);
    if w_sq < 0.0 {
        return 0.5;
    }

    let w = t_hat.signum() * w_sq.sqrt();
    let v = t_hat * k2.sqrt();
    if w.abs() < 1e-12 {
        return 0.5;
    }

    (norm.cdf(-w) + norm.pdf(w) * (1.0 / w - 1.0 / v)).clamp(0.0, 1.0)
}

// Cumulant generating function and derivatives for centered score S = G'(Y-μ)
//   K(t)  = Σ [-t μ_i g_i + log(1 - μ_i + μ_i exp(t g_i))]
//   K'(t) = Σ [-μ_i g_i + μ_i g_i exp(t g_i) / (1 - μ_i + μ_i exp(t g_i))]
//   K''(t)= Σ μ_i(1-μ_i) g_i² exp(t g_i) / (1 - μ_i + μ_i exp(t g_i))²
// The -t·μ·g centering ensures K'(0) = 0, matching R's K_Binary_SPA.

fn cgf(t: f64, mu: &[f64], g: &[f64]) -> f64 {
    mu.iter()
        .zip(g)
        .map(|(&m, &gi)| {
            -t * m * gi + (1.0 - m + m * (t * gi).clamp(-EXP_CLAMP, EXP_CLAMP).exp()).ln()
        })
        .sum()
}

fn cgf_d1(t: f64, mu: &[f64], g: &[f64]) -> f64 {
    mu.iter()
        .zip(g)
        .map(|(&m, &gi)| {
            let e = (t * gi).clamp(-EXP_CLAMP, EXP_CLAMP).exp();
            -m * gi + m * gi * e / (1.0 - m + m * e)
        })
        .sum()
}

fn cgf_d2(t: f64, mu: &[f64], g: &[f64]) -> f64 {
    mu.iter()
        .zip(g)
        .map(|(&m, &gi)| {
            let e = (t * gi).clamp(-EXP_CLAMP, EXP_CLAMP).exp();
            let d = 1.0 - m + m * e;
            m * (1.0 - m) * gi * gi * e / (d * d)
        })
        .sum()
}

/// Bisection root-finding for K'(t*) = q.
///
/// K'(t) is monotonically increasing so bisection always converges.
fn find_root(q: f64, mu: &[f64], g: &[f64]) -> Option<f64> {
    let k1_0 = cgf_d1(0.0, mu, g);
    let diff = q - k1_0;
    if diff.abs() < ROOT_TOL {
        return None;
    }

    let going_right = diff > 0.0;
    let (mut lo, mut hi) = if going_right { (0.0, 1.0) } else { (-1.0, 0.0) };

    for _ in 0..50 {
        if going_right {
            if cgf_d1(hi, mu, g) >= q {
                break;
            }
            hi *= 2.0;
        } else {
            if cgf_d1(lo, mu, g) <= q {
                break;
            }
            lo *= 2.0;
        }
    }

    for _ in 0..ROOT_MAX_ITER {
        let mid = 0.5 * (lo + hi);
        if (hi - lo).abs() < ROOT_TOL {
            return Some(mid);
        }
        if cgf_d1(mid, mu, g) < q {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cauchy_combine_basic() {
        let p = cauchy_combine(&[0.01, 0.05, 0.5]);
        assert!(p > 0.0 && p < 0.05, "Combined p should be small: {p}");
    }

    #[test]
    fn test_cauchy_combine_all_large() {
        let p = cauchy_combine(&[0.5, 0.5, 0.5]);
        assert!((p - 0.5).abs() < 0.01, "All 0.5 should give ~0.5: {p}");
    }

    #[test]
    fn test_cauchy_combine_extreme() {
        let p = cauchy_combine(&[1e-10, 0.99]);
        assert!(p < 0.01, "One very small p should dominate: {p}");
    }

    #[test]
    fn test_cauchy_combine_empty() {
        assert!(cauchy_combine(&[]).is_nan());
    }

    #[test]
    fn test_weighted_vs_equal() {
        let pvals = [0.01, 0.05, 0.5];
        let equal = cauchy_combine(&pvals);
        let weighted = cauchy_combine_weighted(&pvals, &[1.0, 1.0, 1.0]);
        assert!((equal - weighted).abs() < 1e-12);
    }

    #[test]
    fn test_cauchy_combine_extreme_underflow() {
        // Many ultra-small p-values: naive formula returns 0.0 because
        // atan(T)/π saturates at 0.5 in f64. Asymptotic branch must return
        // a positive denormal that respects the harmonic-mean asymptotic.
        let pvals = vec![1e-200; 6];
        let p = cauchy_combine(&pvals);
        assert!(p > 0.0, "must not underflow to zero: {p:e}");
        assert!(p.is_finite(), "must be finite: {p:e}");
        // For equal weights and p_i = p0, asymptotic gives p ≈ p0 (harmonic
        // mean of K identical values is the value itself).
        assert!((p - 1e-200).abs() < 1e-205, "expected ~1e-200, got {p:e}");
    }

    #[test]
    fn test_cauchy_combine_monotone_in_tail() {
        // Combining smaller p-values must yield a smaller (or equal) result.
        let p_a = cauchy_combine(&[1e-100, 1e-100, 1e-100]);
        let p_b = cauchy_combine(&[1e-200, 1e-200, 1e-200]);
        assert!(
            p_b < p_a,
            "smaller inputs must give smaller p: {p_a:e} vs {p_b:e}"
        );
    }

    #[test]
    fn test_cauchy_combine_floors_at_denormal() {
        // Even at absolute floor (one input is 0.0), result must be positive
        // and at least the smallest representable positive double.
        let p = cauchy_combine(&[0.0, 0.5, 0.5]);
        assert!(p > 0.0, "must floor to >0: {p:e}");
        assert!(p >= f64::MIN_POSITIVE * f64::EPSILON);
    }

    #[test]
    fn test_cauchy_combine_harmonic_asymptotic() {
        // Mixed extreme p-values: result should approximate harmonic-mean
        // form K / sum(1/p_i) when all inputs are << 1e-16.
        let pvals = [1e-50, 2e-80, 3e-100];
        let p = cauchy_combine(&pvals);
        let expected = pvals.len() as f64 / pvals.iter().map(|p| 1.0 / p).sum::<f64>();
        let rel_err = (p - expected).abs() / expected;
        assert!(
            rel_err < 1e-3,
            "p={p:e} expected={expected:e} rel_err={rel_err:e}"
        );
    }

    #[test]
    fn test_weighted_emphasis() {
        let pvals = [0.01, 0.99];
        // High weight on the small p-value
        let p_high = cauchy_combine_weighted(&pvals, &[10.0, 1.0]);
        // High weight on the large p-value
        let p_low = cauchy_combine_weighted(&pvals, &[1.0, 10.0]);
        assert!(
            p_high < p_low,
            "Weighting small p should give smaller combined p"
        );
    }

    #[test]
    fn test_single_eigenvalue() {
        // Q ~ λ χ²(1). Exact: P(Q > t) = P(χ²(1) > t/λ)
        let lambda = 2.0;
        let t = 5.0;
        let p = mixture_chisq_pvalue(t, &[lambda]);
        let exact = 1.0 - ChiSquared::new(1.0).unwrap().cdf(t / lambda);
        assert!(
            (p - exact).abs() < 0.01,
            "Single eigenvalue: got {p:.6}, expected {exact:.6}"
        );
    }

    #[test]
    fn test_equal_eigenvalues() {
        // Q ~ λ * χ²(k). Exact: P(Q > t) = P(χ²(k) > t/λ)
        let lambda = 1.5;
        let k = 5;
        let eigenvalues = vec![lambda; k];
        let t = 10.0;
        let p = mixture_chisq_pvalue(t, &eigenvalues);
        let exact = 1.0 - ChiSquared::new(k as f64).unwrap().cdf(t / lambda);
        assert!(
            (p - exact).abs() < 0.05,
            "Equal eigenvalues: got {p:.6}, expected {exact:.6}"
        );
    }

    #[test]
    fn test_null_statistic() {
        let p = mixture_chisq_pvalue(0.0, &[1.0, 2.0, 3.0]);
        assert!(p > 0.5, "Zero statistic → large p: {p}");
    }

    #[test]
    fn test_uniform_under_null() {
        // Under H0, p-values should be approximately uniform
        // Generate Q = Σ λ_j z_j² where z ~ N(0,1) (under null)
        // For a quick check: Q = c1 (the mean) should give p ≈ 0.5
        let eigenvalues = [1.0, 2.0, 3.0];
        let mean = eigenvalues.iter().sum::<f64>();
        let p = mixture_chisq_pvalue(mean, &eigenvalues);
        assert!(p > 0.3 && p < 0.7, "p at mean should be ~0.5: {p}");
    }

    #[test]
    fn cgf_at_zero_is_zero() {
        let mu = vec![0.5, 0.3, 0.8];
        let g = vec![1.0, 0.0, 2.0];
        assert!(cgf(0.0, &mu, &g).abs() < 1e-12);
    }

    #[test]
    fn k1_at_zero_is_zero() {
        // Centered CGF: K'(0) = 0 (the centering term -t*mu*g makes this hold)
        let mu = vec![0.5, 0.3, 0.8];
        let g = vec![1.0, 0.0, 2.0];
        assert!(
            cgf_d1(0.0, &mu, &g).abs() < 1e-12,
            "K'(0) should be 0 for centered CGF"
        );
    }

    #[test]
    fn k2_at_zero_equals_variance() {
        let mu = vec![0.5, 0.3, 0.8];
        let g = vec![1.0, 0.0, 2.0];
        let expected: f64 = mu
            .iter()
            .zip(&g)
            .map(|(m, gi)| m * (1.0 - m) * gi * gi)
            .sum();
        assert!((cgf_d2(0.0, &mu, &g) - expected).abs() < 1e-12);
    }

    #[test]
    fn balanced_matches_normal() {
        // For balanced case-control (mu=0.5), SPA should match normal approximation
        let mu = vec![0.5; 100];
        let g: Vec<f64> = (0..100).map(|i| if i < 10 { 1.0 } else { 0.0 }).collect();
        let score = 3.0; // centered score: sum(G_i * (Y_i - mu_i))
        let p_spa = spa_pvalue(score, &mu, &g);

        let norm = Normal::new(0.0, 1.0).unwrap();
        let var: f64 = mu
            .iter()
            .zip(&g)
            .map(|(m, gi)| m * (1.0 - m) * gi * gi)
            .sum();
        let z = score / var.sqrt(); // centered CGF: E[S] = 0
        let p_normal = 2.0 * norm.cdf(-z.abs());

        assert!(
            (p_spa - p_normal).abs() < 0.05,
            "Balanced: SPA={p_spa:.6}, normal={p_normal:.6}"
        );
    }

    #[test]
    fn score_zero_gives_large_p() {
        // Centered score of 0 means no evidence against H0 → large p-value
        let mu = vec![0.1; 50];
        let g = vec![1.0; 50];
        let p = spa_pvalue(0.0, &mu, &g);
        assert!(p > 0.9, "Score=0 should give large p: {p}");
    }

    #[test]
    fn extreme_imbalance_does_not_panic() {
        let mu = vec![0.001; 1000];
        let g: Vec<f64> = (0..1000).map(|i| if i < 5 { 1.0 } else { 0.0 }).collect();
        let score = 4.0;
        let p = spa_pvalue(score, &mu, &g);
        assert!((0.0..=1.0).contains(&p), "p must be valid: {p}");
    }

    #[test]
    fn zero_genotypes_gives_one() {
        let mu = vec![0.3; 50];
        let g = vec![0.0; 50];
        assert!((spa_pvalue(0.0, &mu, &g) - 1.0).abs() < 1e-10);
    }
}
