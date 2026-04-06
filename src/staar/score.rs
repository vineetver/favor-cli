use faer::Mat;
use statrs::distribution::{Beta, ChiSquared, Continuous, ContinuousCDF};

use super::model::NullModel;
use super::stats;

/// Full STAAR results for one gene: all tests × all annotation channels,
/// plus per-test omnibus and overall omnibus.
///
/// Matches STAARpipeline R naming:
///   Base tests:      Burden(1,25), SKAT(1,1), ACAT-V(1,25), etc.
///   Per-annotation:  Burden(1,25)-CADD, Burden(1,25)-LINSIGHT, etc.
///   Per-test omni:   STAAR-B(1,25), STAAR-S(1,1), STAAR-A(1,25), etc.
///   Cross-test omni: ACAT-O, STAAR-O
#[derive(Debug, Clone)]
pub struct StaarResult {
    /// 6 base test p-values (no annotation weighting, just MAF-based beta weights)
    pub burden_1_25: f64,
    pub burden_1_1: f64,
    pub skat_1_25: f64,
    pub skat_1_1: f64,
    pub acat_v_1_25: f64,
    pub acat_v_1_1: f64,

    /// Per-annotation p-values: [channel][test_index]
    /// test_index: 0=Burden(1,25) 1=Burden(1,1) 2=SKAT(1,25) 3=SKAT(1,1) 4=ACAT-V(1,25) 5=ACAT-V(1,1)
    /// channel order matches AnnotationWeights::DISPLAY_NAMES
    pub per_annotation: Vec<[f64; 6]>,

    /// Per-test omnibus: Cauchy across all annotation channels for one test type
    /// Matches R's STAAR-B(1,25), STAAR-B(1,1), STAAR-S(1,25), etc.
    pub staar_b_1_25: f64,
    pub staar_b_1_1: f64,
    pub staar_s_1_25: f64,
    pub staar_s_1_1: f64,
    pub staar_a_1_25: f64,
    pub staar_a_1_1: f64,

    /// ACAT-O: Cauchy of the 6 base test p-values (no annotation)
    pub acat_o: f64,

    /// STAAR-O: Cauchy of ALL test × ALL annotation p-values (the omnibus)
    pub staar_o: f64,
}

/// Run full STAAR analysis for one gene from raw genotypes.
///
/// Computes U = G'r and K = G'(I-H)G, then delegates to the shared test
/// engine. When `use_spa` is true and the trait is binary, Burden and ACAT-V
/// use saddlepoint approximation; SKAT always uses moment-matching.
pub fn run_staar(
    g: &Mat<f64>,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
    null: &NullModel,
    use_spa: bool,
) -> StaarResult {
    let m = g.ncols();
    if m == 0 {
        return nan_result();
    }

    let u = g.transpose() * &null.residuals;
    let k = null.compute_kernel(g);
    let spa_mu = if use_spa {
        null.fitted_values.as_deref()
    } else {
        None
    };

    staar_tests(
        &u,
        &k,
        annotation_matrix,
        mafs,
        null.sigma2,
        g.nrows(),
        Some(g),
        spa_mu,
    )
}

/// Burden test with SPA: weighted genotype per sample, then saddlepoint p-value.
fn burden_spa(g: &Mat<f64>, u: &Mat<f64>, w: &[f64], mu: &[f64]) -> f64 {
    let n = g.nrows();
    let m = w.len();

    let mut w_g = vec![0.0; n];
    for j in 0..m {
        if w[j] == 0.0 {
            continue;
        }
        for i in 0..n {
            w_g[i] += w[j] * g[(i, j)];
        }
    }

    let score: f64 = (0..m).map(|j| w[j] * u[(j, 0)]).sum();
    stats::spa_pvalue(score, mu, &w_g)
}

/// ACAT-V with SPA: per-variant saddlepoint p-values, Cauchy combined.
fn acat_v_spa(g: &Mat<f64>, u: &Mat<f64>, w: &[f64], mu: &[f64]) -> f64 {
    let n = g.nrows();
    let m = w.len();
    let mut p_values = Vec::with_capacity(m);
    let mut cauchy_weights = Vec::with_capacity(m);
    let mut g_col = vec![0.0; n];

    for j in 0..m {
        if w[j] == 0.0 {
            continue;
        }
        for i in 0..n {
            g_col[i] = g[(i, j)];
        }
        p_values.push(stats::spa_pvalue(u[(j, 0)], mu, &g_col));
        cauchy_weights.push(w[j]);
    }

    if p_values.is_empty() {
        return 1.0;
    }
    stats::cauchy_combine_weighted(&p_values, &cauchy_weights)
}

fn burden(u: &Mat<f64>, k: &Mat<f64>, w: &[f64], sigma2: f64) -> f64 {
    let m = w.len();
    let wu: f64 = (0..m).map(|j| w[j] * u[(j, 0)]).sum();
    let mut wkw = 0.0;
    for j in 0..m {
        if w[j] == 0.0 {
            continue;
        }
        for l in 0..m {
            wkw += w[j] * k[(j, l)] * w[l];
        }
    }
    let var = sigma2 * wkw;
    if var <= 0.0 {
        return 1.0;
    }
    chisq1_pvalue(wu * wu / var)
}

fn skat(u: &Mat<f64>, k: &Mat<f64>, w: &[f64], sigma2: f64, kernel: &mut Mat<f64>) -> f64 {
    let m = w.len();
    let q: f64 = (0..m).map(|j| w[j] * w[j] * u[(j, 0)] * u[(j, 0)]).sum();
    for j in 0..m {
        for l in 0..m {
            kernel[(j, l)] = 0.0;
        }
    }
    for j in 0..m {
        if w[j] == 0.0 {
            continue;
        }
        for l in j..m {
            let val = sigma2 * w[j] * k[(j, l)] * w[l];
            kernel[(j, l)] = val;
            kernel[(l, j)] = val;
        }
    }
    let eigenvalues = symmetric_eigenvalues(kernel);
    stats::mixture_chisq_pvalue(q, &eigenvalues)
}

/// MAC threshold below which variants are grouped into a single Burden-like
/// statistic before Cauchy combination. Matches R STAARpipeline's mac_thres=10.
const MAC_THRESHOLD: f64 = 10.0;

/// ACAT-V with MAC-based grouping of very rare variants.
///
/// `w_acat`: ACAT-V weights (for Cauchy combination weights and per-variant variance).
/// `w_burden`: corresponding Burden weights (for grouped rare variant score).
/// Matches R STAAR_O_SMMAT.cpp: Burden weights for grouped score, ACAT weights
/// for the covariance denominator and Cauchy combination.
fn acat_v(
    u: &Mat<f64>,
    k: &Mat<f64>,
    w_acat: &[f64],
    w_burden: &[f64],
    mafs: &[f64],
    n_samples: usize,
    sigma2: f64,
) -> f64 {
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
            // Common: individual chi-sq(1) with unweighted score variance
            let var_j = sigma2 * k[(j, j)];
            if var_j <= 0.0 {
                continue;
            }
            let t = u[(j, 0)] * u[(j, 0)] / var_j;
            p_values.push(chisq1_pvalue(t));
            cauchy_weights.push(w_acat[j]);
        } else {
            rare_indices.push(j);
        }
    }

    // Group very rare variants: Burden-weighted score, ACAT-weighted covariance
    if !rare_indices.is_empty() {
        let wu: f64 = rare_indices.iter().map(|&j| w_burden[j] * u[(j, 0)]).sum();
        let mut wkw = 0.0;
        for &j in &rare_indices {
            for &l in &rare_indices {
                wkw += w_acat[j] * k[(j, l)] * w_acat[l];
            }
        }
        let var = sigma2 * wkw;
        if var > 0.0 {
            p_values.push(chisq1_pvalue(wu * wu / var));
            let mean_w: f64 =
                rare_indices.iter().map(|&j| w_acat[j]).sum::<f64>() / rare_indices.len() as f64;
            cauchy_weights.push(mean_w);
        }
    }

    if p_values.is_empty() {
        return 1.0;
    }
    stats::cauchy_combine_weighted(&p_values, &cauchy_weights)
}

fn chisq1_pvalue(t: f64) -> f64 {
    if t <= 0.0 || !t.is_finite() {
        return 1.0;
    }
    match ChiSquared::new(1.0) {
        Ok(dist) => 1.0 - dist.cdf(t),
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

/// Run STAAR from pre-computed U and K (for MetaSTAAR).
///
/// Both U and K must be pre-scaled by 1/σ² (matching R MetaSTAAR convention).
/// The test functions receive sigma2=1.0 since the scaling is already applied.
/// SPA is unavailable without raw genotypes; uses chi-squared throughout.
/// `n_total` is the merged sample count across studies — needed for MAC-based
/// ACAT-V grouping of very rare variants.
pub fn run_staar_from_sumstats(
    u: &Mat<f64>,
    k: &Mat<f64>,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
    n_total: usize,
) -> StaarResult {
    if u.nrows() == 0 {
        return nan_result();
    }
    staar_tests(u, k, annotation_matrix, mafs, 1.0, n_total, None, None)
}

/// Shared test engine for both single-study and meta-analysis paths.
/// Computes all 6 base tests, annotation-weighted variants, and omnibus combinations.
#[allow(clippy::too_many_arguments)]
fn staar_tests(
    u: &Mat<f64>,
    k: &Mat<f64>,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
    sigma2: f64,
    n_samples: usize,
    g_for_spa: Option<&Mat<f64>>,
    spa_mu: Option<&[f64]>,
) -> StaarResult {
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

    let run_burden = |w: &[f64]| -> f64 {
        match (g_for_spa, spa_mu) {
            (Some(g), Some(mu)) => burden_spa(g, u, w, mu),
            _ => burden(u, k, w, sigma2),
        }
    };
    let run_acat_v = |w_acat: &[f64], w_burden: &[f64]| -> f64 {
        match (g_for_spa, spa_mu) {
            (Some(g), Some(mu)) => acat_v_spa(g, u, w_acat, mu),
            _ => acat_v(u, k, w_acat, w_burden, mafs, n_samples, sigma2),
        }
    };

    let m = mafs.len();
    let mut kernel_buf = Mat::zeros(m, m);

    let base_burden_1_25 = run_burden(&beta_1_25);
    let base_burden_1_1 = run_burden(&beta_1_1);
    let base_skat_1_25 = skat(u, k, &beta_1_25, sigma2, &mut kernel_buf);
    let base_skat_1_1 = skat(u, k, &beta_1_1, sigma2, &mut kernel_buf);
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

    // Accumulate p-values per test type across annotation channels for STAAR omnibus.
    // Index order: Burden(1,25), Burden(1,1), SKAT(1,25), SKAT(1,1), ACAT-V(1,25), ACAT-V(1,1)
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
            skat(u, k, &ws_1_25, sigma2, &mut kernel_buf),
            skat(u, k, &ws_1_1, sigma2, &mut kernel_buf),
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

// ═══════════════════════════════════════════════════════════════════════
// Beta density weights
// ═══════════════════════════════════════════════════════════════════════

/// Beta density weight: dbeta(maf, a1, a2).
///
/// For beta(1,25): upweights ultra-rare variants (MAF near 0).
/// For beta(1,1): uniform weight = 1.0 for all MAFs in (0,1).
///
/// Returns 0.0 for maf=0 or maf≥0.5 (monomorphic or not minor).
pub fn beta_density_weight(maf: f64, a1: f64, a2: f64) -> f64 {
    if maf <= 0.0 || maf >= 0.5 || !maf.is_finite() {
        return 0.0;
    }
    if (a1 - 1.0).abs() < 1e-10 && (a2 - 1.0).abs() < 1e-10 {
        return 1.0;
    }
    match Beta::new(a1, a2) {
        Ok(dist) => dist.pdf(maf),
        Err(_) => 0.0,
    }
}
