use faer::Mat;
use statrs::distribution::{Beta, Continuous};
use statrs::function::erf::erfc;

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
    /// channel order matches `STAAR_PHRED_CHANNELS` in `column.rs`
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
/// engine. When `use_spa` is true, the trait is binary, AND fitted values
/// are present, Burden and ACAT-V use saddlepoint approximation; SKAT
/// always uses moment-matching. Fallback to the sum-stat path when SPA is
/// requested without fitted values.
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
    let n = g.nrows();
    let sigma2 = null.sigma2;
    let spa_mu = if use_spa {
        null.fitted_values.as_deref()
    } else {
        None
    };

    // Decide SPA-or-not once here instead of pattern-matching inside a
    // closure that fires ~68 times per gene.
    match spa_mu {
        Some(mu) => staar_tests(
            &u,
            &k,
            annotation_matrix,
            mafs,
            sigma2,
            |w| burden_spa(g, &u, w, mu),
            |w_acat, _w_burden| acat_v_spa(g, &u, w_acat, mu),
        ),
        None => staar_tests(
            &u,
            &k,
            annotation_matrix,
            mafs,
            sigma2,
            |w| burden(&u, &k, w, sigma2),
            |w_acat, w_burden| acat_v(&u, &k, w_acat, w_burden, mafs, n, sigma2),
        ),
    }
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
/// `w_acat`: ACAT weights (Cauchy combination weights only).
/// `w_burden`: Burden weights (used for both numerator and denominator of
/// the rare-group score statistic so it's χ²(1) under H0).
/// Matches R STAAR_O.cpp lines 199-216: per-variant common test stat is
/// `x²/Cov` (no weights in stat); rare-group test stat is `(Σ x·w_B)² /
/// Σ Covw[w_B]`. `w_A` is only a Cauchy combination weight.
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

    // Group very rare variants: Burden-weighted score AND Burden-weighted
    // covariance — both numerator and denominator must share weights for the
    // statistic to be χ²(1) under H0. Matches R STAAR_O.cpp lines 199-216:
    // sum0 = Σ x·w_B, sumx = sum0² / Σ Covw[w_B], pchisq(sumx, 1). w_A enters
    // only as the Cauchy combination weight `sumw/n0` on the resulting p-value.
    if !rare_indices.is_empty() {
        let wu: f64 = rare_indices.iter().map(|&j| w_burden[j] * u[(j, 0)]).sum();
        let mut wkw = 0.0;
        for &j in &rare_indices {
            for &l in &rare_indices {
                wkw += w_burden[j] * k[(j, l)] * w_burden[l];
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

/// χ²(1) survival function = `erfc(sqrt(t/2))`.
///
/// Direct closed form rather than `1 - ChiSquared::cdf(t)`. statrs's gamma
/// CDF saturates at 1.0 for t ≳ 1400 (so the survival becomes literal 0.0
/// for any larger statistic), whereas `erfc` retains precision out to
/// t ≈ 2700. Past that we floor at the smallest positive double so the
/// downstream Cauchy combination sees a strictly positive p-value instead
/// of poisoning the omnibus with an exact zero.
pub(crate) fn chisq1_pvalue(t: f64) -> f64 {
    if t <= 0.0 || !t.is_finite() {
        return 1.0;
    }
    let p = erfc((t * 0.5).sqrt());
    if p > 0.0 { p } else { stats::P_FLOOR }
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
    staar_tests(
        u,
        k,
        annotation_matrix,
        mafs,
        1.0,
        |w| burden(u, k, w, 1.0),
        |w_acat, w_burden| acat_v(u, k, w_acat, w_burden, mafs, n_total, 1.0),
    )
}

/// Per-gene scratch for the STAAR test battery: weight vectors (length m),
/// SKAT kernel buffer (m × m), and per-gene omnibus accumulators. Collects
/// the allocations that used to be scattered across `staar_tests` so they
/// can be sized once per gene instead of interleaved with the test loop.
struct GeneScratch {
    beta_1_25: Vec<f64>,
    beta_1_1: Vec<f64>,
    acat_denom: Vec<f64>,
    wa_base_1_25: Vec<f64>,
    wa_base_1_1: Vec<f64>,
    wb_1_25: Vec<f64>,
    wb_1_1: Vec<f64>,
    ws_1_25: Vec<f64>,
    ws_1_1: Vec<f64>,
    wa_1_25: Vec<f64>,
    wa_1_1: Vec<f64>,
    kernel_buf: Mat<f64>,
    by_test: [Vec<f64>; 6],
    per_annotation: Vec<[f64; 6]>,
    all_p: Vec<f64>,
}

impl GeneScratch {
    fn with_capacity(m: usize, n_channels: usize) -> Self {
        let mzeros = || vec![0.0; m];
        Self {
            beta_1_25: mzeros(),
            beta_1_1: mzeros(),
            acat_denom: mzeros(),
            wa_base_1_25: mzeros(),
            wa_base_1_1: mzeros(),
            wb_1_25: mzeros(),
            wb_1_1: mzeros(),
            ws_1_25: mzeros(),
            ws_1_1: mzeros(),
            wa_1_25: mzeros(),
            wa_1_1: mzeros(),
            kernel_buf: Mat::zeros(m, m),
            by_test: std::array::from_fn(|_| Vec::with_capacity(1 + n_channels)),
            per_annotation: Vec::with_capacity(n_channels),
            all_p: Vec::with_capacity(6 + n_channels * 6),
        }
    }
}

/// Shared test engine for both single-study and meta-analysis paths.
/// Computes all 6 base tests, annotation-weighted variants, and omnibus combinations.
///
/// Burden and ACAT-V dispatch is passed in as closures so the SPA vs
/// sum-stat branch is decided exactly once — at the caller — instead of
/// pattern-matching per channel per call.
fn staar_tests<BF, AF>(
    u: &Mat<f64>,
    k: &Mat<f64>,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
    sigma2: f64,
    run_burden: BF,
    run_acat_v: AF,
) -> StaarResult
where
    BF: Fn(&[f64]) -> f64,
    AF: Fn(&[f64], &[f64]) -> f64,
{
    let m = mafs.len();
    let n_channels = annotation_matrix.len();
    let mut s = GeneScratch::with_capacity(m, n_channels);

    for (j, &maf) in mafs.iter().enumerate() {
        s.beta_1_25[j] = beta_density_weight(maf, 1.0, 25.0);
        s.beta_1_1[j] = beta_density_weight(maf, 1.0, 1.0);
        let d = beta_density_weight(maf, 0.5, 0.5);
        s.acat_denom[j] = if d > 0.0 { d * d } else { 1.0 };
        s.wa_base_1_25[j] = s.beta_1_25[j] * s.beta_1_25[j] / s.acat_denom[j];
        s.wa_base_1_1[j] = s.beta_1_1[j] * s.beta_1_1[j] / s.acat_denom[j];
    }

    let base_burden_1_25 = run_burden(&s.beta_1_25);
    let base_burden_1_1 = run_burden(&s.beta_1_1);
    let base_skat_1_25 = skat(u, k, &s.beta_1_25, sigma2, &mut s.kernel_buf);
    let base_skat_1_1 = skat(u, k, &s.beta_1_1, sigma2, &mut s.kernel_buf);
    let base_acat_v_1_25 = run_acat_v(&s.wa_base_1_25, &s.beta_1_25);
    let base_acat_v_1_1 = run_acat_v(&s.wa_base_1_1, &s.beta_1_1);

    let acat_o = stats::cauchy_combine(&[
        base_burden_1_25,
        base_burden_1_1,
        base_skat_1_25,
        base_skat_1_1,
        base_acat_v_1_25,
        base_acat_v_1_1,
    ]);

    // Accumulate p-values per test type across annotation channels for STAAR omnibus.
    // Index order: Burden(1,25), Burden(1,1), SKAT(1,25), SKAT(1,1), ACAT-V(1,25), ACAT-V(1,1)
    s.by_test[0].push(base_burden_1_25);
    s.by_test[1].push(base_burden_1_1);
    s.by_test[2].push(base_skat_1_25);
    s.by_test[3].push(base_skat_1_1);
    s.by_test[4].push(base_acat_v_1_25);
    s.by_test[5].push(base_acat_v_1_1);

    for channel_weights in annotation_matrix {
        // Eight parallel slices indexed by j; a single range loop beats chaining
        // seven zips. Keep the indexing form; clippy's range-loop lint disagrees.
        #[allow(clippy::needless_range_loop)]
        for j in 0..m {
            let a = channel_weights[j];
            let a_sqrt = a.sqrt();
            s.wb_1_25[j] = s.beta_1_25[j] * a;
            s.wb_1_1[j] = s.beta_1_1[j] * a;
            s.ws_1_25[j] = s.beta_1_25[j] * a_sqrt;
            s.ws_1_1[j] = s.beta_1_1[j] * a_sqrt;
            s.wa_1_25[j] = a * s.beta_1_25[j] * s.beta_1_25[j] / s.acat_denom[j];
            s.wa_1_1[j] = a * s.beta_1_1[j] * s.beta_1_1[j] / s.acat_denom[j];
        }

        let p = [
            run_burden(&s.wb_1_25),
            run_burden(&s.wb_1_1),
            skat(u, k, &s.ws_1_25, sigma2, &mut s.kernel_buf),
            skat(u, k, &s.ws_1_1, sigma2, &mut s.kernel_buf),
            run_acat_v(&s.wa_1_25, &s.wb_1_25),
            run_acat_v(&s.wa_1_1, &s.wb_1_1),
        ];
        for (bucket, pv) in s.by_test.iter_mut().zip(p.iter()) {
            bucket.push(*pv);
        }
        s.per_annotation.push(p);
    }

    let staar_b_1_25 = stats::cauchy_combine(&s.by_test[0]);
    let staar_b_1_1 = stats::cauchy_combine(&s.by_test[1]);
    let staar_s_1_25 = stats::cauchy_combine(&s.by_test[2]);
    let staar_s_1_1 = stats::cauchy_combine(&s.by_test[3]);
    let staar_a_1_25 = stats::cauchy_combine(&s.by_test[4]);
    let staar_a_1_1 = stats::cauchy_combine(&s.by_test[5]);

    s.all_p.extend_from_slice(&[
        base_burden_1_25,
        base_burden_1_1,
        base_skat_1_25,
        base_skat_1_1,
        base_acat_v_1_25,
        base_acat_v_1_1,
    ]);
    for p in &s.per_annotation {
        s.all_p.extend_from_slice(p);
    }
    let staar_o = stats::cauchy_combine(&s.all_p);

    StaarResult {
        burden_1_25: base_burden_1_25,
        burden_1_1: base_burden_1_1,
        skat_1_25: base_skat_1_25,
        skat_1_1: base_skat_1_1,
        acat_v_1_25: base_acat_v_1_25,
        acat_v_1_1: base_acat_v_1_1,
        per_annotation: s.per_annotation,
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
