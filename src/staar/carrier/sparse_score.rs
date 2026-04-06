//! Sparse-native STAAR score computation from carrier lists.
//!
//! Computes G'r, G'X, and G'G directly from carrier lists in O(total_MAC)
//! instead of O(n_samples × n_variants). For a gene with 20 variants at
//! MAC=5 in 50K samples, this is ~8,000× less work.
//!
//! The speedup is proportional to n_samples/MAC — it GROWS with cohort size.

use std::collections::HashMap;

use faer::Mat;

use super::reader::{CarrierEntry, CarrierList};
use crate::staar::model::NullModel;
use crate::staar::score::{self, StaarResult};

/// Pre-parsed analysis vectors from the null model, laid out for
/// cache-friendly carrier-indexed access. Computed once, shared
/// read-only across all genes and threads.
///
/// Arrays are COMPACT at `n_pheno` size (only samples with phenotype data).
/// `vcf_to_pheno` maps VCF sample indices → phenotype array indices.
/// Samples without phenotype are skipped in scoring, not multiplied by zero.
pub struct AnalysisVectors {
    /// Residuals from null model. Length = `n_pheno`.
    pub residuals: Vec<f64>,
    /// Covariate matrix, ROW-major: `x[i * k .. (i+1) * k]`. Length = `n_pheno * k`.
    pub x_row_major: Vec<f64>,
    /// `(X'X)^{-1}` from null model, row-major `[k × k]`.
    pub xtx_inv: Vec<f64>,
    /// σ² from null model.
    pub sigma2: f64,
    /// Number of covariates (including intercept).
    pub k: usize,
    /// Number of samples with phenotype data.
    pub n_pheno: usize,
    /// Total number of VCF samples (for dense matrix construction, e.g. SPA).
    pub n_vcf_total: usize,
    /// Working weights W_i = μ_i(1-μ_i) for binary traits. Empty for continuous.
    /// Length = `n_pheno` when present.
    pub working_weights: Vec<f64>,
    /// Fitted values μ_i for SPA. Empty if SPA not used.
    /// Length = `n_pheno` when present.
    pub fitted_values: Vec<f64>,
    /// Maps VCF sample index → phenotype index. `None` = no phenotype (skip).
    /// Length = `n_vcf_total`.
    pub vcf_to_pheno: Vec<Option<u32>>,
}

impl AnalysisVectors {
    /// Build from a fitted NullModel with compact arrays and VCF→pheno remapping.
    /// Arrays stay at `n_pheno` size. No zero-padding.
    pub fn from_null_model(null: &NullModel, pheno_mask: &[bool]) -> Self {
        let n_pheno = null.n_samples;
        let n_vcf_total = pheno_mask.len();
        let k = null.x_matrix.ncols();

        // Compact arrays at n_pheno size
        let mut residuals = Vec::with_capacity(n_pheno);
        let mut x_row_major = Vec::with_capacity(n_pheno * k);
        let mut working_weights = if null.working_weights.is_some() {
            Vec::with_capacity(n_pheno)
        } else {
            Vec::new()
        };
        let mut fitted_values = if null.fitted_values.is_some() {
            Vec::with_capacity(n_pheno)
        } else {
            Vec::new()
        };

        // Build VCF→pheno mapping
        let mut vcf_to_pheno = vec![None; n_vcf_total];
        let mut j = 0u32;
        for (vcf_idx, &has_pheno) in pheno_mask.iter().enumerate() {
            if has_pheno {
                vcf_to_pheno[vcf_idx] = Some(j);
                let ji = j as usize;
                residuals.push(null.residuals[(ji, 0)]);
                for c in 0..k {
                    x_row_major.push(null.x_matrix[(ji, c)]);
                }
                if let Some(ref ww) = null.working_weights {
                    working_weights.push(ww[ji]);
                }
                if let Some(ref fv) = null.fitted_values {
                    fitted_values.push(fv[ji]);
                }
                j += 1;
            }
        }
        debug_assert_eq!(j as usize, n_pheno);

        // Row-major (X'X)^{-1}
        let mut xtx_inv = vec![0.0; k * k];
        for i in 0..k {
            for c in 0..k {
                xtx_inv[i * k + c] = null.xtx_inv[(i, c)];
            }
        }

        Self {
            residuals,
            x_row_major,
            xtx_inv,
            sigma2: null.sigma2,
            k,
            n_pheno,
            n_vcf_total,
            working_weights,
            fitted_values,
            vcf_to_pheno,
        }
    }

    pub fn is_binary(&self) -> bool {
        !self.working_weights.is_empty()
    }
}

/// Score one gene using sparse carrier-indexed computation.
///
/// Returns (U/σ², K/σ²) — pre-scaled to the MetaSTAAR sumstats convention so
/// callers can hand them directly to `score::run_staar_from_sumstats`. For
/// binary traits σ² is 1 and the scaling is a no-op.
///
/// For typical rare-variant genes (20 variants × MAC=5 × k=6 covariates),
/// total work is ~3,300 multiply-adds vs ~27,000,000 for the dense path.
pub fn score_gene_sparse(
    carriers: &[CarrierList],
    analysis: &AnalysisVectors,
) -> (Mat<f64>, Mat<f64>) {
    let m = carriers.len();
    let k = analysis.k;
    let is_binary = analysis.is_binary();

    // ── Phase 1: Single pass over carrier lists ──────────────────────
    // Compute U = G'r and GtX = G'X (or G'WX for binary) simultaneously.
    // Uses vcf_to_pheno remapping: skip samples without phenotype data.
    let mut u = vec![0.0f64; m];
    let mut gtx = vec![0.0f64; m * k];

    for (j, clist) in carriers.iter().enumerate() {
        let u_j = &mut u[j];
        let gtx_row = &mut gtx[j * k..(j + 1) * k];

        for &CarrierEntry { sample_idx, dosage } in &clist.entries {
            if dosage == 255 {
                continue; // missing
            }
            let pi = match analysis.vcf_to_pheno[sample_idx as usize] {
                Some(idx) => idx as usize,
                None => continue, // no phenotype — skip, don't multiply by zero
            };
            let d = dosage as f64;

            // Score vector: G'r (residuals already y - μ for binary)
            *u_j += d * analysis.residuals[pi];

            // G'X for continuous, G'WX for binary
            let x_row = &analysis.x_row_major[pi * k..(pi + 1) * k];
            if is_binary {
                let w_i = analysis.working_weights[pi];
                for c in 0..k {
                    gtx_row[c] += d * w_i * x_row[c];
                }
            } else {
                for c in 0..k {
                    gtx_row[c] += d * x_row[c];
                }
            }
        }
    }

    // ── Phase 2: G'G via inverted carrier index ──────────────────────
    // Build sample → variant map, accumulate G'G from each sample's contribution.
    // For binary: G'WG with working weights.
    //
    // Most samples carry 0 or 1 variant → zero off-diagonal contribution.
    // Only compound hets contribute off-diagonal terms.

    // pheno_idx → list of (variant_local_idx, dosage)
    let total_carriers: usize = carriers.iter().map(|c| c.len()).sum();
    let mut sample_variants: HashMap<u32, Vec<(u16, u8)>> =
        HashMap::with_capacity(total_carriers.min(1024));

    for (j, clist) in carriers.iter().enumerate() {
        for &CarrierEntry { sample_idx, dosage } in &clist.entries {
            if dosage == 255 {
                continue;
            }
            if let Some(pi) = analysis.vcf_to_pheno[sample_idx as usize] {
                sample_variants
                    .entry(pi)
                    .or_default()
                    .push((j as u16, dosage));
            }
        }
    }

    let mut gtg = vec![0.0f64; m * m];

    for (&pheno_idx, variant_dosages) in &sample_variants {
        let vd = variant_dosages.as_slice();
        let w = if is_binary {
            analysis.working_weights[pheno_idx as usize]
        } else {
            1.0
        };

        // Diagonal
        for &(vi, di) in vd {
            let i = vi as usize;
            gtg[i * m + i] += (di as f64) * (di as f64) * w;
        }

        // Off-diagonal (only when sample carries >1 variant)
        if vd.len() > 1 {
            for a in 0..vd.len() {
                for b in (a + 1)..vd.len() {
                    let (vi, di) = vd[a];
                    let (vj, dj) = vd[b];
                    let (i, j) = (vi as usize, vj as usize);
                    let product = (di as f64) * (dj as f64) * w;
                    gtg[i * m + j] += product;
                    gtg[j * m + i] += product;
                }
            }
        }
    }

    // ── Phase 3: Assemble kernel K = G'G − G'X (X'X)⁻¹ X'G ─────────
    // For binary: K = G'WG − G'WX (X'WX)⁻¹ X'WG
    // All small matrix operations: O(m × k² + m²)

    // tmp = (X'X)⁻¹ × (G'X)^T → [k × m]
    let mut tmp = vec![0.0f64; k * m];
    for i in 0..k {
        for j in 0..m {
            let mut sum = 0.0;
            for l in 0..k {
                sum += analysis.xtx_inv[i * k + l] * gtx[j * k + l];
            }
            tmp[i * m + j] = sum;
        }
    }

    // kernel[i,j] = gtg[i,j] − dot(gtx[i,:], tmp[:,j])
    let mut kernel_flat = gtg;
    for i in 0..m {
        for j in 0..m {
            let mut correction = 0.0;
            for l in 0..k {
                correction += gtx[i * k + l] * tmp[l * m + j];
            }
            kernel_flat[i * m + j] -= correction;
        }
    }

    // Scale to MetaSTAAR sumstats convention: U/σ², K/σ².
    let inv_s2 = 1.0 / analysis.sigma2;
    let u_mat = Mat::from_fn(m, 1, |i, _| u[i] * inv_s2);
    let k_mat = Mat::from_fn(m, m, |i, j| kernel_flat[i * m + j] * inv_s2);

    (u_mat, k_mat)
}

/// Run full STAAR analysis for a set of carrier lists using sparse scoring.
///
/// For SPA: builds dense G directly (needs raw genotypes for saddlepoint approx).
/// For non-SPA: computes U and K from carrier lists, then uses summary-stat tests.
pub fn run_staar_sparse(
    carriers: &[CarrierList],
    analysis: &AnalysisVectors,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
    use_spa: bool,
) -> StaarResult {
    let m = carriers.len();
    if m == 0 {
        return score::run_staar_from_sumstats(
            &Mat::zeros(0, 1),
            &Mat::zeros(0, 0),
            annotation_matrix,
            mafs,
            analysis.n_pheno,
        );
    }

    if use_spa && !analysis.fitted_values.is_empty() {
        let g = carriers_to_dense(carriers, analysis);
        score::run_staar(&g, annotation_matrix, mafs, &null_model_from_analysis(analysis), use_spa)
    } else {
        let (u, k) = score_gene_sparse(carriers, analysis);
        score::run_staar_from_sumstats(&u, &k, annotation_matrix, mafs, analysis.n_pheno)
    }
}

/// Slice pre-computed summary statistics (U vector and K matrix) for a subset
/// of variant indices. Used to avoid recomputing U/K for each mask within a gene.
///
/// Given full-gene U\[m×1\] and K\[m×m\], returns (U_sub\[s×1\], K_sub\[s×s\])
/// where s = indices.len(). All indices must be < m.
pub fn slice_sumstats(
    u: &Mat<f64>,
    k: &Mat<f64>,
    indices: &[usize],
) -> (Mat<f64>, Mat<f64>) {
    let s = indices.len();
    let u_sub = Mat::from_fn(s, 1, |i, _| u[(indices[i], 0)]);
    let k_sub = Mat::from_fn(s, s, |i, j| k[(indices[i], indices[j])]);
    (u_sub, k_sub)
}

/// Convert carrier lists to a dense faer::Mat for SPA / AI-STAAR paths.
/// Uses vcf_to_pheno remapping to build n_pheno × m matrix.
pub(crate) fn carriers_to_dense(carriers: &[CarrierList], analysis: &AnalysisVectors) -> Mat<f64> {
    let m = carriers.len();
    let mut g = Mat::zeros(analysis.n_pheno, m);
    for (j, clist) in carriers.iter().enumerate() {
        for &CarrierEntry { sample_idx, dosage } in &clist.entries {
            if dosage != 255 {
                if let Some(pi) = analysis.vcf_to_pheno[sample_idx as usize] {
                    g[(pi as usize, j)] = dosage as f64;
                }
            }
        }
    }
    g
}

/// Reconstruct a NullModel from AnalysisVectors for the SPA / AI-STAAR paths.
/// Uses compact n_pheno-sized arrays.
pub(crate) fn null_model_from_analysis(analysis: &AnalysisVectors) -> NullModel {
    let n = analysis.n_pheno;
    let k = analysis.k;

    let mut residuals = Mat::zeros(n, 1);
    for i in 0..n {
        residuals[(i, 0)] = analysis.residuals[i];
    }

    let mut x_matrix = Mat::zeros(n, k);
    for i in 0..n {
        for j in 0..k {
            x_matrix[(i, j)] = analysis.x_row_major[i * k + j];
        }
    }

    let mut xtx_inv = Mat::zeros(k, k);
    for i in 0..k {
        for j in 0..k {
            xtx_inv[(i, j)] = analysis.xtx_inv[i * k + j];
        }
    }

    NullModel {
        residuals,
        x_matrix,
        xtx_inv,
        sigma2: analysis.sigma2,
        n_samples: n,
        fitted_values: if analysis.fitted_values.is_empty() {
            None
        } else {
            Some(analysis.fitted_values.clone())
        },
        working_weights: if analysis.working_weights.is_empty() {
            None
        } else {
            Some(analysis.working_weights.clone())
        },
    }
}

/// Compute U/σ² from carrier lists (Phase 1 of score_gene_sparse).
/// O(total_MAC), no K. Used for genes exceeding MAX_K_VARIANTS.
pub fn compute_u_only(
    carriers: &[CarrierList],
    analysis: &AnalysisVectors,
) -> Vec<f64> {
    let m = carriers.len();
    let mut u = vec![0.0f64; m];
    for (j, clist) in carriers.iter().enumerate() {
        for &CarrierEntry { sample_idx, dosage } in &clist.entries {
            if dosage == 255 {
                continue;
            }
            if let Some(pi) = analysis.vcf_to_pheno[sample_idx as usize] {
                u[j] += dosage as f64 * analysis.residuals[pi as usize];
            }
        }
    }
    let inv_s2 = 1.0 / analysis.sigma2;
    for v in &mut u {
        *v *= inv_s2;
    }
    u
}

/// Slice K from a flat row-major array by local indices.
/// The flat array represents an [m × m] matrix; indices select a subset.
pub fn slice_sumstats_flat(
    k_flat: &[f64],
    m: usize,
    indices: &[usize],
) -> Mat<f64> {
    let s = indices.len();
    Mat::from_fn(s, s, |i, j| k_flat[indices[i] * m + indices[j]])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::staar::model;
    use faer::Mat;

    /// Build a simple test scenario with known values and verify
    /// sparse-vs-dense parity.
    #[test]
    fn sparse_dense_parity_continuous() {
        let n = 100; // samples
        let m = 5; // variants
        let k = 3; // covariates (intercept + 2)

        // Deterministic pseudo-random via LCG
        let mut rng_state = 42u64;
        let mut next_f64 = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        // Build covariate matrix with intercept
        let mut x = Mat::zeros(n, k);
        for i in 0..n {
            x[(i, 0)] = 1.0; // intercept
            for j in 1..k {
                x[(i, j)] = next_f64() * 2.0 - 1.0;
            }
        }

        // Build sparse genotype matrix (rare variants: ~5% carriers)
        let mut g_dense = Mat::zeros(n, m);
        let mut carrier_lists: Vec<Vec<CarrierEntry>> = vec![Vec::new(); m];
        for j in 0..m {
            for i in 0..n {
                if next_f64() < 0.05 {
                    let dosage = if next_f64() < 0.9 { 1 } else { 2 };
                    g_dense[(i, j)] = dosage as f64;
                    carrier_lists[j].push(CarrierEntry {
                        sample_idx: i as u32,
                        dosage: dosage as u8,
                    });
                }
            }
        }

        // Build phenotype
        let mut y = Mat::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = x[(i, 1)] * 0.5 + next_f64() * 0.1;
        }

        // Fit null model
        let null = model::fit_glm(&y, &x);

        let inv_s2 = 1.0 / null.sigma2;
        let u_dense = g_dense.transpose() * &null.residuals;
        let k_dense = null.compute_kernel(&g_dense);

        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; null.n_samples]);
        let carriers: Vec<CarrierList> = carrier_lists
            .into_iter()
            .map(|entries| CarrierList { entries })
            .collect();
        let (u_sparse, k_sparse) = score_gene_sparse(&carriers, &analysis);

        for j in 0..m {
            let expected = u_dense[(j, 0)] * inv_s2;
            let diff = (expected - u_sparse[(j, 0)]).abs();
            assert!(diff < 1e-10, "U[{j}]: expected={expected} got={} diff={diff}", u_sparse[(j, 0)]);
        }
        for i in 0..m {
            for j in 0..m {
                let expected = k_dense[(i, j)] * inv_s2;
                let diff = (expected - k_sparse[(i, j)]).abs();
                assert!(diff < 1e-10, "K[{i},{j}]: expected={expected} got={} diff={diff}", k_sparse[(i, j)]);
            }
        }
    }

    #[test]
    fn sparse_dense_parity_binary() {
        let n = 200;
        let m = 4;
        let k = 2;

        let mut rng_state = 123u64;
        let mut next_f64 = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut x = Mat::zeros(n, k);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            x[(i, 1)] = next_f64() * 2.0 - 1.0;
        }

        let mut g_dense = Mat::zeros(n, m);
        let mut carrier_lists: Vec<Vec<CarrierEntry>> = vec![Vec::new(); m];
        for j in 0..m {
            for i in 0..n {
                if next_f64() < 0.03 {
                    g_dense[(i, j)] = 1.0;
                    carrier_lists[j].push(CarrierEntry {
                        sample_idx: i as u32,
                        dosage: 1,
                    });
                }
            }
        }

        let mut y = Mat::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = if next_f64() < 0.3 { 1.0 } else { 0.0 };
        }

        let null = model::fit_logistic(&y, &x, 25);

        let u_dense = g_dense.transpose() * &null.residuals;
        let k_dense = null.compute_kernel(&g_dense);

        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; null.n_samples]);
        let carriers: Vec<CarrierList> = carrier_lists
            .into_iter()
            .map(|entries| CarrierList { entries })
            .collect();
        let (u_sparse, k_sparse) = score_gene_sparse(&carriers, &analysis);

        for j in 0..m {
            let diff = (u_dense[(j, 0)] - u_sparse[(j, 0)]).abs();
            assert!(
                diff < 1e-10,
                "Binary U[{j}] mismatch: dense={} sparse={} diff={diff}",
                u_dense[(j, 0)],
                u_sparse[(j, 0)]
            );
        }
        for i in 0..m {
            for j in 0..m {
                let diff = (k_dense[(i, j)] - k_sparse[(i, j)]).abs();
                assert!(
                    diff < 1e-10,
                    "Binary K[{i},{j}] mismatch: dense={} sparse={} diff={diff}",
                    k_dense[(i, j)],
                    k_sparse[(i, j)]
                );
            }
        }
    }

    /// Deterministic synthetic dataset shared by the e2e p-value parity tests.
    /// LCG keeps the build reproducible without pulling in an rng dependency.
    struct SyntheticContinuous {
        g: Mat<f64>,
        carriers: Vec<CarrierList>,
        y: Mat<f64>,
        x: Mat<f64>,
        mafs: Vec<f64>,
        ann: Vec<Vec<f64>>,
    }

    fn synthetic_continuous() -> SyntheticContinuous {
        let n: usize = 200;
        let m: usize = 6;
        let k: usize = 3;
        let mut state = 0xdeadbeefu64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((state >> 33) as f64) / ((1u64 << 31) as f64)
        };

        let mut x = Mat::zeros(n, k);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            x[(i, 1)] = next() * 2.0 - 1.0;
            x[(i, 2)] = next() * 2.0 - 1.0;
        }

        let mut g = Mat::zeros(n, m);
        let mut carriers: Vec<CarrierList> =
            (0..m).map(|_| CarrierList { entries: Vec::new() }).collect();
        let mut mac = vec![0u32; m];
        for j in 0..m {
            for i in 0..n {
                if next() < 0.05 {
                    let dosage: u8 = if next() < 0.9 { 1 } else { 2 };
                    g[(i, j)] = dosage as f64;
                    carriers[j].entries.push(CarrierEntry { sample_idx: i as u32, dosage });
                    mac[j] += dosage as u32;
                }
            }
            // Guarantee >=2 carriers so beta-weight masks aren't degenerate.
            if carriers[j].entries.len() < 2 {
                let i = (j * 7) % n;
                if g[(i, j)] == 0.0 {
                    g[(i, j)] = 1.0;
                    carriers[j].entries.push(CarrierEntry { sample_idx: i as u32, dosage: 1 });
                    mac[j] += 1;
                }
            }
        }
        let mafs: Vec<f64> = mac.iter().map(|&c| c as f64 / (2.0 * n as f64)).collect();

        let mut y = Mat::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = 0.7 * x[(i, 1)] - 0.3 * x[(i, 2)] + (next() - 0.5) * 0.4;
        }

        let ann: Vec<Vec<f64>> = (0..2).map(|_| (0..m).map(|_| next()).collect()).collect();
        SyntheticContinuous { g, carriers, y, x, mafs, ann }
    }

    fn assert_staar_eq(actual: &StaarResult, expected: &StaarResult) {
        const TOL: f64 = 1e-9;
        const TOL_SKAT: f64 = 1e-6;
        let pairs = [
            ("burden_1_25", actual.burden_1_25, expected.burden_1_25, TOL),
            ("burden_1_1", actual.burden_1_1, expected.burden_1_1, TOL),
            ("skat_1_25", actual.skat_1_25, expected.skat_1_25, TOL_SKAT),
            ("skat_1_1", actual.skat_1_1, expected.skat_1_1, TOL_SKAT),
            ("acat_v_1_25", actual.acat_v_1_25, expected.acat_v_1_25, TOL),
            ("acat_v_1_1", actual.acat_v_1_1, expected.acat_v_1_1, TOL),
            ("acat_o", actual.acat_o, expected.acat_o, TOL),
            ("staar_o", actual.staar_o, expected.staar_o, TOL_SKAT),
            ("staar_b_1_25", actual.staar_b_1_25, expected.staar_b_1_25, TOL),
            ("staar_b_1_1", actual.staar_b_1_1, expected.staar_b_1_1, TOL),
            ("staar_s_1_25", actual.staar_s_1_25, expected.staar_s_1_25, TOL_SKAT),
            ("staar_s_1_1", actual.staar_s_1_1, expected.staar_s_1_1, TOL_SKAT),
            ("staar_a_1_25", actual.staar_a_1_25, expected.staar_a_1_25, TOL),
            ("staar_a_1_1", actual.staar_a_1_1, expected.staar_a_1_1, TOL),
        ];
        for (name, a, e, tol) in pairs {
            assert!((a - e).abs() < tol, "{name}: actual={a:e} expected={e:e}");
        }
        assert_eq!(actual.per_annotation.len(), expected.per_annotation.len());
        for (ci, (a_ann, e_ann)) in actual.per_annotation.iter().zip(&expected.per_annotation).enumerate() {
            for (ti, (&a, &e)) in a_ann.iter().zip(e_ann).enumerate() {
                let tol = if ti == 2 || ti == 3 { TOL_SKAT } else { TOL };
                assert!((a - e).abs() < tol, "per_annotation[{ci}][{ti}]: actual={a:e} expected={e:e}");
            }
        }
    }

    /// Build the U/K oracle by hand: U=G'r and K=G'(I-H)G scaled by 1/σ², then
    /// fed to run_staar_from_sumstats. This entry point is the path validated
    /// against R STAAR by `staar_continuous_matches_r`.
    fn oracle(
        g: &Mat<f64>, null: &crate::staar::model::NullModel,
        ann: &[Vec<f64>], mafs: &[f64],
    ) -> StaarResult {
        let inv_s2 = 1.0 / null.sigma2;
        let m = mafs.len();
        let u = g.transpose() * &null.residuals;
        let k = null.compute_kernel(g);
        let u_scaled = Mat::from_fn(m, 1, |i, _| u[(i, 0)] * inv_s2);
        let k_scaled = Mat::from_fn(m, m, |i, j| k[(i, j)] * inv_s2);
        score::run_staar_from_sumstats(&u_scaled, &k_scaled, ann, mafs, null.n_samples)
    }

    /// Raw genotype path (`score::run_staar`) must produce identical p-values
    /// to the U/K oracle for every Burden/SKAT/ACAT-V/omnibus statistic.
    #[test]
    fn raw_path_matches_oracle() {
        let d = synthetic_continuous();
        let null = model::fit_glm(&d.y, &d.x);
        let raw = score::run_staar(&d.g, &d.ann, &d.mafs, &null, false);
        let expected = oracle(&d.g, &null, &d.ann, &d.mafs);
        assert_staar_eq(&raw, &expected);
    }

    /// Sparse carrier-list path (`run_staar_sparse`) must produce identical
    /// p-values to the U/K oracle. Extends sparse-dense parity from U/K to
    /// the full set of test statistics, including the STAAR-O omnibus.
    #[test]
    fn sparse_path_matches_oracle() {
        let d = synthetic_continuous();
        let null = model::fit_glm(&d.y, &d.x);
        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; null.n_samples]);
        let sparse = run_staar_sparse(&d.carriers, &analysis, &d.ann, &d.mafs, false);
        let expected = oracle(&d.g, &null, &d.ann, &d.mafs);
        assert_staar_eq(&sparse, &expected);
    }

    #[test]
    fn singleton_variant() {
        let n = 50;
        let k = 2;

        let mut x = Mat::zeros(n, k);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            x[(i, 1)] = i as f64 / n as f64;
        }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = (i as f64 / n as f64) * 0.5;
        }

        let null = model::fit_glm(&y, &x);
        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; null.n_samples]);

        // Single variant with MAC=1
        let carriers = vec![
            CarrierList {
                entries: vec![CarrierEntry { sample_idx: 10, dosage: 1 }],
            },
            CarrierList {
                entries: vec![CarrierEntry { sample_idx: 20, dosage: 1 }],
            },
        ];

        let (u, k) = score_gene_sparse(&carriers, &analysis);
        assert_eq!(u.nrows(), 2);
        assert_eq!(k.nrows(), 2);
        assert_eq!(k.ncols(), 2);
        let expected = analysis.residuals[10] / analysis.sigma2;
        assert!((u[(0, 0)] - expected).abs() < 1e-10);
    }

    #[test]
    fn empty_gene() {
        let analysis = AnalysisVectors {
            residuals: vec![0.0; 10],
            x_row_major: vec![1.0; 20],
            xtx_inv: vec![1.0; 4],
            sigma2: 1.0,
            k: 2,
            n_pheno: 10,
            n_vcf_total: 10,
            working_weights: Vec::new(),
            fitted_values: Vec::new(),
            vcf_to_pheno: (0..10).map(|i| Some(i as u32)).collect(),
        };

        let carriers: Vec<CarrierList> = Vec::new();

        let (u, k) = score_gene_sparse(&carriers, &analysis);
        assert_eq!(u.nrows(), 0);
        assert_eq!(k.nrows(), 0);
    }

    /// Full mask: all variants selected. U/K from slice must equal full U/K.
    #[test]
    fn full_mask_equals_unsliced() {
        let n = 50;
        let m = 4;
        let k = 2;
        let mut rng_state = 99u64;
        let mut next_f64 = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut x = Mat::zeros(n, k);
        for i in 0..n { x[(i, 0)] = 1.0; for j in 1..k { x[(i, j)] = next_f64(); } }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n { y[(i, 0)] = next_f64(); }

        let mut carrier_lists: Vec<Vec<CarrierEntry>> = vec![Vec::new(); m];
        for cl in &mut carrier_lists {
            for i in 0..n {
                if next_f64() < 0.08 {
                    cl.push(CarrierEntry { sample_idx: i as u32, dosage: 1 });
                }
            }
        }

        let null = model::fit_glm(&y, &x);
        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; n]);
        let carriers: Vec<CarrierList> = carrier_lists.into_iter()
            .map(|entries| CarrierList { entries }).collect();
        let (u_full, k_full) = score_gene_sparse(&carriers, &analysis);

        // Full mask: all indices
        let all_indices: Vec<usize> = (0..m).collect();
        let (u_sliced, k_sliced) = slice_sumstats(&u_full, &k_full, &all_indices);

        for i in 0..m {
            assert!((u_full[(i, 0)] - u_sliced[(i, 0)]).abs() < 1e-14,
                "Full mask U mismatch at {i}");
        }
        for i in 0..m {
            for j in 0..m {
                assert!((k_full[(i, j)] - k_sliced[(i, j)]).abs() < 1e-14,
                    "Full mask K mismatch at ({i},{j})");
            }
        }
    }

    /// Random subset: slice U/K, then compare against recomputed U/K
    /// from the carrier subset only.
    #[test]
    fn random_subset_slice_equals_recompute() {
        let n = 80;
        let m = 6;
        let k = 2;
        let mut rng_state = 777u64;
        let mut next_f64 = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut x = Mat::zeros(n, k);
        for i in 0..n { x[(i, 0)] = 1.0; for j in 1..k { x[(i, j)] = next_f64(); } }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n { y[(i, 0)] = next_f64(); }

        let mut carrier_lists: Vec<Vec<CarrierEntry>> = vec![Vec::new(); m];
        for cl in &mut carrier_lists {
            for i in 0..n {
                if next_f64() < 0.06 {
                    cl.push(CarrierEntry { sample_idx: i as u32, dosage: 1 });
                }
            }
        }

        let null = model::fit_glm(&y, &x);
        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; n]);
        let carriers: Vec<CarrierList> = carrier_lists.into_iter()
            .map(|entries| CarrierList { entries }).collect();

        // Full gene U/K
        let (u_full, k_full) = score_gene_sparse(&carriers, &analysis);

        // Random subset: variants [1, 3, 5]
        let subset = vec![1usize, 3, 5];
        let (u_sliced, k_sliced) = slice_sumstats(&u_full, &k_full, &subset);

        // Recompute from carrier subset
        let subset_carriers: Vec<CarrierList> = subset.iter()
            .map(|&i| carriers[i].clone()).collect();
        let (u_recomputed, k_recomputed) = score_gene_sparse(&subset_carriers, &analysis);

        let s = subset.len();
        for i in 0..s {
            let diff = (u_sliced[(i, 0)] - u_recomputed[(i, 0)]).abs();
            assert!(diff < 1e-12, "Subset U[{i}] diff={diff}: sliced={} recomputed={}",
                u_sliced[(i, 0)], u_recomputed[(i, 0)]);
        }
        for i in 0..s {
            for j in 0..s {
                let diff = (k_sliced[(i, j)] - k_recomputed[(i, j)]).abs();
                assert!(diff < 1e-10, "Subset K[{i},{j}] diff={diff}: sliced={} recomputed={}",
                    k_sliced[(i, j)], k_recomputed[(i, j)]);
            }
        }
    }

    /// Dense carrier variant: a variant with many carriers (50% carrier rate).
    /// Tests accumulation accuracy with large sums.
    #[test]
    fn dense_carrier_variant() {
        let n = 200;
        let k = 2;
        let mut rng_state = 555u64;
        let mut next_f64 = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut x = Mat::zeros(n, k);
        for i in 0..n { x[(i, 0)] = 1.0; for j in 1..k { x[(i, j)] = next_f64(); } }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n { y[(i, 0)] = next_f64(); }

        // One dense variant (50% carriers) + one sparse
        let mut g_dense = Mat::zeros(n, 2);
        let mut carrier_lists: Vec<Vec<CarrierEntry>> = vec![Vec::new(); 2];
        for i in 0..n {
            if next_f64() < 0.5 { // dense
                let d = if next_f64() < 0.8 { 1 } else { 2 };
                g_dense[(i, 0)] = d as f64;
                carrier_lists[0].push(CarrierEntry { sample_idx: i as u32, dosage: d as u8 });
            }
            if next_f64() < 0.02 { // sparse
                g_dense[(i, 1)] = 1.0;
                carrier_lists[1].push(CarrierEntry { sample_idx: i as u32, dosage: 1 });
            }
        }

        let null = model::fit_glm(&y, &x);
        let inv_s2 = 1.0 / null.sigma2;
        let u_dense = g_dense.transpose() * &null.residuals;
        let k_dense = null.compute_kernel(&g_dense);

        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; n]);
        let carriers: Vec<CarrierList> = carrier_lists.into_iter()
            .map(|entries| CarrierList { entries }).collect();
        let (u_sparse, k_sparse) = score_gene_sparse(&carriers, &analysis);

        for j in 0..2 {
            let diff = (u_dense[(j, 0)] * inv_s2 - u_sparse[(j, 0)]).abs();
            assert!(diff < 1e-10, "Dense carrier U[{j}] diff={diff}");
        }
        for i in 0..2 {
            for j in 0..2 {
                let diff = (k_dense[(i, j)] * inv_s2 - k_sparse[(i, j)]).abs();
                assert!(diff < 1e-10, "Dense carrier K[{i},{j}] diff={diff}");
            }
        }
    }

    /// Permutation invariance: scoring the same carriers in a different order
    /// must produce identical U and K (after accounting for the permutation).
    #[test]
    fn permutation_invariance() {
        let n = 60;
        let m = 4;
        let k = 2;
        let mut rng_state = 321u64;
        let mut next_f64 = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut x = Mat::zeros(n, k);
        for i in 0..n { x[(i, 0)] = 1.0; for j in 1..k { x[(i, j)] = next_f64(); } }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n { y[(i, 0)] = next_f64(); }

        let mut carrier_lists: Vec<Vec<CarrierEntry>> = vec![Vec::new(); m];
        for cl in &mut carrier_lists {
            for i in 0..n {
                if next_f64() < 0.07 {
                    cl.push(CarrierEntry { sample_idx: i as u32, dosage: 1 });
                }
            }
        }

        let null = model::fit_glm(&y, &x);
        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; n]);
        let carriers: Vec<CarrierList> = carrier_lists.into_iter()
            .map(|entries| CarrierList { entries }).collect();

        // Original order: [0, 1, 2, 3]
        let (u_orig, k_orig) = score_gene_sparse(&carriers, &analysis);

        // Permuted order: [2, 0, 3, 1]
        let perm = [2usize, 0, 3, 1];
        let permuted: Vec<CarrierList> = perm.iter().map(|&i| carriers[i].clone()).collect();
        let (u_perm, k_perm) = score_gene_sparse(&permuted, &analysis);

        // U and K must match under the same permutation
        for (new_j, &orig_j) in perm.iter().enumerate() {
            let diff = (u_orig[(orig_j, 0)] - u_perm[(new_j, 0)]).abs();
            assert!(diff < 1e-14, "Permutation U mismatch: orig[{orig_j}]={} perm[{new_j}]={}",
                u_orig[(orig_j, 0)], u_perm[(new_j, 0)]);
        }
        for (ni, &oi) in perm.iter().enumerate() {
            for (nj, &oj) in perm.iter().enumerate() {
                let diff = (k_orig[(oi, oj)] - k_perm[(ni, nj)]).abs();
                assert!(diff < 1e-12,
                    "Permutation K mismatch: orig[{oi},{oj}]={} perm[{ni},{nj}]={}",
                    k_orig[(oi, oj)], k_perm[(ni, nj)]);
            }
        }
    }

    /// Duplicate variant idempotence: [v, v, v] must not corrupt results.
    /// Scoring with duplicate carrier lists should produce results where the
    /// duplicated entries are mathematically consistent (each duplicate gets
    /// the same U value, K is rank-deficient but not garbage).
    #[test]
    fn duplicate_variant_idempotence() {
        let n = 50;
        let k = 2;
        let mut rng_state = 654u64;
        let mut next_f64 = || -> f64 {
            rng_state = rng_state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (rng_state >> 33) as f64 / (1u64 << 31) as f64
        };

        let mut x = Mat::zeros(n, k);
        for i in 0..n { x[(i, 0)] = 1.0; for j in 1..k { x[(i, j)] = next_f64(); } }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n { y[(i, 0)] = next_f64(); }

        let mut entries = Vec::new();
        for i in 0..n {
            if next_f64() < 0.1 {
                entries.push(CarrierEntry { sample_idx: i as u32, dosage: 1 });
            }
        }

        let null = model::fit_glm(&y, &x);
        let analysis = AnalysisVectors::from_null_model(&null, &vec![true; n]);

        let cl = CarrierList { entries };

        // Single variant
        let (u_single, _k_single) = score_gene_sparse(std::slice::from_ref(&cl), &analysis);

        // Triplicated: [v, v, v]
        let (u_triple, k_triple) = score_gene_sparse(
            &[cl.clone(), cl.clone(), cl.clone()], &analysis,
        );

        // All three U values must be identical to the single
        for j in 0..3 {
            let diff = (u_single[(0, 0)] - u_triple[(j, 0)]).abs();
            assert!(diff < 1e-14,
                "Duplicate U[{j}] should equal single: single={} triple={}",
                u_single[(0, 0)], u_triple[(j, 0)]);
        }

        // K must be symmetric and all entries identical (since all columns of G are the same)
        for i in 0..3 {
            for j in 0..3 {
                let diff = (k_triple[(0, 0)] - k_triple[(i, j)]).abs();
                assert!(diff < 1e-12,
                    "Duplicate K[{i},{j}] should equal K[0,0]: {} vs {}",
                    k_triple[(i, j)], k_triple[(0, 0)]);
            }
        }
    }
}
