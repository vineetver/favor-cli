//! AI-STAAR. Mirrors `STAAR/R/AI_STAAR.R` and the weight generator in
//! `STAARpipeline/R/staar2aistaar_nullmodel.R`.

#![allow(clippy::needless_range_loop)]

use faer::Mat;

use super::carrier::sparse_score::{self, AnalysisVectors};
use super::model::NullModel;
use super::score::{self, StaarResult};
use super::stats;
use crate::store::cohort::variants::{CarrierEntry, CarrierList};

#[derive(Clone)]
pub struct AncestryInfo {
    pub group: Vec<usize>,
    pub n_pops: usize,
    /// n_pops × (B+1). Column 0 is all-ones, columns 1..=B are |N(0,1)|.
    pub pop_weights_1_1: Vec<Vec<f64>>,
    pub pop_weights_1_25: Vec<Vec<f64>>,
}

impl AncestryInfo {
    pub fn new(group: Vec<usize>, n_pops: usize, n_base_tests: usize, seed: u64) -> Self {
        let n_cols = n_base_tests + 1;
        let mut pop_weights_1_1 = vec![vec![1.0; n_cols]; n_pops];
        let mut pop_weights_1_25 = vec![vec![1.0; n_cols]; n_pops];
        let mut rng = Lcg::new(seed);
        for k in 0..n_pops {
            for b in 1..n_cols {
                pop_weights_1_1[k][b] = rng.next_normal().abs();
            }
        }
        for k in 0..n_pops {
            for b in 1..n_cols {
                pop_weights_1_25[k][b] = rng.next_normal().abs();
            }
        }
        Self {
            group,
            n_pops,
            pop_weights_1_1,
            pop_weights_1_25,
        }
    }
}

/// `pop_mafs[v][p] = MAC_p / (2 * n_p)`.
pub fn pop_mafs_from_carriers(
    carriers: &[CarrierList],
    analysis: &AnalysisVectors,
    ancestry: &AncestryInfo,
) -> Vec<Vec<f64>> {
    let m = carriers.len();
    let n_pops = ancestry.n_pops;

    let mut pop_n = vec![0usize; n_pops];
    for &g in &ancestry.group {
        pop_n[g] += 1;
    }

    let mut pop_mac = vec![vec![0u64; n_pops]; m];
    for (j, clist) in carriers.iter().enumerate() {
        for &CarrierEntry { sample_idx, dosage } in &clist.entries {
            if dosage == 255 {
                continue;
            }
            let pi = match analysis.vcf_to_pheno[sample_idx as usize] {
                Some(idx) => idx as usize,
                None => continue,
            };
            pop_mac[j][ancestry.group[pi]] += dosage as u64;
        }
    }

    pop_mac
        .into_iter()
        .map(|mac_per_pop| {
            mac_per_pop
                .into_iter()
                .enumerate()
                .map(|(p, mac)| {
                    if pop_n[p] == 0 {
                        0.0
                    } else {
                        mac as f64 / (2.0 * pop_n[p] as f64)
                    }
                })
                .collect()
        })
        .collect()
}

pub fn run_ai_staar_gene(
    carriers: &[CarrierList],
    analysis: &AnalysisVectors,
    ancestry: &AncestryInfo,
    annotation_matrix: &[Vec<f64>],
    use_spa: bool,
) -> StaarResult {
    let g = sparse_score::carriers_to_dense(carriers, analysis);
    let pop_mafs = pop_mafs_from_carriers(carriers, analysis, ancestry);
    let null = sparse_score::null_model_from_analysis(analysis);
    run_ai_staar(&g, annotation_matrix, &pop_mafs, ancestry, &null, use_spa)
}

/// `pop_mafs[v][p]` = unfolded MAF of variant v in population p.
pub fn run_ai_staar(
    g: &Mat<f64>,
    annotation_matrix: &[Vec<f64>],
    pop_mafs: &[Vec<f64>],
    ancestry: &AncestryInfo,
    null: &NullModel,
    use_spa: bool,
) -> StaarResult {
    let m = g.ncols();
    let n = g.nrows();
    let n_pops = ancestry.n_pops;
    let n_channels = annotation_matrix.len();

    if m == 0 || n_pops == 0 || ancestry.pop_weights_1_1.is_empty() {
        let pooled = pooled_mafs_from_g(g);
        return score::run_staar(g, annotation_matrix, &pooled, null, use_spa);
    }

    // a_p[k] = dbeta(mean folded MAF in pop k, 1, 25), else 0.
    let mut a_p = vec![0.0f64; n_pops];
    for k in 0..n_pops {
        let mut sum = 0.0;
        for v in 0..m {
            let maf = pop_mafs.get(v).and_then(|r| r.get(k)).copied().unwrap_or(0.0);
            sum += maf.min(1.0 - maf);
        }
        let mean = sum / m as f64;
        a_p[k] = if mean > 0.0 {
            score::beta_density_weight(mean, 1.0, 25.0)
        } else {
            0.0
        };
    }

    let mafs = pooled_mafs_from_g(g);

    let n_cols = ancestry.pop_weights_1_1[0].len();
    let mut all_runs: Vec<Vec<f64>> = Vec::with_capacity(2 * n_cols);
    let mut g1 = Mat::zeros(n, m);
    let mut g2 = Mat::zeros(n, m);

    for b in 0..n_cols {
        let w_b_1: Vec<f64> = (0..n_pops).map(|k| ancestry.pop_weights_1_1[k][b]).collect();
        let w_b_2: Vec<f64> = (0..n_pops)
            .map(|k| a_p[k] * ancestry.pop_weights_1_25[k][b])
            .collect();

        for i in 0..n {
            let pop = ancestry.group[i];
            let s1 = w_b_1[pop];
            let s2 = w_b_2[pop];
            for j in 0..m {
                let v = g[(i, j)];
                g1[(i, j)] = s1 * v;
                g2[(i, j)] = s2 * v;
            }
        }

        let r1 = score::run_staar(&g1, annotation_matrix, &mafs, null, use_spa);
        let r2 = score::run_staar(&g2, annotation_matrix, &mafs, null, use_spa);
        all_runs.push(flatten_staar(&r1, n_channels));
        all_runs.push(flatten_staar(&r2, n_channels));
    }

    let len = 6 * (1 + n_channels);
    let mut aggregated = vec![0.0f64; len];
    let mut row = Vec::with_capacity(all_runs.len());
    for i in 0..len {
        row.clear();
        row.extend(all_runs.iter().map(|v| v[i]));
        aggregated[i] = stats::cauchy_combine(&row);
    }

    unflatten_staar(&aggregated, n_channels)
}

fn pooled_mafs_from_g(g: &Mat<f64>) -> Vec<f64> {
    let n = g.nrows();
    let m = g.ncols();
    (0..m)
        .map(|j| {
            let mut mac = 0.0;
            for i in 0..n {
                mac += g[(i, j)];
            }
            (mac / (2.0 * n as f64)).clamp(1e-10, 0.499)
        })
        .collect()
}

// Block order matches AI_STAAR.R indexing into pvalues_aggregate:
// SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25), ACAT-V(1,1).
// Each block: MAF-only p-value first, then one entry per annotation channel.
// Length = 6 * (1 + n_channels).
fn flatten_staar(result: &StaarResult, n_channels: usize) -> Vec<f64> {
    let na = 1 + n_channels;
    let mut v = Vec::with_capacity(6 * na);
    // per_annotation index: 0=B(1,25) 1=B(1,1) 2=S(1,25) 3=S(1,1) 4=A(1,25) 5=A(1,1)
    v.push(result.skat_1_25);
    for ch in 0..n_channels {
        v.push(result.per_annotation[ch][2]);
    }
    v.push(result.skat_1_1);
    for ch in 0..n_channels {
        v.push(result.per_annotation[ch][3]);
    }
    v.push(result.burden_1_25);
    for ch in 0..n_channels {
        v.push(result.per_annotation[ch][0]);
    }
    v.push(result.burden_1_1);
    for ch in 0..n_channels {
        v.push(result.per_annotation[ch][1]);
    }
    v.push(result.acat_v_1_25);
    for ch in 0..n_channels {
        v.push(result.per_annotation[ch][4]);
    }
    v.push(result.acat_v_1_1);
    for ch in 0..n_channels {
        v.push(result.per_annotation[ch][5]);
    }
    v
}

// Inverse of flatten_staar. STAAR-O = CCT over the full vector;
// ACAT-O = CCT over the 6 base p-values; per-test omni = CCT over each block.
fn unflatten_staar(flat: &[f64], n_channels: usize) -> StaarResult {
    let na = 1 + n_channels;

    let skat_1_25 = flat[0];
    let skat_1_1 = flat[na];
    let burden_1_25 = flat[2 * na];
    let burden_1_1 = flat[3 * na];
    let acat_v_1_25 = flat[4 * na];
    let acat_v_1_1 = flat[5 * na];

    let mut per_annotation: Vec<[f64; 6]> = Vec::with_capacity(n_channels);
    for ch in 0..n_channels {
        per_annotation.push([
            flat[2 * na + 1 + ch],
            flat[3 * na + 1 + ch],
            flat[1 + ch],
            flat[na + 1 + ch],
            flat[4 * na + 1 + ch],
            flat[5 * na + 1 + ch],
        ]);
    }

    let staar_s_1_25 = stats::cauchy_combine(&flat[0..na]);
    let staar_s_1_1 = stats::cauchy_combine(&flat[na..2 * na]);
    let staar_b_1_25 = stats::cauchy_combine(&flat[2 * na..3 * na]);
    let staar_b_1_1 = stats::cauchy_combine(&flat[3 * na..4 * na]);
    let staar_a_1_25 = stats::cauchy_combine(&flat[4 * na..5 * na]);
    let staar_a_1_1 = stats::cauchy_combine(&flat[5 * na..6 * na]);

    let acat_o = stats::cauchy_combine(&[
        skat_1_25,
        skat_1_1,
        burden_1_25,
        burden_1_1,
        acat_v_1_25,
        acat_v_1_1,
    ]);
    let staar_o = stats::cauchy_combine(flat);

    StaarResult {
        burden_1_25,
        burden_1_1,
        skat_1_25,
        skat_1_1,
        acat_v_1_25,
        acat_v_1_1,
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

// LCG + Box-Muller for ensemble-weight generation.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.wrapping_add(0x9E3779B97F4A7C15),
        }
    }
    fn next_u64(&mut self) -> u64 {
        self.state = self
            .state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }
    fn next_unit(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) / ((1u64 << 53) as f64)
    }
    fn next_normal(&mut self) -> f64 {
        let u1 = self.next_unit().max(f64::MIN_POSITIVE);
        let u2 = self.next_unit();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

#[cfg(test)]
mod tests {
    use super::super::model;
    use super::*;

    #[test]
    fn pop_mafs_match_per_pop_mac() {
        let analysis = AnalysisVectors {
            residuals: vec![0.0; 6],
            x_row_major: vec![1.0; 6],
            xtx_inv: vec![1.0],
            sigma2: 1.0,
            k: 1,
            n_pheno: 6,
            n_vcf_total: 6,
            working_weights: Vec::new(),
            fitted_values: Vec::new(),
            vcf_to_pheno: (0..6).map(|i| Some(i as u32)).collect(),
            kinship: None,
        };
        let ancestry = AncestryInfo::new(vec![0, 0, 0, 1, 1, 1], 2, 0, 7590);
        let carriers = vec![
            CarrierList {
                entries: vec![
                    CarrierEntry { sample_idx: 0, dosage: 1 },
                    CarrierEntry { sample_idx: 2, dosage: 2 },
                ],
            },
            CarrierList {
                entries: vec![
                    CarrierEntry { sample_idx: 1, dosage: 1 },
                    CarrierEntry { sample_idx: 4, dosage: 1 },
                ],
            },
        ];

        let pop_mafs = pop_mafs_from_carriers(&carriers, &analysis, &ancestry);

        // variant 0: pop0 MAC=3 (1+2), pop1 MAC=0 → 3/6, 0/6
        assert!((pop_mafs[0][0] - 0.5).abs() < 1e-12);
        assert!(pop_mafs[0][1].abs() < 1e-12);
        // variant 1: pop0 MAC=1, pop1 MAC=1 → 1/6, 1/6
        assert!((pop_mafs[1][0] - 1.0 / 6.0).abs() < 1e-12);
        assert!((pop_mafs[1][1] - 1.0 / 6.0).abs() < 1e-12);
    }

    #[test]
    fn pop_mafs_skip_samples_without_phenotype() {
        let analysis = AnalysisVectors {
            residuals: vec![0.0; 3],
            x_row_major: vec![1.0; 3],
            xtx_inv: vec![1.0],
            sigma2: 1.0,
            k: 1,
            n_pheno: 3,
            n_vcf_total: 4,
            working_weights: Vec::new(),
            fitted_values: Vec::new(),
            vcf_to_pheno: vec![Some(0), None, Some(1), Some(2)],
            kinship: None,
        };
        let ancestry = AncestryInfo::new(vec![0, 0, 1], 2, 0, 7590);
        let carriers = vec![CarrierList {
            entries: vec![
                CarrierEntry { sample_idx: 0, dosage: 1 },
                CarrierEntry { sample_idx: 1, dosage: 2 }, // dropped: no phenotype
                CarrierEntry { sample_idx: 3, dosage: 1 },
            ],
        }];
        let pop_mafs = pop_mafs_from_carriers(&carriers, &analysis, &ancestry);
        // pop0: 2 samples, MAC=1 → 1/4. pop1: 1 sample, MAC=1 → 1/2.
        assert!((pop_mafs[0][0] - 0.25).abs() < 1e-12);
        assert!((pop_mafs[0][1] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn ensemble_weights_first_column_is_ones() {
        let info = AncestryInfo::new(vec![0, 0, 1, 1], 2, 5, 7590);
        assert_eq!(info.pop_weights_1_1.len(), 2);
        assert_eq!(info.pop_weights_1_25.len(), 2);
        assert_eq!(info.pop_weights_1_1[0].len(), 6);
        for k in 0..2 {
            assert_eq!(info.pop_weights_1_1[k][0], 1.0);
            assert_eq!(info.pop_weights_1_25[k][0], 1.0);
        }
        for k in 0..2 {
            for b in 1..6 {
                assert!(info.pop_weights_1_1[k][b] >= 0.0);
                assert!(info.pop_weights_1_25[k][b] >= 0.0);
            }
        }
    }

    #[test]
    fn ensemble_weights_are_deterministic() {
        let a = AncestryInfo::new(vec![0, 1], 2, 4, 7590);
        let b = AncestryInfo::new(vec![0, 1], 2, 4, 7590);
        assert_eq!(a.pop_weights_1_1, b.pop_weights_1_1);
        assert_eq!(a.pop_weights_1_25, b.pop_weights_1_25);
    }

    // B=0 → two STAAR runs (G, G row-scaled by a_p) Cauchy-combined element-wise.
    #[test]
    fn run_ai_staar_b0_matches_two_run_cauchy() {
        let n = 30;
        let m = 4;
        let n_pops = 2;
        let n_channels = 3;

        let mut state = 0xfeedfaceu64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((state >> 33) as f64) / ((1u64 << 31) as f64)
        };

        let mut x = Mat::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            x[(i, 1)] = next() * 2.0 - 1.0;
        }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = 0.3 * x[(i, 1)] + (next() - 0.5) * 0.4;
        }
        let null = model::fit_glm(&y, &x);

        let mut g = Mat::zeros(n, m);
        for j in 0..m {
            for i in 0..n {
                if next() < 0.15 {
                    g[(i, j)] = 1.0;
                }
            }
        }
        let group: Vec<usize> = (0..n).map(|i| if i < n / 2 { 0 } else { 1 }).collect();
        let ancestry = AncestryInfo::new(group.clone(), n_pops, 0, 7590);

        let ann: Vec<Vec<f64>> = (0..n_channels)
            .map(|c| (0..m).map(|j| 0.2 + 0.1 * (c + j) as f64).collect())
            .collect();

        let mut pop_n = vec![0usize; n_pops];
        for &p in &group {
            pop_n[p] += 1;
        }
        let mut pop_mafs: Vec<Vec<f64>> = vec![vec![0.0; n_pops]; m];
        for j in 0..m {
            for i in 0..n {
                pop_mafs[j][group[i]] += g[(i, j)];
            }
            for k in 0..n_pops {
                pop_mafs[j][k] /= 2.0 * pop_n[k] as f64;
            }
        }

        let actual = run_ai_staar(&g, &ann, &pop_mafs, &ancestry, &null, false);

        let mut a_p = vec![0.0f64; n_pops];
        for k in 0..n_pops {
            let mut sum = 0.0;
            for j in 0..m {
                let maf = pop_mafs[j][k];
                sum += maf.min(1.0 - maf);
            }
            let mean = sum / m as f64;
            a_p[k] = if mean > 0.0 {
                score::beta_density_weight(mean, 1.0, 25.0)
            } else {
                0.0
            };
        }
        let mafs: Vec<f64> = (0..m)
            .map(|j| {
                let mut mac = 0.0;
                for i in 0..n {
                    mac += g[(i, j)];
                }
                (mac / (2.0 * n as f64)).clamp(1e-10, 0.499)
            })
            .collect();

        let r1 = score::run_staar(&g, &ann, &mafs, &null, false);
        let mut g2 = Mat::zeros(n, m);
        for i in 0..n {
            let s = a_p[group[i]];
            for j in 0..m {
                g2[(i, j)] = s * g[(i, j)];
            }
        }
        let r2 = score::run_staar(&g2, &ann, &mafs, &null, false);

        let f1 = flatten_staar(&r1, n_channels);
        let f2 = flatten_staar(&r2, n_channels);
        assert_eq!(f1.len(), 6 * (1 + n_channels));
        let mut aggregated = vec![0.0; f1.len()];
        for i in 0..f1.len() {
            aggregated[i] = stats::cauchy_combine(&[f1[i], f2[i]]);
        }
        let expected = unflatten_staar(&aggregated, n_channels);

        assert!((actual.staar_o - expected.staar_o).abs() < 1e-12);
        assert!((actual.acat_o - expected.acat_o).abs() < 1e-12);
        assert!((actual.staar_b_1_25 - expected.staar_b_1_25).abs() < 1e-12);
        assert!((actual.staar_b_1_1 - expected.staar_b_1_1).abs() < 1e-12);
        assert!((actual.staar_s_1_25 - expected.staar_s_1_25).abs() < 1e-12);
        assert!((actual.staar_s_1_1 - expected.staar_s_1_1).abs() < 1e-12);
        assert!((actual.staar_a_1_25 - expected.staar_a_1_25).abs() < 1e-12);
        assert!((actual.staar_a_1_1 - expected.staar_a_1_1).abs() < 1e-12);
        assert!((actual.burden_1_25 - expected.burden_1_25).abs() < 1e-12);
        assert!((actual.skat_1_1 - expected.skat_1_1).abs() < 1e-12);
        assert_eq!(actual.per_annotation.len(), n_channels);
        for ch in 0..n_channels {
            for t in 0..6 {
                assert!((actual.per_annotation[ch][t] - expected.per_annotation[ch][t]).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn run_ai_staar_gene_matches_dense() {
        let n = 24;
        let m = 3;
        let n_pops = 2;
        let n_channels = 2;

        let mut state = 0xabcdef01u64;
        let mut next = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((state >> 33) as f64) / ((1u64 << 31) as f64)
        };
        let mut x = Mat::zeros(n, 2);
        for i in 0..n {
            x[(i, 0)] = 1.0;
            x[(i, 1)] = next() * 2.0 - 1.0;
        }
        let mut y = Mat::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = 0.5 * x[(i, 1)] + (next() - 0.5) * 0.3;
        }
        let null = model::fit_glm(&y, &x);
        let pheno_mask = vec![true; n];
        let analysis = AnalysisVectors::from_null_model(&null, &pheno_mask).unwrap();

        let mut carriers: Vec<CarrierList> =
            (0..m).map(|_| CarrierList { entries: Vec::new() }).collect();
        for j in 0..m {
            for i in 0..n {
                if next() < 0.2 {
                    carriers[j].entries.push(CarrierEntry { sample_idx: i as u32, dosage: 1 });
                }
            }
            if carriers[j].entries.len() < 2 {
                carriers[j].entries.push(CarrierEntry { sample_idx: (j * 3) as u32, dosage: 1 });
                carriers[j].entries.push(CarrierEntry { sample_idx: (j * 3 + 1) as u32, dosage: 1 });
            }
        }

        let group: Vec<usize> = (0..n).map(|i| if i < n / 2 { 0 } else { 1 }).collect();
        let ancestry = AncestryInfo::new(group, n_pops, 3, 7590);
        let ann: Vec<Vec<f64>> = (0..n_channels)
            .map(|c| (0..m).map(|j| 0.3 + 0.1 * (c + j) as f64).collect())
            .collect();

        let from_carriers = run_ai_staar_gene(&carriers, &analysis, &ancestry, &ann, false);

        let g = sparse_score::carriers_to_dense(&carriers, &analysis);
        let pop_mafs = pop_mafs_from_carriers(&carriers, &analysis, &ancestry);
        let null2 = sparse_score::null_model_from_analysis(&analysis);
        let from_dense = run_ai_staar(&g, &ann, &pop_mafs, &ancestry, &null2, false);

        assert!((from_carriers.staar_o - from_dense.staar_o).abs() < 1e-12);
        assert!((from_carriers.acat_o - from_dense.acat_o).abs() < 1e-12);
        assert!((from_carriers.burden_1_25 - from_dense.burden_1_25).abs() < 1e-12);
        assert!((from_carriers.skat_1_25 - from_dense.skat_1_25).abs() < 1e-12);
    }

    #[test]
    fn flatten_unflatten_round_trip() {
        let n_channels = 2;
        let na = 1 + n_channels;
        let mut flat: Vec<f64> = (0..6 * na).map(|i| 0.01 * (i as f64 + 1.0)).collect();
        for v in flat.iter_mut() {
            *v = v.min(0.99);
        }
        let result = unflatten_staar(&flat, n_channels);
        let again = flatten_staar(&result, n_channels);
        assert_eq!(flat.len(), again.len());
        for i in 0..flat.len() {
            assert!((flat[i] - again[i]).abs() < 1e-12, "mismatch at {i}");
        }
    }
}
