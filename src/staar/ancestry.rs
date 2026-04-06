use faer::Mat;

use super::model::NullModel;
use super::score::{self, StaarResult};
use super::stats;

/// Ancestry-informed STAAR: uses population-specific allele frequencies
/// as weights to leverage allelic heterogeneity across ancestries.
///
/// For each population group, variants rare in that population get higher
/// weight, preventing rare-in-one/common-in-another dilution that reduces
/// power in multi-ancestry cohorts.
///
/// Reference: Li et al. (2024), AI-STAAR in xihaoli/STAAR
#[allow(dead_code)] // fields read by run_ai_staar when --ancestry-col is used
/// Per-population allele frequencies and group assignments.
pub struct AncestryInfo {
    /// Population group label for each sample (0-indexed).
    pub group: Vec<usize>,
    /// Number of distinct populations.
    pub n_pops: usize,
}

/// Run AI-STAAR for one gene.
///
/// Computes two scenarios:
///   Scenario 1: ancestry-ensemble weights (population-specific beta weights)
///   Scenario 2: annotation × ancestry weights
/// Final p-value: Cauchy combination of all tests across both scenarios.
///
/// `pop_mafs[variant][population]` = MAF in that population.
#[allow(dead_code)] // wired when --ai-staar lands (v0.2.0)
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

    if m == 0 || n_pops == 0 {
        let zero_mafs = vec![0.01; m];
        return score::run_staar(g, annotation_matrix, &zero_mafs, null, use_spa);
    }

    // Population sample counts for weighted pooling
    let mut pop_n = vec![0usize; n_pops];
    for &p in &ancestry.group { pop_n[p] += 1; }

    // Scenario 1: per-population genotype matrices with population-specific MAFs.
    // Zero out samples not in population p; run_staar applies beta(maf_p, 1, 25)
    // weighting internally via the population-specific MAFs.
    let mut scenario1_results: Vec<StaarResult> = Vec::with_capacity(n_pops);
    let mut g_pop = Mat::zeros(n, m);

    for p in 0..n_pops {
        for j in 0..m {
            for i in 0..n {
                g_pop[(i, j)] = if ancestry.group[i] == p { g[(i, j)] } else { 0.0 };
            }
        }
        let pop_mafs_p: Vec<f64> = pop_mafs.iter()
            .map(|v| v.get(p).copied().unwrap_or(0.01).clamp(1e-10, 0.499))
            .collect();
        scenario1_results.push(score::run_staar(
            &g_pop, annotation_matrix, &pop_mafs_p, null, use_spa,
        ));
    }

    // Scenario 2: standard STAAR with sample-size-weighted pooled MAF
    let total_n: f64 = pop_n.iter().sum::<usize>() as f64;
    let pooled_mafs: Vec<f64> = pop_mafs
        .iter()
        .map(|maf_vec| {
            let weighted: f64 = maf_vec.iter().enumerate()
                .map(|(p, &maf)| maf * pop_n[p] as f64)
                .sum();
            (weighted / total_n).clamp(1e-10, 0.499)
        })
        .collect();
    let scenario2 = score::run_staar(g, annotation_matrix, &pooled_mafs, null, use_spa);

    // Combine: Cauchy across all scenarios' STAAR-O values
    let mut all_staar_o: Vec<f64> = scenario1_results.iter().map(|r| r.staar_o).collect();
    all_staar_o.push(scenario2.staar_o);

    let mut combined = scenario2;
    combined.staar_o = stats::cauchy_combine(&all_staar_o);
    combined.acat_o = stats::cauchy_combine(
        &scenario1_results
            .iter()
            .map(|r| r.acat_o)
            .chain(std::iter::once(combined.acat_o))
            .collect::<Vec<_>>(),
    );

    combined
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::model;

    #[test]
    fn ai_staar_two_populations() {
        let n = 40;
        let m = 4;

        let mut g = Mat::zeros(n, m);
        for j in 0..m {
            for i in 0..2 {
                g[(i + j * 2, j)] = 1.0;
            }
        }

        let ann = vec![vec![0.5; m]; 3];

        // Two populations: first 20 samples = pop 0, last 20 = pop 1
        let ancestry = AncestryInfo {
            group: (0..n).map(|i| if i < 20 { 0 } else { 1 }).collect(),
            n_pops: 2,
        };

        // Population-specific MAFs: variant 0 rare in pop0 but common in pop1
        let pop_mafs = vec![
            vec![0.001, 0.05],
            vec![0.005, 0.005],
            vec![0.003, 0.003],
            vec![0.002, 0.002],
        ];

        let mut y = Mat::zeros(n, 1);
        for i in 0..n {
            y[(i, 0)] = if i % 5 == 0 { 1.0 } else { 0.0 };
        }
        let mut x = Mat::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }

        let null = model::fit_glm(&y, &x);
        let result = run_ai_staar(&g, &ann, &pop_mafs, &ancestry, &null, false);

        assert!(result.staar_o >= 0.0 && result.staar_o <= 1.0);
        assert!(result.burden_1_25 >= 0.0 && result.burden_1_25 <= 1.0);
    }
}
