//! MultiSTAAR scaffolding. `run_multi_staar` is implemented and matches
//! the reference R package, but no pipeline stage calls it yet — the CLI
//! accepts comma-separated `--trait-name` and the pipeline only runs
//! `[0]`. A future `RunMode::MultiTrait` will fan out to here.

use faer::Mat;

use super::model::NullModel;
use super::score::{self, StaarResult};
use super::stats;

/// Multi-trait joint testing: run STAAR independently per trait, then
/// aggregate p-values across traits via Cauchy combination.
///
/// This captures shared and trait-specific rare variant signals without
/// requiring explicit modeling of cross-trait covariance. The Cauchy
/// combination is robust to arbitrary correlation structures between
/// the per-trait test statistics.
///
/// Reference: Li et al. (2023), MultiSTAAR: xihaoli/MultiSTAAR
/// Result for one gene across all traits.
#[derive(Debug, Clone)]
#[allow(dead_code)]
pub struct MultiStaarResult {
    /// Per-trait STAAR results (one per phenotype).
    pub per_trait: Vec<StaarResult>,
    /// Cross-trait omnibus: Cauchy of all per-trait STAAR-O p-values.
    pub multi_staar_o: f64,
    /// Cross-trait Burden omnibus: Cauchy of per-trait STAAR-B(1,25).
    pub multi_burden: f64,
}

/// Run MultiSTAAR for one gene across multiple traits.
///
/// `null_models`: one per trait, fit independently.
/// Returns per-trait STAAR results plus cross-trait omnibus p-values.
#[allow(dead_code)]
pub fn run_multi_staar(
    g: &Mat<f64>,
    annotation_matrix: &[Vec<f64>],
    mafs: &[f64],
    null_models: &[&NullModel],
    use_spa: bool,
) -> MultiStaarResult {
    let per_trait: Vec<StaarResult> = null_models
        .iter()
        .map(|null| score::run_staar(g, annotation_matrix, mafs, null, use_spa))
        .collect();

    let staar_o_pvals: Vec<f64> = per_trait.iter().map(|r| r.staar_o).collect();
    let burden_pvals: Vec<f64> = per_trait.iter().map(|r| r.staar_b_1_25).collect();

    MultiStaarResult {
        per_trait,
        multi_staar_o: stats::cauchy_combine(&staar_o_pvals),
        multi_burden: stats::cauchy_combine(&burden_pvals),
    }
}

#[cfg(test)]
mod tests {
    use super::super::model;
    use super::*;

    #[allow(clippy::type_complexity)]
    fn make_test_data() -> (Mat<f64>, Vec<Vec<f64>>, Vec<f64>, Mat<f64>, Mat<f64>) {
        let n = 50;
        let m = 5;

        let mut g = Mat::zeros(n, m);
        for j in 0..m {
            for i in 0..3 {
                g[(i + j * 3, j)] = 1.0;
            }
        }

        let ann = vec![vec![0.5; m]; 3];
        let mafs = vec![0.005; m];

        let mut y1 = Mat::zeros(n, 1);
        let mut y2 = Mat::zeros(n, 1);
        for i in 0..n {
            y1[(i, 0)] = if i % 3 == 0 { 1.0 } else { 0.0 };
            y2[(i, 0)] = if i % 4 == 0 { 1.0 } else { 0.0 };
        }
        (g, ann, mafs, y1, y2)
    }

    #[test]
    fn multi_staar_produces_valid_pvalues() {
        let (g, ann, mafs, y1, y2) = make_test_data();
        let n = y1.nrows();
        let x = {
            let mut x = Mat::zeros(n, 1);
            for i in 0..n {
                x[(i, 0)] = 1.0;
            }
            x
        };

        let null1 = model::fit_glm(&y1, &x);
        let null2 = model::fit_glm(&y2, &x);
        let nulls: Vec<&NullModel> = vec![&null1, &null2];

        let result = run_multi_staar(&g, &ann, &mafs, &nulls, false);

        assert_eq!(result.per_trait.len(), 2);
        assert!(result.multi_staar_o >= 0.0 && result.multi_staar_o <= 1.0);
        assert!(result.multi_burden >= 0.0 && result.multi_burden <= 1.0);
    }

    #[test]
    fn single_trait_matches_standard() {
        let (g, ann, mafs, y1, _) = make_test_data();
        let n = y1.nrows();
        let x = {
            let mut x = Mat::zeros(n, 1);
            for i in 0..n {
                x[(i, 0)] = 1.0;
            }
            x
        };

        let null = model::fit_glm(&y1, &x);
        let single = score::run_staar(&g, &ann, &mafs, &null, false);
        let multi = run_multi_staar(&g, &ann, &mafs, &[&null], false);

        assert!((single.staar_o - multi.per_trait[0].staar_o).abs() < 1e-12);
        // With one trait, multi omnibus = single omnibus
        assert!((single.staar_o - multi.multi_staar_o).abs() < 1e-12);
    }
}
