//! Phenotype loading, null model fitting (GLM/logistic), and covariate management.

use std::collections::HashMap;
use std::path::Path;

use arrow::array::{Array, Float64Array};
use faer::prelude::*;
use faer::Mat;

use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::output::Output;
use crate::staar;
use crate::staar::genotype::GenotypeResult;
use crate::staar::kinship::KinshipState;

pub struct PhenotypeData {
    pub y: Mat<f64>,
    pub x: Mat<f64>,
    pub trait_type: staar::TraitType,
    pub n: usize,
    /// True for VCF samples that have phenotype data. Length = n_total VCF samples.
    /// Used to expand null model residuals back to VCF-indexed arrays.
    pub pheno_mask: Vec<bool>,
    /// Per-sample population labels and ensemble weights, when `--ancestry-col`
    /// is set. Same row order as `y` / `x`.
    pub ancestry: Option<staar::ancestry::AncestryInfo>,
}

/// Common aliases for the sample ID column.
const ID_ALIASES: &[&str] = &[
    "IID", "iid", "FID", "fid", "sample_id", "SampleID", "SAMPLE_ID",
    "sample", "Sample", "SAMPLE", "ID", "id", "SubjectID", "subject_id",
];

pub(crate) fn resolve_column(
    requested: &str,
    actual_cols: &[String],
    column_map: &HashMap<String, String>,
    role: &str,
) -> Result<String, CohortError> {
    if let Some(mapped) = column_map.get(requested) {
        if actual_cols.contains(mapped) {
            return Ok(mapped.clone());
        }
        return Err(CohortError::Input(format!(
            "Column map '{requested}={mapped}' but '{mapped}' not in phenotype. Available: {}",
            actual_cols.join(", ")
        )));
    }

    if actual_cols.contains(&requested.to_string()) {
        return Ok(requested.to_string());
    }

    let lower = requested.to_lowercase();
    let matches: Vec<&String> = actual_cols.iter()
        .filter(|c| c.to_lowercase() == lower)
        .collect();
    if matches.len() == 1 {
        return Ok(matches[0].clone());
    }

    Err(CohortError::Input(format!(
        "{role} '{requested}' not in phenotype. Available: {}",
        actual_cols.join(", ")
    )))
}

/// Intern an iterator of categorical labels into `(order, assignments)`
/// where `order` is the de-duplicated label list in first-seen order and
/// `assignments[i]` is the index of the i-th label in `order`. Shared
/// between `load_phenotype`'s ancestry path and
/// `kinship::load::load_groups`' group partition.
pub(crate) fn intern_labels<'a, I: IntoIterator<Item = &'a str>>(
    labels: I,
) -> (Vec<String>, Vec<usize>) {
    let mut order: Vec<String> = Vec::new();
    let mut assignments: Vec<usize> = Vec::new();
    for label in labels {
        let idx = order
            .iter()
            .position(|l| l == label)
            .unwrap_or_else(|| {
                order.push(label.to_string());
                order.len() - 1
            });
        assignments.push(idx);
    }
    (order, assignments)
}

pub(crate) fn resolve_id_column(
    actual_cols: &[String],
    column_map: &HashMap<String, String>,
) -> String {
    if let Some(mapped) = column_map.get("id") {
        if actual_cols.contains(mapped) {
            return mapped.clone();
        }
    }

    let first = &actual_cols[0];
    let first_lower = first.to_lowercase();
    if ID_ALIASES.iter().any(|a| a.to_lowercase() == first_lower) {
        return first.clone();
    }

    for alias in ID_ALIASES {
        if let Some(col) = actual_cols.iter().find(|c| c.as_str() == *alias) {
            return col.clone();
        }
    }

    first.clone()
}

/// Extract a string from any Arrow column (handles Utf8, Utf8View, etc.).
pub(crate) fn arrow_str(col: &dyn Array, row: usize) -> String {
    arrow::util::display::array_value_to_string(col, row).unwrap_or_default()
}

/// Extract a f64 from any Arrow numeric column (handles Float64, Int64, etc.).
fn arrow_f64(col: &dyn Array, row: usize) -> f64 {
    if col.is_null(row) { return f64::NAN; }
    if let Some(a) = col.as_any().downcast_ref::<Float64Array>() {
        return a.value(row);
    }
    arrow::util::display::array_value_to_string(col, row)
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(f64::NAN)
}

#[allow(clippy::too_many_arguments)]
pub fn load_phenotype(
    engine: &DfEngine,
    phenotype: &Path,
    covariates: &[String],
    geno: &GenotypeResult,
    trait_name: &str,
    ancestry_col: Option<&str>,
    ai_base_tests: usize,
    ai_seed: u64,
    column_map: &HashMap<String, String>,
    out: &dyn Output,
) -> Result<PhenotypeData, CohortError> {
    out.status("Loading phenotype...");

    engine.register_csv("_pheno", phenotype, b'\t')?;
    let pheno_cols = engine.table_columns("_pheno")?;

    let trait_col = resolve_column(trait_name, &pheno_cols, column_map, "Trait")?;
    if trait_col != trait_name {
        out.status(&format!("  Resolved trait '{trait_name}' -> column '{trait_col}'"));
    }

    let mut resolved_covs: Vec<String> = Vec::with_capacity(covariates.len());
    for cov in covariates {
        let resolved = resolve_column(cov, &pheno_cols, column_map, "Covariate")?;
        if resolved != *cov {
            out.status(&format!("  Resolved covariate '{cov}' -> column '{resolved}'"));
        }
        resolved_covs.push(resolved);
    }

    let distinct_count = engine.query_scalar(&format!(
        "SELECT COUNT(DISTINCT \"{trait_col}\") FROM _pheno"
    ))?;
    let trait_type = if distinct_count <= 2 {
        staar::TraitType::Binary
    } else {
        staar::TraitType::Continuous
    };
    out.status(&format!("  Trait '{trait_col}' -> {:?}", trait_type));

    let id_col = resolve_id_column(&pheno_cols, column_map);
    let n_cov = resolved_covs.len();
    let cov_select = if n_cov == 0 { String::new() }
        else { format!(", {}", resolved_covs.iter()
            .map(|c| format!("CAST(p.\"{c}\" AS DOUBLE)"))
            .collect::<Vec<_>>().join(", ")) };

    let resolved_ancestry: Option<String> = match ancestry_col {
        Some(name) => Some(resolve_column(name, &pheno_cols, column_map, "Ancestry")?),
        None => None,
    };
    let ancestry_select = match resolved_ancestry.as_deref() {
        Some(c) => format!(", CAST(p.\"{c}\" AS VARCHAR)"),
        None => String::new(),
    };

    let sample_list = geno.sample_names.iter()
        .map(|s| format!("'{s}'")).collect::<Vec<_>>().join(",");

    // Include ID column — we reorder rows to match genotype sample order below.
    let pheno_sql = format!(
        "SELECT p.\"{id_col}\", CAST(p.\"{trait_col}\" AS DOUBLE) {cov_select} {ancestry_select} \
         FROM _pheno p \
         WHERE p.\"{id_col}\" IN ({sample_list}) AND p.\"{trait_col}\" IS NOT NULL"
    );

    let batches = engine.collect(&pheno_sql)?;

    // Collect phenotype keyed by sample ID — SQL row order is irrelevant.
    let ancestry_col_idx = 2 + n_cov;
    let mut pheno_map: HashMap<String, (f64, Vec<f64>, Option<String>)> = HashMap::new();
    for batch in &batches {
        for row in 0..batch.num_rows() {
            let id = arrow_str(batch.column(0).as_ref(), row);
            let y_val = arrow_f64(batch.column(1).as_ref(), row);
            let covs: Vec<f64> = (0..n_cov)
                .map(|j| arrow_f64(batch.column(2 + j).as_ref(), row))
                .collect();
            let pop_label = if resolved_ancestry.is_some() {
                let raw = arrow_str(batch.column(ancestry_col_idx).as_ref(), row);
                if raw.is_empty() { None } else { Some(raw) }
            } else {
                None
            };
            pheno_map.insert(id, (y_val, covs, pop_label));
        }
    }

    // Build y, x in genotype sample order — CRITICAL for correct p-values.
    // Dosages are indexed by VCF header position; phenotype must match exactly.
    // Samples without phenotype get NaN y and zero covariates — they are excluded
    // from the null model but keep their VCF index position so carrier indices
    // stay valid.
    let n_total = geno.sample_names.len();
    let mut y_vec = Vec::with_capacity(n_total);
    let mut x_vecs: Vec<Vec<f64>> = vec![Vec::with_capacity(n_total); n_cov];
    let mut pheno_mask = Vec::with_capacity(n_total);
    let mut pop_labels: Vec<Option<String>> = Vec::with_capacity(n_total);
    for name in &geno.sample_names {
        match pheno_map.get(name.as_str()) {
            Some((y_val, covs, pop_label)) => {
                y_vec.push(*y_val);
                for (j, c) in covs.iter().enumerate() {
                    x_vecs[j].push(*c);
                }
                pheno_mask.push(true);
                if resolved_ancestry.is_some() {
                    pop_labels.push(pop_label.clone());
                }
            }
            None => {
                pheno_mask.push(false);
            }
        }
    }

    let n_missing = pheno_mask.iter().filter(|v| !*v).count();
    if n_missing > 0 {
        out.warn(&format!(
            "{n_missing} of {n_total} genotype samples not in phenotype — excluded from analysis"
        ));
    }

    let n = n_total - n_missing;
    out.status(&format!("  {} samples with phenotype + genotype", n));
    if n < 10 {
        return Err(CohortError::Analysis(format!(
            "Only {n} samples with both phenotype and genotype data (need >= 10). \
             Check that sample IDs in '{}' match the VCF header.",
            phenotype.display()
        )));
    }

    // NaN guards — NaN trait values silently corrupt the null model.
    let nan_y = y_vec.iter().filter(|v| v.is_nan()).count();
    if nan_y > 0 {
        return Err(CohortError::Input(format!(
            "{nan_y} samples have NaN trait values in column '{trait_col}'. \
             Remove or impute these before running STAAR."
        )));
    }
    for (j, cov_vals) in x_vecs.iter().enumerate() {
        let nan_count = cov_vals.iter().filter(|v| v.is_nan()).count();
        if nan_count > 0 {
            return Err(CohortError::Input(format!(
                "{nan_count} samples have missing values for covariate '{}'. \
                 Remove samples with missing covariates or impute before running STAAR.",
                resolved_covs[j]
            )));
        }
    }

    let mut y = Mat::zeros(n, 1);
    let mut x = Mat::zeros(n, 1 + n_cov);
    for i in 0..n {
        y[(i, 0)] = y_vec[i];
        x[(i, 0)] = 1.0; // intercept
        for j in 0..n_cov {
            x[(i, j + 1)] = x_vecs[j][i];
        }
    }

    let ancestry = match resolved_ancestry.as_deref() {
        None => None,
        Some(col) => {
            let missing = pop_labels.iter().filter(|p| p.is_none()).count();
            if missing > 0 {
                return Err(CohortError::Input(format!(
                    "{missing} samples have missing values for ancestry column '{col}'. \
                     Remove or impute before running --ancestry-col."
                )));
            }
            let (order, group) =
                intern_labels(pop_labels.iter().flatten().map(|s| s.as_str()));
            if order.len() < 2 {
                return Err(CohortError::Input(format!(
                    "Ancestry column '{col}' has only {} distinct value(s); \
                     --ancestry-col needs at least 2 populations.",
                    order.len()
                )));
            }
            out.status(&format!(
                "  Ancestry: {} populations ({}), AI-STAAR B={ai_base_tests}",
                order.len(),
                order.join(", "),
            ));
            Some(staar::ancestry::AncestryInfo::new(
                group,
                order.len(),
                ai_base_tests,
                ai_seed,
            ))
        }
    };

    Ok(PhenotypeData { y, x, trait_type, n, pheno_mask, ancestry })
}

pub fn load_known_loci(
    engine: &DfEngine,
    geno: &GenotypeResult,
    loci_path: &Path,
    n_samples: usize,
    out: &dyn Output,
) -> Result<Mat<f64>, CohortError> {
    use crate::types::Chromosome;

    let content = std::fs::read_to_string(loci_path)
        .map_err(|e| CohortError::Resource(format!(
            "Cannot read known loci file '{}': {e}", loci_path.display()
        )))?;
    let loci: Vec<(String, i32)> = content.lines()
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .filter_map(|l| {
            let parts: Vec<&str> = l.split(':').collect();
            if parts.len() >= 2 {
                let chrom: Chromosome = parts[0].parse().ok()?;
                parts[1].parse::<i32>().ok().map(|pos| (chrom.label(), pos))
            } else { None }
        })
        .collect();

    if loci.is_empty() {
        return Err(CohortError::Input(format!(
            "Known loci file '{}' is empty or unparseable.", loci_path.display()
        )));
    }

    out.status(&format!("  Loading {} known loci for conditional analysis...", loci.len()));

    let mut by_chrom: std::collections::HashMap<&str, Vec<(usize, i32)>> = std::collections::HashMap::new();
    for (col, (chrom, pos)) in loci.iter().enumerate() {
        by_chrom.entry(chrom.as_str()).or_default().push((col, *pos));
    }

    let mut x_cond = Mat::zeros(n_samples, loci.len());

    for (chrom, chrom_loci) in &by_chrom {
        let geno_path = format!("{}/chromosome={chrom}/data.parquet", geno.output_dir.display());
        let positions: Vec<u32> = chrom_loci.iter().map(|(_, p)| *p as u32).collect();
        let flat = crate::staar::genotype::load(engine, &geno_path, &positions, n_samples)?;

        let pos_to_cols: std::collections::HashMap<i32, Vec<usize>> = {
            let mut m: std::collections::HashMap<i32, Vec<usize>> = std::collections::HashMap::new();
            for &(col, pos) in chrom_loci {
                m.entry(pos).or_default().push(col);
            }
            m
        };

        for (flat_idx, &pos) in positions.iter().enumerate() {
            if let Some(cols) = pos_to_cols.get(&(pos as i32)) {
                for &col in cols {
                    let base = flat_idx * n_samples;
                    for i in 0..n_samples {
                        x_cond[(i, col)] = flat[base + i];
                    }
                }
            }
        }
    }
    Ok(x_cond)
}

pub fn augment_covariates(x: &Mat<f64>, x_cond: &Mat<f64>) -> Mat<f64> {
    let n = x.nrows();
    let k_orig = x.ncols();
    let k_cond = x_cond.ncols();
    let mut x_new = Mat::zeros(n, k_orig + k_cond);
    for i in 0..n {
        for j in 0..k_orig { x_new[(i, j)] = x[(i, j)]; }
        for j in 0..k_cond { x_new[(i, k_orig + j)] = x_cond[(i, j)]; }
    }
    x_new
}

/// Fitted null model — immutable after construction.
///
/// Stores components of the projection matrix, NOT the full n×n matrix.
/// For continuous traits: P = I - X(X'X)^{-1}X'
/// For binary traits:     P = W - WX(X'WX)^{-1}X'W  where W = diag(μ(1-μ))
///
/// The kernel K = G'PG is computed on-the-fly which is O(n*k*m)
/// instead of O(n^2*m) and avoids materializing a 50K×50K matrix.
pub struct NullModel {
    pub residuals: Mat<f64>,
    pub x_matrix: Mat<f64>,
    /// (X'X)^{-1} for continuous, (X'WX)^{-1} for binary, (X'Σ⁻¹X)⁻¹ when kinship is set.
    pub xtx_inv: Mat<f64>,
    pub sigma2: f64,
    pub n_samples: usize,
    /// Fitted probabilities μ_i = P(Y_i = 1) under the null. Binary traits only.
    pub fitted_values: Option<Vec<f64>>,
    /// Working weights W_i = μ_i(1 - μ_i). Binary traits only.
    pub working_weights: Option<Vec<f64>>,
    /// Kinship-aware variance components and projection. Set when `--kinship`
    /// is in use; the score path then dispatches to the kinship-aware kernel.
    pub kinship: Option<KinshipState>,
}

impl NullModel {
    /// Compute (I - H) * G for continuous traits: G - X(X'X)^{-1}X'G
    fn project(&self, g: &Mat<f64>) -> Mat<f64> {
        let n = g.nrows();
        let m = g.ncols();
        assert_eq!(n, self.n_samples, "G rows must match n_samples");

        let xt_g = self.x_matrix.transpose() * g;
        let xtx_inv_xt_g = &self.xtx_inv * &xt_g;
        let mut h_g = &self.x_matrix * &xtx_inv_xt_g;

        for j in 0..m {
            for i in 0..n {
                h_g[(i, j)] = g[(i, j)] - h_g[(i, j)];
            }
        }
        h_g
    }

    /// Compute the score test covariance kernel K = G'PG.
    ///
    /// For continuous: K = G'(I-H)G where H = X(X'X)^{-1}X'.
    /// For binary:     K = G'WG - G'WX(X'WX)^{-1}X'WG where W = diag(μ(1-μ)).
    pub fn compute_kernel(&self, g: &Mat<f64>) -> Mat<f64> {
        let n = g.nrows();
        let m = g.ncols();
        assert_eq!(n, self.n_samples, "G rows must match n_samples");

        if let Some(ref w) = self.working_weights {
            // Binary: K = G'WG - G'WX(X'WX)^{-1}X'WG
            let mut wg = Mat::zeros(n, m);
            for j in 0..m {
                for i in 0..n {
                    wg[(i, j)] = w[i] * g[(i, j)];
                }
            }
            let gtwg = g.transpose() * &wg;
            let xtwg = self.x_matrix.transpose() * &wg;
            let inv_xtwg = &self.xtx_inv * &xtwg;
            let x_inv_xtwg = &self.x_matrix * &inv_xtwg;
            let correction = wg.transpose() * &x_inv_xtwg;

            let mut k = gtwg;
            for i in 0..m {
                for j in 0..m {
                    k[(i, j)] -= correction[(i, j)];
                }
            }
            k
        } else {
            // Continuous: K = G'(I-H)G
            let pg = self.project(g);
            g.transpose() * &pg
        }
    }
}

/// Fit GLM null model for continuous trait via QR decomposition.
///
/// y = X * beta + epsilon
///
/// Caller must prepend an intercept column to X if desired.
pub fn fit_glm(y: &Mat<f64>, x: &Mat<f64>) -> NullModel {
    let n = y.nrows();
    let k = x.ncols();
    assert_eq!(x.nrows(), n, "X rows must match y rows");
    assert!(n > k, "Need more samples than covariates");

    // Solve via normal equations: beta = (X'X)^-1 X'y
    // For n >> k (typical: 50K samples, 10 covariates), this is efficient and stable.
    let xtx = x.transpose() * x;           // k × k
    let xty = x.transpose() * y;           // k × 1
    let eye_k = Mat::<f64>::identity(k, k);
    let xtx_inv = xtx.col_piv_qr().solve(&eye_k); // (X'X)^-1
    let beta = &xtx_inv * &xty;            // beta = (X'X)^-1 X'y

    // Residuals: r = y - X * beta
    let y_hat = x * &beta;
    let mut residuals = Mat::zeros(n, 1);
    for i in 0..n {
        residuals[(i, 0)] = y[(i, 0)] - y_hat[(i, 0)];
    }

    // Residual variance: sigma2 = ||r||^2 / (n - k)
    let rss: f64 = (0..n).map(|i| residuals[(i, 0)].powi(2)).sum();
    let sigma2 = rss / (n - k) as f64;

    NullModel {
        residuals,
        x_matrix: x.to_owned(),
        xtx_inv,
        sigma2,
        n_samples: n,
        fitted_values: None,
        working_weights: None,
        kinship: None,
    }
}

/// Fit logistic null model for binary trait via IRLS.
///
/// y = 0/1, logit(P(y=1)) = X * beta
pub fn fit_logistic(y: &Mat<f64>, x: &Mat<f64>, max_iter: usize) -> NullModel {
    let n = y.nrows();
    let k = x.ncols();
    let max_iter = if max_iter == 0 { 25 } else { max_iter };

    let mut beta = Mat::zeros(k, 1);
    // Reusable IRLS scratch — allocated once before the loop, zeroed each
    // iteration. Allocating fresh per iteration was costing ~2,700 allocs
    // per fit on a typical k=10 / max_iter=25 run.
    let mut mu = Mat::<f64>::zeros(n, 1);
    let mut w_diag = vec![0.0_f64; n];
    let mut xtwx: Mat<f64> = Mat::zeros(k, k);
    let mut xtwy: Mat<f64> = Mat::zeros(k, 1);

    for _ in 0..max_iter {
        // mu = sigmoid(X * beta)
        let eta = x * &beta;
        for i in 0..n {
            let eta_i = eta[(i, 0)].clamp(-500.0, 500.0); // prevent exp overflow
            let p = 1.0 / (1.0 + (-eta_i).exp());
            mu[(i, 0)] = p;
            w_diag[i] = p * (1.0 - p);
        }

        // Working response: z = eta + (y - mu) / w
        // Weighted least squares: X' * W * X * delta = X' * W * z
        for j in 0..k {
            xtwy[(j, 0)] = 0.0;
            for l in 0..k {
                xtwx[(j, l)] = 0.0;
            }
        }
        for i in 0..n {
            let wi = w_diag[i].max(1e-10);
            let zi = eta[(i, 0)] + (y[(i, 0)] - mu[(i, 0)]) / wi;
            for j in 0..k {
                xtwy[(j, 0)] += x[(i, j)] * wi * zi;
                for l in 0..k {
                    xtwx[(j, l)] += x[(i, j)] * wi * x[(i, l)];
                }
            }
        }

        let new_beta: Mat<f64> = xtwx.col_piv_qr().solve(&xtwy);

        let mut max_diff = 0.0f64;
        for j in 0..k {
            let diff: f64 = new_beta[(j, 0)] - beta[(j, 0)];
            max_diff = max_diff.max(diff.abs());
        }
        beta = new_beta;
        if max_diff < 1e-8 {
            break;
        }
    }

    // Final residuals, fitted probabilities, and working weights.
    // Score test uses U = G'r where r = Y - μ (raw score residuals, NOT working).
    // SPA needs the fitted μ_i to compute the exact CGF of the score distribution.
    // The W-weighted projection (X'WX)^{-1} is needed for correct SKAT covariance.
    let eta = x * &beta;
    let mut residuals = Mat::zeros(n, 1);
    let mut fitted = Vec::with_capacity(n);
    let mut w_final = Vec::with_capacity(n);
    for j in 0..k {
        for l in 0..k {
            xtwx[(j, l)] = 0.0;
        }
    }
    for i in 0..n {
        let eta_i = eta[(i, 0)].clamp(-500.0, 500.0);
        let mu_i = 1.0 / (1.0 + (-eta_i).exp());
        let wi = mu_i * (1.0 - mu_i);
        fitted.push(mu_i);
        w_final.push(wi);
        residuals[(i, 0)] = y[(i, 0)] - mu_i;
        for j in 0..k {
            for l in 0..k {
                xtwx[(j, l)] += x[(i, j)] * wi * x[(i, l)];
            }
        }
    }

    let eye_k = Mat::<f64>::identity(k, k);
    let xtwx_inv = xtwx.col_piv_qr().solve(&eye_k);

    NullModel {
        residuals,
        x_matrix: x.to_owned(),
        xtx_inv: xtwx_inv,
        sigma2: 1.0,
        n_samples: n,
        fitted_values: Some(fitted),
        working_weights: Some(w_final),
        kinship: None,
    }
}
