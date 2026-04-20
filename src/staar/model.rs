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
pub(crate) fn arrow_f64(col: &dyn Array, row: usize) -> f64 {
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
    trait_names: &[String],
    ancestry_col: Option<&str>,
    ai_base_tests: usize,
    ai_seed: u64,
    column_map: &HashMap<String, String>,
    out: &dyn Output,
) -> Result<PhenotypeData, CohortError> {
    out.status("Loading phenotype...");
    assert!(
        !trait_names.is_empty(),
        "load_phenotype requires at least one trait column; \
         build_config should have rejected an empty --trait-name list"
    );

    engine.register_csv("_pheno", phenotype, b'\t')?;
    let pheno_cols = engine.table_columns("_pheno")?;

    // Resolve every requested trait column up front so the SQL below can
    // reference them by position. `k=1` is the single-trait case expressed
    // as a length-1 vec; there is no separate loader.
    let mut trait_cols: Vec<String> = Vec::with_capacity(trait_names.len());
    for t in trait_names {
        let resolved = resolve_column(t, &pheno_cols, column_map, "Trait")?;
        if resolved != *t {
            out.status(&format!("  Resolved trait '{t}' -> column '{resolved}'"));
        }
        trait_cols.push(resolved);
    }
    let k_traits = trait_cols.len();

    let mut resolved_covs: Vec<String> = Vec::with_capacity(covariates.len());
    for cov in covariates {
        let resolved = resolve_column(cov, &pheno_cols, column_map, "Covariate")?;
        if resolved != *cov {
            out.status(&format!("  Resolved covariate '{cov}' -> column '{resolved}'"));
        }
        resolved_covs.push(resolved);
    }

    // Trait type detection runs against the first resolved column and
    // applies to the entire joint multi-trait run. Multi-trait joint STAAR
    // now supports all-binary traits via per-trait logistic IRLS +
    // cross-trait residual correlation (see `multi::MultiNull::fit_binary`);
    // we enforce that every trait in a multi-trait run shares the same
    // detected family, since the joint null cannot mix gaussian and binomial.
    let primary_trait = &trait_cols[0];
    let distinct_count = engine.query_scalar(&format!(
        "SELECT COUNT(DISTINCT \"{primary_trait}\") FROM _pheno"
    ))?;
    let trait_type = if distinct_count <= 2 {
        staar::TraitType::Binary
    } else {
        staar::TraitType::Continuous
    };
    out.status(&format!("  Trait '{primary_trait}' -> {:?}", trait_type));
    if k_traits > 1 {
        for extra in &trait_cols[1..] {
            let dc = engine.query_scalar(&format!(
                "SELECT COUNT(DISTINCT \"{extra}\") FROM _pheno"
            ))?;
            let tt = if dc <= 2 {
                staar::TraitType::Binary
            } else {
                staar::TraitType::Continuous
            };
            if tt != trait_type {
                return Err(CohortError::Input(format!(
                    "Multi-trait joint STAAR requires every --trait-name to \
                     share the same family, but '{primary_trait}' is {:?} \
                     and '{extra}' is {:?}. Run each family separately.",
                    trait_type, tt,
                )));
            }
        }
    }

    let id_col = resolve_id_column(&pheno_cols, column_map);
    let n_cov = resolved_covs.len();

    // Column layout in the SQL result, in order: id, k trait casts, n_cov
    // covariate casts, optional ancestry. The row collection loop below
    // reads columns by fixed positional offsets that depend on k_traits.
    let trait_select = trait_cols
        .iter()
        .map(|c| format!("CAST(p.\"{c}\" AS DOUBLE)"))
        .collect::<Vec<_>>()
        .join(", ");
    let cov_select = if n_cov == 0 {
        String::new()
    } else {
        format!(
            ", {}",
            resolved_covs
                .iter()
                .map(|c| format!("CAST(p.\"{c}\" AS DOUBLE)"))
                .collect::<Vec<_>>()
                .join(", ")
        )
    };

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

    // Complete-cases semantics: a sample is dropped if ANY requested trait
    // is NULL. For `k=1` this is the original single-trait filter; for
    // `k>1` it matches the R MultiSTAAR convention of fitting the joint
    // null model only on samples with observations for every trait.
    let trait_not_null = trait_cols
        .iter()
        .map(|c| format!("p.\"{c}\" IS NOT NULL"))
        .collect::<Vec<_>>()
        .join(" AND ");

    // Include ID column — we reorder rows to match genotype sample order below.
    let pheno_sql = format!(
        "SELECT p.\"{id_col}\", {trait_select} {cov_select} {ancestry_select} \
         FROM _pheno p \
         WHERE p.\"{id_col}\" IN ({sample_list}) AND {trait_not_null}"
    );

    let batches = engine.collect(&pheno_sql)?;

    // One row of the joined phenotype keyed by sample ID. Named fields
    // beat a three-tuple both for readability and because clippy flags
    // `HashMap<String, (Vec<f64>, Vec<f64>, Option<String>)>` as a type
    // complexity warning.
    struct PhenoRow {
        traits: Vec<f64>,
        covs: Vec<f64>,
        ancestry: Option<String>,
    }

    // Collect phenotype keyed by sample ID — SQL row order is irrelevant.
    // Column offsets: id=0, traits=1..=k, covariates=(1+k)..(1+k+n_cov),
    // ancestry=(1+k+n_cov).
    let trait_base = 1;
    let cov_base = trait_base + k_traits;
    let ancestry_col_idx = cov_base + n_cov;
    let mut pheno_map: HashMap<String, PhenoRow> = HashMap::new();
    for batch in &batches {
        for row in 0..batch.num_rows() {
            let id = arrow_str(batch.column(0).as_ref(), row);
            let traits: Vec<f64> = (0..k_traits)
                .map(|c| arrow_f64(batch.column(trait_base + c).as_ref(), row))
                .collect();
            let covs: Vec<f64> = (0..n_cov)
                .map(|j| arrow_f64(batch.column(cov_base + j).as_ref(), row))
                .collect();
            let ancestry = if resolved_ancestry.is_some() {
                let raw = arrow_str(batch.column(ancestry_col_idx).as_ref(), row);
                if raw.is_empty() { None } else { Some(raw) }
            } else {
                None
            };
            pheno_map.insert(id, PhenoRow { traits, covs, ancestry });
        }
    }

    // Build y, x in genotype sample order — CRITICAL for correct p-values.
    // Dosages are indexed by VCF header position; phenotype must match exactly.
    // Samples without phenotype are excluded from the null model but keep
    // their VCF index position so carrier indices stay valid.
    let n_total = geno.sample_names.len();
    let mut y_rows: Vec<Vec<f64>> = Vec::with_capacity(n_total);
    let mut x_vecs: Vec<Vec<f64>> = vec![Vec::with_capacity(n_total); n_cov];
    let mut pheno_mask = Vec::with_capacity(n_total);
    let mut pop_labels: Vec<Option<String>> = Vec::with_capacity(n_total);
    for name in &geno.sample_names {
        match pheno_map.get(name.as_str()) {
            Some(row) => {
                y_rows.push(row.traits.clone());
                for (j, c) in row.covs.iter().enumerate() {
                    x_vecs[j].push(*c);
                }
                pheno_mask.push(true);
                if resolved_ancestry.is_some() {
                    pop_labels.push(row.ancestry.clone());
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

    // NaN guards — NaN trait values silently corrupt the null model. Scan
    // each trait column and report the first one with NaN by name so the
    // operator knows which column to fix.
    for (c, col_name) in trait_cols.iter().enumerate() {
        let nan_count = y_rows.iter().filter(|row| row[c].is_nan()).count();
        if nan_count > 0 {
            return Err(CohortError::Input(format!(
                "{nan_count} samples have NaN trait values in column '{col_name}'. \
                 Remove or impute these before running STAAR."
            )));
        }
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

    let mut y = Mat::zeros(n, k_traits);
    let mut x = Mat::zeros(n, 1 + n_cov);
    for i in 0..n {
        for c in 0..k_traits {
            y[(i, c)] = y_rows[i][c];
        }
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
    /// SCANG-side state populated lazily once a run needs Monte Carlo
    /// thresholds. See `crate::staar::scang::ScangExt`.
    pub scang: Option<crate::staar::scang::ScangExt>,
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
        scang: None,
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
        scang: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::DfEngine;
    use crate::output::Progress;
    use crate::resource::Resources;
    use crate::staar::TraitType;
    use std::io::Write;
    use std::path::PathBuf;

    /// Test sink: Output trait impl that swallows every message. The
    /// loader calls `status`/`warn` for progress reporting; tests don't
    /// care about the text, only about the returned `PhenotypeData`.
    struct SilentOut;
    impl crate::output::Output for SilentOut {
        fn status(&self, _msg: &str) {}
        fn success(&self, _msg: &str) {}
        fn warn(&self, _msg: &str) {}
        fn error(&self, _err: &CohortError) {}
        fn result_json(&self, _data: &serde_json::Value) {}
        fn table(&self, _headers: &[&str], _rows: &[Vec<String>]) {}
        fn progress(&self, _total: u64, _label: &str) -> Progress {
            Progress::noop()
        }
    }

    fn fake_geno(sample_ids: &[&str]) -> GenotypeResult {
        GenotypeResult {
            sample_names: sample_ids.iter().map(|s| (*s).to_string()).collect(),
            output_dir: PathBuf::from("/unused-in-load-phenotype"),
        }
    }

    /// Write a minimal phenotype TSV with `IID, trait columns, age, sex`
    /// to a tempdir and return the path plus tempdir (caller keeps the
    /// dir alive for the duration of the test).
    fn write_pheno_tsv(
        dir: &tempfile::TempDir,
        header: &[&str],
        rows: &[Vec<&str>],
    ) -> PathBuf {
        let path = dir.path().join("pheno.tsv");
        let mut f = std::fs::File::create(&path).unwrap();
        writeln!(f, "{}", header.join("\t")).unwrap();
        for row in rows {
            writeln!(f, "{}", row.join("\t")).unwrap();
        }
        path
    }

    fn test_engine() -> DfEngine {
        DfEngine::new(&Resources::detect()).expect("DfEngine for loader tests")
    }

    /// 12 samples, two continuous traits, shared covariates. The k=2 load
    /// must produce a `(12, 2)` Y with the values in genotype sample order
    /// (which is the identity permutation here) and two rows of X.
    #[test]
    fn load_phenotype_two_traits_reads_both_columns() {
        let dir = tempfile::tempdir().unwrap();
        let sample_ids: Vec<String> = (1..=12).map(|i| format!("s{i}")).collect();
        let rows: Vec<Vec<String>> = sample_ids
            .iter()
            .enumerate()
            .map(|(i, id)| {
                let bmi = 20.0 + i as f64; // 20..31, 12 distinct → Continuous
                let height = 150.0 + 3.0 * i as f64;
                let age = 30.0 + i as f64;
                let sex = (i % 2) as f64;
                vec![
                    id.clone(),
                    format!("{bmi}"),
                    format!("{height}"),
                    format!("{age}"),
                    format!("{sex}"),
                ]
            })
            .collect();
        let row_refs: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(|s| s.as_str()).collect())
            .collect();
        let pheno_path = write_pheno_tsv(
            &dir,
            &["IID", "BMI", "HEIGHT", "age", "sex"],
            &row_refs,
        );
        let geno = fake_geno(
            &sample_ids
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
        );

        let pheno = load_phenotype(
            &test_engine(),
            &pheno_path,
            &["age".into(), "sex".into()],
            &geno,
            &["BMI".into(), "HEIGHT".into()],
            None,
            5,
            42,
            &HashMap::new(),
            &SilentOut,
        )
        .expect("multi-trait load should succeed");

        assert_eq!(pheno.n, 12);
        assert_eq!(pheno.trait_type, TraitType::Continuous);
        assert_eq!(pheno.y.nrows(), 12);
        assert_eq!(pheno.y.ncols(), 2);
        assert_eq!(pheno.x.ncols(), 3); // intercept + age + sex
        for i in 0..12 {
            assert_eq!(pheno.y[(i, 0)], 20.0 + i as f64);
            assert_eq!(pheno.y[(i, 1)], 150.0 + 3.0 * i as f64);
            assert_eq!(pheno.x[(i, 0)], 1.0);
            assert_eq!(pheno.x[(i, 1)], 30.0 + i as f64);
            assert_eq!(pheno.x[(i, 2)], (i % 2) as f64);
        }
        assert!(pheno.pheno_mask.iter().all(|v| *v));
    }

    /// Complete-cases semantics: a row with NULL in any trait is excluded
    /// by the SQL WHERE, so `pheno_mask` is false at that VCF index and
    /// `n` is reduced by one. The surviving rows must stay dense.
    #[test]
    fn load_phenotype_two_traits_drops_incomplete_cases() {
        let dir = tempfile::tempdir().unwrap();
        // s6 has a missing HEIGHT — blank field in the TSV becomes NULL
        // after register_csv, which the `IS NOT NULL AND` filter drops.
        let rows: Vec<Vec<String>> = (1..=12)
            .map(|i| {
                let id = format!("s{i}");
                let bmi = 20.0 + i as f64;
                let height = if i == 6 {
                    String::new()
                } else {
                    format!("{}", 150.0 + 3.0 * i as f64)
                };
                let age = 30.0 + i as f64;
                let sex = (i % 2) as f64;
                vec![
                    id,
                    format!("{bmi}"),
                    height,
                    format!("{age}"),
                    format!("{sex}"),
                ]
            })
            .collect();
        let row_refs: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(|s| s.as_str()).collect())
            .collect();
        let pheno_path = write_pheno_tsv(
            &dir,
            &["IID", "BMI", "HEIGHT", "age", "sex"],
            &row_refs,
        );
        let sample_ids: Vec<String> = (1..=12).map(|i| format!("s{i}")).collect();
        let geno = fake_geno(
            &sample_ids
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
        );

        let pheno = load_phenotype(
            &test_engine(),
            &pheno_path,
            &["age".into(), "sex".into()],
            &geno,
            &["BMI".into(), "HEIGHT".into()],
            None,
            5,
            42,
            &HashMap::new(),
            &SilentOut,
        )
        .expect("incomplete-case drop should not be a hard error");

        assert_eq!(pheno.n, 11, "sample with NULL HEIGHT must be excluded");
        assert_eq!(pheno.y.nrows(), 11);
        assert_eq!(pheno.y.ncols(), 2);
        assert_eq!(pheno.pheno_mask.len(), 12);
        assert!(
            !pheno.pheno_mask[5],
            "s6 (VCF index 5) should be masked out after complete-cases filter"
        );
        assert!(
            pheno.pheno_mask.iter().enumerate().all(|(i, &v)| v == (i != 5))
        );
        for i in 0..11 {
            assert!(pheno.y[(i, 0)].is_finite());
            assert!(pheno.y[(i, 1)].is_finite());
        }
    }

    /// Joint multi-trait STAAR requires every --trait-name to share the
    /// same family. Mixing a binary trait and a continuous trait in the
    /// same run must reject with an Input error naming both columns.
    #[test]
    fn load_phenotype_rejects_mixed_family_in_multi_mode() {
        let dir = tempfile::tempdir().unwrap();
        let rows: Vec<Vec<String>> = (1..=12)
            .map(|i| {
                let id = format!("s{i}");
                let bmi = (i % 2) as f64;
                let height = 150.0 + 3.0 * i as f64;
                let age = 30.0 + i as f64;
                let sex = (i % 2) as f64;
                vec![
                    id,
                    format!("{bmi}"),
                    format!("{height}"),
                    format!("{age}"),
                    format!("{sex}"),
                ]
            })
            .collect();
        let row_refs: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(|s| s.as_str()).collect())
            .collect();
        let pheno_path = write_pheno_tsv(
            &dir,
            &["IID", "BMI", "HEIGHT", "age", "sex"],
            &row_refs,
        );
        let sample_ids: Vec<String> = (1..=12).map(|i| format!("s{i}")).collect();
        let geno = fake_geno(
            &sample_ids
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
        );

        let result = load_phenotype(
            &test_engine(),
            &pheno_path,
            &["age".into(), "sex".into()],
            &geno,
            &["BMI".into(), "HEIGHT".into()],
            None,
            5,
            42,
            &HashMap::new(),
            &SilentOut,
        );
        match result {
            Ok(_) => panic!("mixed binary + continuous family must be rejected"),
            Err(CohortError::Input(msg)) => {
                assert!(msg.contains("BMI"), "error must name BMI: {msg}");
                assert!(msg.contains("HEIGHT"), "error must name HEIGHT: {msg}");
                assert!(
                    msg.contains("family"),
                    "error must explain the shared-family rule: {msg}"
                );
            }
            Err(other) => panic!("expected Input error, got {other:?}"),
        }
    }

    /// Multi-trait all-binary must now load successfully and produce a
    /// `(n, k)` y matrix with binary-coded phenotypes.
    #[test]
    fn load_phenotype_accepts_all_binary_multi() {
        let dir = tempfile::tempdir().unwrap();
        // Two binary traits: case1 and case2, each balanced 0/1.
        let rows: Vec<Vec<String>> = (1..=12)
            .map(|i| {
                let id = format!("s{i}");
                let case1 = (i % 2) as f64;
                let case2 = ((i + 1) % 2) as f64;
                let age = 30.0 + i as f64;
                let sex = (i % 2) as f64;
                vec![
                    id,
                    format!("{case1}"),
                    format!("{case2}"),
                    format!("{age}"),
                    format!("{sex}"),
                ]
            })
            .collect();
        let row_refs: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(|s| s.as_str()).collect())
            .collect();
        let pheno_path = write_pheno_tsv(
            &dir,
            &["IID", "case1", "case2", "age", "sex"],
            &row_refs,
        );
        let sample_ids: Vec<String> = (1..=12).map(|i| format!("s{i}")).collect();
        let geno = fake_geno(
            &sample_ids
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
        );

        let pheno = load_phenotype(
            &test_engine(),
            &pheno_path,
            &["age".into(), "sex".into()],
            &geno,
            &["case1".into(), "case2".into()],
            None,
            5,
            42,
            &HashMap::new(),
            &SilentOut,
        )
        .expect("all-binary multi must load");
        assert_eq!(pheno.trait_type, staar::TraitType::Binary);
        assert_eq!(pheno.y.ncols(), 2);
        assert_eq!(pheno.y.nrows(), 12);
        for i in 0..12 {
            assert!(pheno.y[(i, 0)] == 0.0 || pheno.y[(i, 0)] == 1.0);
            assert!(pheno.y[(i, 1)] == 0.0 || pheno.y[(i, 1)] == 1.0);
        }
    }

    /// Regression: the k=1 path must stay identical. A single continuous
    /// trait loads as a `(n, 1)` Y with the same values the old
    /// `&str`-taking loader would have produced.
    #[test]
    fn load_phenotype_single_trait_stays_single_column() {
        let dir = tempfile::tempdir().unwrap();
        let sample_ids: Vec<String> = (1..=12).map(|i| format!("s{i}")).collect();
        let rows: Vec<Vec<String>> = sample_ids
            .iter()
            .enumerate()
            .map(|(i, id)| {
                vec![
                    id.clone(),
                    format!("{}", 20.0 + i as f64),
                    format!("{}", 30.0 + i as f64),
                ]
            })
            .collect();
        let row_refs: Vec<Vec<&str>> = rows
            .iter()
            .map(|r| r.iter().map(|s| s.as_str()).collect())
            .collect();
        let pheno_path = write_pheno_tsv(&dir, &["IID", "BMI", "age"], &row_refs);
        let geno = fake_geno(
            &sample_ids
                .iter()
                .map(|s| s.as_str())
                .collect::<Vec<_>>(),
        );

        let pheno = load_phenotype(
            &test_engine(),
            &pheno_path,
            &["age".into()],
            &geno,
            &["BMI".into()],
            None,
            5,
            42,
            &HashMap::new(),
            &SilentOut,
        )
        .expect("single-trait load should succeed");

        assert_eq!(pheno.n, 12);
        assert_eq!(pheno.y.nrows(), 12);
        assert_eq!(pheno.y.ncols(), 1);
        assert_eq!(pheno.trait_type, TraitType::Continuous);
        for i in 0..12 {
            assert_eq!(pheno.y[(i, 0)], 20.0 + i as f64);
        }
    }
}
