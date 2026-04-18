//! Randomized PCA and carrier-indexed genotype-matrix operations.
//!
//! Implements per-chromosome G*v (postmultiply) and G'*v (premultiply)
//! on the sparse carrier representation in `sparse_g.bin`, matching
//! FastSparseGRM cppFunct.cpp:postmultiply (lines 252-283) and
//! premultiply (lines 331-366). The randomized SVD follows
//! runPCA.R:drpca (lines 2-78).
//!
//! All operations accumulate per-chromosome so the caller controls
//! memory: one ChromosomeView open at a time, sum across chromosomes.

use faer::Mat;

use crate::error::CohortError;
use crate::store::cohort::variants::CarrierEntry;
use crate::store::cohort::{ChromosomeView, CohortHandle, CohortManifest};
use crate::types::Chromosome;

use super::types::PcaScores;

/// Per-variant allele frequency and inverse standard deviation for a
/// subset of samples. `mu[v] = allele_count / (2 * n_nonmissing)`,
/// `inv_sd[v] = 1 / sqrt(2 * mu * (1 - mu))`. Variants with mu=0 or
/// mu=1 get inv_sd=0 so they contribute nothing to G*v.
pub struct VariantStats {
    pub mu: Vec<f64>,
    pub inv_sd: Vec<f64>,
}

/// Compute allele frequencies from a subset of samples on one chromosome.
/// Walks each variant's carrier list once.
pub fn allele_freq_chrom(
    view: &ChromosomeView<'_>,
    sample_set: &[bool],
) -> Result<VariantStats, CohortError> {
    let index = view.index()?;
    let n_variants = index.len();
    let n_subset: usize = sample_set.iter().filter(|&&b| b).count();

    let mut mu = vec![0.0f64; n_variants];
    let mut inv_sd = vec![0.0f64; n_variants];

    for v in 0..n_variants {
        let carriers = view.sparse_g()?.load_variant(v as u32);
        let mut allele_count = 0u64;
        let mut n_missing = 0u64;
        for &CarrierEntry { sample_idx, dosage } in &carriers.entries {
            let si = sample_idx as usize;
            if si >= sample_set.len() || !sample_set[si] {
                continue;
            }
            if dosage == 255 {
                n_missing += 1;
                continue;
            }
            allele_count += dosage as u64;
        }
        let n_obs = n_subset as u64 - n_missing;
        if n_obs == 0 {
            continue;
        }
        let p = allele_count as f64 / (2.0 * n_obs as f64);
        mu[v] = p;
        let var = 2.0 * p * (1.0 - p);
        if var > 1e-10 {
            inv_sd[v] = 1.0 / var.sqrt();
        }
    }

    Ok(VariantStats { mu, inv_sd })
}

/// G_chrom * v: (p_chrom × L) result.
///
/// For each variant, the centered-and-scaled genotype is
/// `(dosage - 2*mu) * inv_sd`. Non-carriers (dosage=0) contribute
/// `-2*mu*inv_sd * v[sample]` which aggregates as a constant shift per
/// variant: `-2*mu*inv_sd * sum(v[subset])`. Carrier contributions
/// deviate from this baseline by `(dosage - 0) * inv_sd * v[sample]`.
///
/// O(total_carriers_chrom × L) not O(n × p × L).
pub fn postmultiply_chrom(
    view: &ChromosomeView<'_>,
    v: &Mat<f64>,
    sample_set: &[bool],
    stats: &VariantStats,
) -> Result<Mat<f64>, CohortError> {
    let n_variants = stats.mu.len();
    let l = v.ncols();

    let mut col_sums = vec![0.0f64; l];
    for (si, &in_set) in sample_set.iter().enumerate() {
        if in_set {
            for c in 0..l {
                col_sums[c] += v[(si, c)];
            }
        }
    }

    let mut result = Mat::<f64>::zeros(n_variants, l);

    for vi in 0..n_variants {
        let mu_v = stats.mu[vi];
        let isd = stats.inv_sd[vi];
        if isd == 0.0 {
            continue;
        }
        let carriers = view.sparse_g()?.load_variant(vi as u32);

        let mut dosage_v_sum = vec![0.0f64; l];
        let mut missing_v_sum = vec![0.0f64; l];

        for &CarrierEntry { sample_idx, dosage } in &carriers.entries {
            let si = sample_idx as usize;
            if si >= sample_set.len() || !sample_set[si] {
                continue;
            }
            if dosage == 255 {
                for c in 0..l {
                    missing_v_sum[c] += v[(si, c)];
                }
                continue;
            }
            for c in 0..l {
                dosage_v_sum[c] += dosage as f64 * v[(si, c)];
            }
        }

        // result[vi, c] = sum_over_samples((g - 2*mu) * inv_sd * v[s, c])
        //               = inv_sd * (dosage_v_sum - 2*mu * (col_sums - missing_v_sum))
        // Missing samples are imputed to mean → contribute 0 to centered genotype.
        for c in 0..l {
            let obs_v_sum = col_sums[c] - missing_v_sum[c];
            result[(vi, c)] = isd * (dosage_v_sum[c] - 2.0 * mu_v * obs_v_sum);
        }
    }

    Ok(result)
}

/// G_chrom' * v: accumulates into `result` (n_samples × L).
///
/// For each variant, the centered-and-scaled contribution to sample s is
/// `(g_s - 2*mu) * inv_sd * v[snp]`. Non-carriers get `-2*mu*inv_sd*v[snp]`
/// as a constant shift applied to all samples in the set; carriers get an
/// additional `dosage * inv_sd * v[snp]`. Missing samples get nothing.
///
/// Precomputes `mu_ratio_sum[c] = sum_snps(2*mu*inv_sd*v[snp,c])` once,
/// then scatters carrier deviations per variant. O(total_carriers_chrom × L).
pub fn premultiply_chrom(
    view: &ChromosomeView<'_>,
    v: &Mat<f64>,
    sample_set: &[bool],
    stats: &VariantStats,
    result: &mut Mat<f64>,
) -> Result<(), CohortError> {
    let n_variants = stats.mu.len();
    let l = v.ncols();

    let mut mu_ratio_sum = vec![0.0f64; l];
    for vi in 0..n_variants {
        let factor = 2.0 * stats.mu[vi] * stats.inv_sd[vi];
        if factor == 0.0 {
            continue;
        }
        for c in 0..l {
            mu_ratio_sum[c] += factor * v[(vi, c)];
        }
    }

    // Baseline: every sample in set gets -mu_ratio_sum (non-carrier shift).
    for (si, &in_set) in sample_set.iter().enumerate() {
        if in_set {
            for c in 0..l {
                result[(si, c)] -= mu_ratio_sum[c];
            }
        }
    }

    // Carrier deviations: add dosage*inv_sd*v per carrier, and undo the
    // -2*mu*inv_sd*v baseline for missing samples (they should contribute 0).
    for vi in 0..n_variants {
        let isd = stats.inv_sd[vi];
        if isd == 0.0 {
            continue;
        }
        let mu_shift = 2.0 * stats.mu[vi] * isd;
        let carriers = view.sparse_g()?.load_variant(vi as u32);
        for &CarrierEntry { sample_idx, dosage } in &carriers.entries {
            let si = sample_idx as usize;
            if si >= sample_set.len() || !sample_set[si] {
                continue;
            }
            if dosage == 255 {
                // Missing: undo the baseline shift (this sample should get 0).
                for c in 0..l {
                    result[(si, c)] += mu_shift * v[(vi, c)];
                }
                continue;
            }
            let d_isd = dosage as f64 * isd;
            for c in 0..l {
                result[(si, c)] += d_isd * v[(vi, c)];
            }
        }
    }

    Ok(())
}

/// Full randomized PCA across all chromosomes.
///
/// Block Krylov iteration matching runPCA.R:drpca (lines 2-78).
/// Each iteration accumulates G*h into a SNP-space Krylov matrix H,
/// then the SVD of H gives a high-quality eigenspace approximation.
/// Per-chromosome accumulation keeps peak memory to one chromosome.
pub fn randomized_pca(
    cohort: &CohortHandle<'_>,
    manifest: &CohortManifest,
    unrelated_mask: &[bool],
    _all_sample_mask: &[bool],
    n_pcs: usize,
    n_iter: usize,
) -> Result<PcaScores, CohortError> {
    let n_unrel: usize = unrelated_mask.iter().filter(|&&b| b).count();
    let l = 2 * n_pcs;

    let mut chrom_stats: Vec<(Chromosome, VariantStats)> = Vec::new();
    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
        let view = cohort.chromosome(&chrom)?;
        let stats = allele_freq_chrom(&view, unrelated_mask)?;
        chrom_stats.push((chrom, stats));
    }

    let total_variants: usize = chrom_stats.iter().map(|(_, s)| s.mu.len()).sum();

    let mut rng = super::super::scang::Xorshift64::new(42);
    let mut h = Mat::<f64>::from_fn(n_unrel, l, |_, _| {
        super::super::scang::standard_normal(&mut rng)
    });

    let unrel_compact: Vec<Option<usize>> = {
        let mut map = Vec::with_capacity(unrelated_mask.len());
        let mut next = 0usize;
        for &b in unrelated_mask {
            if b {
                map.push(Some(next));
                next += 1;
            } else {
                map.push(None);
            }
        }
        map
    };

    // Krylov accumulation: H = [x_1 | x_2 | … | x_{n_iter}] where
    // x_k = G * h_k (SNP-space). Matches drpca line 33: H<-cbind(H,x).
    let krylov_cols = n_iter * l;
    let mut krylov = Mat::<f64>::zeros(total_variants, krylov_cols);

    for iter in 0..n_iter {
        let h_full = expand_compact(&h, unrelated_mask, &unrel_compact);
        let mut x = Mat::<f64>::zeros(total_variants, l);
        let mut offset = 0usize;
        for (chrom, stats) in &chrom_stats {
            let view = cohort.chromosome(chrom)?;
            let x_chrom = postmultiply_chrom(&view, &h_full, unrelated_mask, stats)?;
            let p = stats.mu.len();
            for vi in 0..p {
                for c in 0..l {
                    x[(offset + vi, c)] = x_chrom[(vi, c)];
                }
            }
            offset += p;
        }

        let col_base = iter * l;
        for vi in 0..total_variants {
            for c in 0..l {
                krylov[(vi, col_base + c)] = x[(vi, c)];
            }
        }

        let mut h_new_full = Mat::<f64>::zeros(unrelated_mask.len(), l);
        offset = 0;
        for (chrom, stats) in &chrom_stats {
            let view = cohort.chromosome(chrom)?;
            let p = stats.mu.len();
            let x_slice = Mat::from_fn(p, l, |i, c| x[(offset + i, c)]);
            premultiply_chrom(&view, &x_slice, unrelated_mask, stats, &mut h_new_full)?;
            offset += p;
        }

        h = compact_rows(&h_new_full, unrelated_mask, &unrel_compact);

        for c in 0..l {
            let mut norm = 0.0f64;
            for i in 0..n_unrel {
                norm += h[(i, c)] * h[(i, c)];
            }
            let inv = if norm > 0.0 { 1.0 / norm.sqrt() } else { 0.0 };
            for i in 0..n_unrel {
                h[(i, c)] *= inv;
            }
        }
    }

    // SVD of accumulated Krylov matrix (SNP-space). Matches drpca
    // line 42: init<-svd(H,nv=0)$u.
    let svd_k = krylov.thin_svd().map_err(|e| {
        CohortError::Analysis(format!("PCA Krylov SVD failed: {e:?}"))
    })?;
    let u_snp = svd_k.U();
    let n_krylov_cols = u_snp.ncols();

    // T = G' * U_snp using ALL Krylov columns, then truncate via SVD.
    // Matches drpca lines 44-46: T<-premultiply(init,...); svd(T,nu=nd).
    let mut t_mat = Mat::<f64>::zeros(unrelated_mask.len(), n_krylov_cols);
    let mut offset = 0usize;
    for (chrom, stats) in &chrom_stats {
        let view = cohort.chromosome(chrom)?;
        let p = stats.mu.len();
        let u_chrom = Mat::from_fn(p, n_krylov_cols, |i, c| u_snp[(offset + i, c)]);
        premultiply_chrom(&view, &u_chrom, unrelated_mask, stats, &mut t_mat)?;
        offset += p;
    }
    let t_compact = compact_rows(&t_mat, unrelated_mask, &unrel_compact);

    let svd_t = t_compact.thin_svd().map_err(|e| {
        CohortError::Analysis(format!("PCA final SVD failed: {e:?}"))
    })?;
    let u_final = svd_t.U();
    let s_diag = svd_t.S();
    let s_col = s_diag.column_vector();
    let nd = n_pcs.min(u_final.ncols());

    // Singular values of T = G' * U_snp ≈ Σ_G. Eigenvalues = Σ_G².
    let eigenvalues: Vec<f64> = (0..nd).map(|c| s_col[c] * s_col[c]).collect();

    let mut unrel_to_global = vec![0usize; n_unrel];
    {
        let mut ci = 0usize;
        for (gi, &b) in unrelated_mask.iter().enumerate() {
            if b {
                unrel_to_global[ci] = gi;
                ci += 1;
            }
        }
    }

    let n_full = unrelated_mask.len();
    let mut scores = Mat::<f64>::zeros(n_full, nd);
    for ci in 0..n_unrel {
        let gi = unrel_to_global[ci];
        for c in 0..nd {
            scores[(gi, c)] = u_final[(ci, c)];
        }
    }

    // Related projection: v_res = G * u_final / d, then
    // scores_rel = G_rel' * v_res / d = G_rel' * G * u_final / d².
    let related_mask: Vec<bool> = unrelated_mask.iter().map(|&b| !b).collect();
    let n_related: usize = related_mask.iter().filter(|&&b| b).count();
    if n_related > 0 {
        let u_final_full = expand_compact(
            &Mat::from_fn(n_unrel, nd, |i, c| u_final[(i, c)]),
            unrelated_mask,
            &unrel_compact,
        );
        let mut rel_accum = Mat::<f64>::zeros(n_full, nd);
        for (chrom, stats) in &chrom_stats {
            let view = cohort.chromosome(chrom)?;
            let g_u = postmultiply_chrom(&view, &u_final_full, unrelated_mask, stats)?;
            premultiply_chrom(&view, &g_u, &related_mask, stats, &mut rel_accum)?;
        }
        for (gi, &is_rel) in related_mask.iter().enumerate() {
            if is_rel {
                for c in 0..nd {
                    let d = s_col[c];
                    if d.abs() > 1e-12 {
                        scores[(gi, c)] = rel_accum[(gi, c)] / (d * d);
                    }
                }
            }
        }
    }

    Ok(PcaScores { scores, eigenvalues })
}

fn expand_compact(
    compact: &Mat<f64>,
    mask: &[bool],
    compact_idx: &[Option<usize>],
) -> Mat<f64> {
    let n_full = mask.len();
    let l = compact.ncols();
    Mat::from_fn(n_full, l, |i, c| {
        compact_idx[i].map_or(0.0, |ci| compact[(ci, c)])
    })
}

fn compact_rows(
    full: &Mat<f64>,
    mask: &[bool],
    compact_idx: &[Option<usize>],
) -> Mat<f64> {
    let n_compact: usize = mask.iter().filter(|&&b| b).count();
    let l = full.ncols();
    let mut out = Mat::<f64>::zeros(n_compact, l);
    for (i, &in_set) in mask.iter().enumerate() {
        if in_set {
            if let Some(ci) = compact_idx[i] {
                for c in 0..l {
                    out[(ci, c)] = full[(i, c)];
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn lookup_tables_in_scang_module_accessible() {
        // Smoke test that the xorshift + standard_normal from scang are reachable.
        let mut rng = crate::staar::scang::Xorshift64::new(1);
        let v = crate::staar::scang::standard_normal(&mut rng);
        assert!(v.is_finite());
    }
}
