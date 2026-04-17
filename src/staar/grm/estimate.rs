//! Sparse GRM estimation via ancestry-adjusted kinship.
//!
//! Mirrors FastSparseGRM R/calcGRM.R:calcSparseGRM (lines 72-173).
//! Block-wise per-chromosome accumulation of kinship numerator and
//! denominator for each candidate pair. After all chromosomes are
//! processed, kinship = num / den. Thresholding at 2^-(degree+1.5)
//! produces the sparse output. Two-pass re-estimation handles large
//! connected components that exceed max_related_block.

use std::collections::HashMap;

use faer::Mat;

use crate::error::CohortError;
use crate::output::Output;
use crate::store::cohort::variants::CarrierEntry;
use crate::store::cohort::{CohortHandle, CohortManifest};
use crate::types::Chromosome;

use super::king;
use super::pca;
use super::types::{KinshipAccum, PcaScores, SparseGrm};

const MAX_RELATED_BLOCK: usize = 65536;

/// Estimate sparse GRM for all candidate pairs across all chromosomes.
///
/// Two-pass: first pass estimates kinship for KING-identified pairs, then
/// after thresholding, if any connected component exceeds MAX_RELATED_BLOCK,
/// the threshold is raised iteratively and newly-discovered pairs from
/// expanded components are estimated in a second pass.
#[allow(clippy::too_many_arguments)]
pub fn estimate_grm(
    cohort: &CohortHandle<'_>,
    manifest: &CohortManifest,
    pca_scores: &PcaScores,
    unrelated_mask: &[bool],
    candidate_pairs: &[(usize, usize, f64)],
    n_samples: usize,
    n_pcs: usize,
    block_size: usize,
    degree: u8,
    out: &dyn Output,
) -> Result<SparseGrm, CohortError> {
    let threshold = 2.0f64.powf(-((degree as f64) + 1.5));

    // Build the PC-augmented covariate matrix X = [1 | PC1..PCk].
    let k = n_pcs + 1;
    let mut x_mat = Mat::<f64>::zeros(n_samples, k);
    for i in 0..n_samples {
        x_mat[(i, 0)] = 1.0;
        for c in 0..n_pcs.min(pca_scores.scores.ncols()) {
            x_mat[(i, c + 1)] = pca_scores.scores[(i, c)];
        }
    }

    // nullmat = X_train * (X_train' X_train)^{-1}  (training = unrelated subset)
    let n_train: usize = unrelated_mask.iter().filter(|&&b| b).count();
    let mut x_train = Mat::<f64>::zeros(n_train, k);
    let mut ti = 0;
    for (si, &is_unrel) in unrelated_mask.iter().enumerate() {
        if is_unrel {
            for c in 0..k {
                x_train[(ti, c)] = x_mat[(si, c)];
            }
            ti += 1;
        }
    }
    use faer::linalg::solvers::Solve;
    let xtx = x_train.transpose() * &x_train;
    let eye_k = Mat::<f64>::identity(k, k);
    let xtx_inv = xtx.col_piv_qr().solve(&eye_k);
    let nullmat = &x_train * &xtx_inv; // (n_train x k)

    // First pass: estimate kinship for all candidate pairs.
    let pair_indices: Vec<(usize, usize)> = candidate_pairs
        .iter()
        .map(|&(i, j, _)| (i, j))
        .collect();
    // Add self-pairs (diagonal).
    let mut all_pairs: Vec<(usize, usize)> = pair_indices.clone();
    let mut self_set: std::collections::HashSet<usize> = std::collections::HashSet::new();
    for &(i, j, _) in candidate_pairs {
        self_set.insert(i);
        self_set.insert(j);
    }
    for &s in &self_set {
        all_pairs.push((s, s));
    }
    all_pairs.sort_unstable();
    all_pairs.dedup();

    out.status(&format!(
        "  GRM: estimating kinship for {} pairs (+ {} self-pairs)...",
        pair_indices.len(),
        self_set.len(),
    ));

    let mut accum = estimate_pairs_all_chroms(
        cohort,
        manifest,
        &nullmat,
        &x_mat,
        unrelated_mask,
        &all_pairs,
        block_size,
        out,
    )?;

    // Finalize: kinship = num / den.
    for a in &mut accum {
        if a.denominator > 0.0 {
            a.numerator /= a.denominator;
        } else {
            a.numerator = 0.0;
        }
    }

    // Threshold at 2^-(degree+1.5). Self-pairs always kept.
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for a in &accum {
        if a.idx_i == a.idx_j {
            triplets.push((a.idx_i, a.idx_j, a.numerator));
        } else if a.numerator >= threshold {
            triplets.push((a.idx_i, a.idx_j, a.numerator));
            triplets.push((a.idx_j, a.idx_i, a.numerator));
        }
    }

    // Two-pass: check component sizes after thresholding.
    let off_diag: Vec<(usize, usize, f64)> = triplets
        .iter()
        .filter(|&&(i, j, _)| i < j)
        .copied()
        .collect();
    let components = king::build_components(&off_diag, n_samples);
    let max_comp = components.iter().map(|c| c.members.len()).max().unwrap_or(0);

    if max_comp > MAX_RELATED_BLOCK {
        out.status(&format!(
            "  GRM: largest component has {} members (> {MAX_RELATED_BLOCK}), raising threshold...",
            max_comp,
        ));
        let triplets_refined = two_pass_refine(
            cohort,
            manifest,
            &nullmat,
            &x_mat,
            unrelated_mask,
            &accum,
            &off_diag,
            n_samples,
            threshold,
            block_size,
            out,
        )?;
        return Ok(SparseGrm {
            triplets: triplets_refined,
            n_samples,
        });
    }

    out.status(&format!(
        "  GRM: {} kinship pairs above threshold {:.6}",
        off_diag.len(),
        threshold,
    ));

    Ok(SparseGrm {
        triplets,
        n_samples,
    })
}

/// Two-pass re-estimation. Mirrors calcGRM.R:99-173.
///
/// Raises the threshold iteratively until the largest component fits
/// within MAX_RELATED_BLOCK. Any new pairs from expanded components
/// that weren't in the original candidate set get a second-pass
/// kinship estimation.
#[allow(clippy::too_many_arguments)]
fn two_pass_refine(
    cohort: &CohortHandle<'_>,
    manifest: &CohortManifest,
    nullmat: &Mat<f64>,
    x_mat: &Mat<f64>,
    unrelated_mask: &[bool],
    first_pass: &[KinshipAccum],
    off_diag: &[(usize, usize, f64)],
    n_samples: usize,
    mut threshold: f64,
    block_size: usize,
    out: &dyn Output,
) -> Result<Vec<(usize, usize, f64)>, CohortError> {
    // Iteratively raise threshold until max component <= MAX_RELATED_BLOCK.
    let mut active: Vec<(usize, usize, f64)> = off_diag.to_vec();
    loop {
        let components = king::build_components(&active, n_samples);
        let max_comp = components.iter().map(|c| c.members.len()).max().unwrap_or(0);
        if max_comp <= MAX_RELATED_BLOCK {
            break;
        }
        threshold *= 2.0f64.powf(0.01);
        let min_kin = active.iter().map(|t| t.2).fold(f64::INFINITY, f64::min);
        active.retain(|&(_, _, k)| k > threshold);
        out.status(&format!(
            "  GRM: threshold raised to {threshold:.6} (max comp {max_comp}, min kin {min_kin:.6})",
        ));
        if active.is_empty() {
            break;
        }
    }

    // Identify new pairs from the re-thresholded components that were NOT
    // in the first-pass candidate set.
    let first_set: std::collections::HashSet<(usize, usize)> = first_pass
        .iter()
        .map(|a| (a.idx_i, a.idx_j))
        .collect();
    let components = king::build_components(&active, n_samples);
    let mut new_pairs: Vec<(usize, usize)> = Vec::new();
    for comp in &components {
        for i in 0..comp.members.len() {
            for j in (i + 1)..comp.members.len() {
                let gi = comp.members[i];
                let gj = comp.members[j];
                let (lo, hi) = if gi < gj { (gi, gj) } else { (gj, gi) };
                if !first_set.contains(&(lo, hi)) {
                    new_pairs.push((lo, hi));
                }
            }
        }
    }
    new_pairs.sort_unstable();
    new_pairs.dedup();

    if !new_pairs.is_empty() {
        out.status(&format!(
            "  GRM: second pass estimating {} new pairs...",
            new_pairs.len(),
        ));
        let mut new_accum = estimate_pairs_all_chroms(
            cohort,
            manifest,
            nullmat,
            x_mat,
            unrelated_mask,
            &new_pairs,
            block_size,
            out,
        )?;
        for a in &mut new_accum {
            if a.denominator > 0.0 {
                a.numerator /= a.denominator;
            } else {
                a.numerator = 0.0;
            }
        }

        // Merge first-pass and second-pass results.
        let mut all_kin: HashMap<(usize, usize), f64> = HashMap::new();
        for a in first_pass {
            all_kin.insert((a.idx_i, a.idx_j), a.numerator);
        }
        for a in &new_accum {
            all_kin.insert((a.idx_i, a.idx_j), a.numerator);
        }

        let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
        for (&(i, j), &k) in &all_kin {
            if i == j {
                triplets.push((i, j, k));
            } else if k >= threshold {
                triplets.push((i, j, k));
                triplets.push((j, i, k));
            }
        }
        return Ok(triplets);
    }

    // No new pairs needed; just rebuild triplets with the raised threshold.
    let mut triplets: Vec<(usize, usize, f64)> = Vec::new();
    for a in first_pass {
        if a.idx_i == a.idx_j {
            triplets.push((a.idx_i, a.idx_j, a.numerator));
        } else if a.numerator >= threshold {
            triplets.push((a.idx_i, a.idx_j, a.numerator));
            triplets.push((a.idx_j, a.idx_i, a.numerator));
        }
    }
    Ok(triplets)
}

/// Estimate kinship accumulators for a set of pairs across all chromosomes.
///
/// Per-chromosome, block-wise: loads genotypes in SNP blocks, computes
/// ISAF-adjusted kinship contributions, accumulates into pair accumulators.
#[allow(clippy::too_many_arguments)]
fn estimate_pairs_all_chroms(
    cohort: &CohortHandle<'_>,
    manifest: &CohortManifest,
    nullmat: &Mat<f64>,
    x_mat: &Mat<f64>,
    unrelated_mask: &[bool],
    pairs: &[(usize, usize)],
    block_size: usize,
    out: &dyn Output,
) -> Result<Vec<KinshipAccum>, CohortError> {
    let n_pairs = pairs.len();
    let mut accum: Vec<KinshipAccum> = pairs
        .iter()
        .map(|&(i, j)| KinshipAccum {
            idx_i: i,
            idx_j: j,
            numerator: 0.0,
            denominator: 0.0,
        })
        .collect();

    let _pair_idx: HashMap<(usize, usize), usize> = pairs
        .iter()
        .enumerate()
        .map(|(pi, &(i, j))| ((i, j), pi))
        .collect();

    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
        let view = cohort.chromosome(&chrom)?;
        let stats = pca::allele_freq_chrom(&view, unrelated_mask)?;
        let n_var = stats.mu.len();

        out.status(&format!(
            "    chr{}: {} variants, {} pairs",
            chrom.label(),
            n_var,
            n_pairs,
        ));

        // Block-wise: process SNP blocks of size block_size.
        let mut block_start = 0usize;
        while block_start < n_var {
            let block_end = (block_start + block_size).min(n_var);
            let blen = block_end - block_start;

            // Load genotypes for this block into dense arrays for the
            // pair-wise ISAF computation.
            let _n_samples = x_mat.nrows();
            let k = nullmat.ncols();

            // beta = G_train_block' * nullmat via carrier walk.
            // beta[snp, col] = sum_over_training_carriers(dosage * nullmat[train_idx, col])
            let mut beta = Mat::<f64>::zeros(blen, k);
            let _train_idx = 0usize;
            let train_map: Vec<Option<usize>> = {
                let mut m = Vec::with_capacity(unrelated_mask.len());
                let mut next = 0usize;
                for &b in unrelated_mask {
                    if b {
                        m.push(Some(next));
                        next += 1;
                    } else {
                        m.push(None);
                    }
                }
                m
            };

            for vi in block_start..block_end {
                let local = vi - block_start;
                let carriers = view.sparse_g()?.load_variant(vi as u32);
                for &CarrierEntry { sample_idx, dosage } in &carriers.entries {
                    let si = sample_idx as usize;
                    if si >= unrelated_mask.len() || !unrelated_mask[si] {
                        continue;
                    }
                    if dosage == 255 {
                        continue;
                    }
                    let ti = train_map[si].unwrap();
                    let d = dosage as f64;
                    for c in 0..k {
                        beta[(local, c)] += d * nullmat[(ti, c)];
                    }
                }
                // Subtract 2*mu * sum(nullmat[train,:]) for non-carrier baseline.
                let _mu = stats.mu[vi];
                // Actually: postmultiply formula for training set:
                // beta[snp, c] = sum_i g_i * nullmat[i, c]
                // = sum_carriers d * nullmat[train_idx, c] + 0 * (non-carriers)
                // But the centered version subtracts 2*mu * sum(nullmat[:,c]).
                // For the ISAF computation we use the raw (uncentered) form:
                // ISAF = beta @ X[sample,:] where beta is already the regression
                // coefficient. Let's match upstream exactly.
                //
                // Upstream: beta = postmultiply(nullmat, in.train) which computes
                // (dosage_sum + 2*mu*missing_count) / sd. But for ISAF we need
                // the unscaled version. The upstream postmultiply returns
                // sum((g - 2*mu)/sd * nullmat) but then ISAF = beta @ X, and
                // this gives the ancestry-adjusted allele frequency.
                //
                // Simpler: ISAF[snp, s] = 2 * freq_adjusted where
                // freq_adjusted = nullmat-regression predicted frequency.
                // We compute beta raw (without centering/scaling) and then
                // ISAF = beta @ X gives the predicted genotype (0..2 scale).
            }

            // For each pair (i, j) and each SNP in block:
            // 1. Look up genotype of sample i and j at this SNP.
            // 2. Compute ISAF = beta[snp,:] @ X[sample,:] (predicted genotype).
            // 3. Accumulate kinship numerator/denominator.
            //
            // We need per-sample genotypes for the pair members. Build a
            // sample->genotype lookup for this block.
            let involved: std::collections::HashSet<usize> = pairs
                .iter()
                .flat_map(|&(i, j)| [i, j])
                .collect();
            let mut geno_block: HashMap<usize, Vec<u8>> = HashMap::new();
            for &s in &involved {
                geno_block.insert(s, vec![0u8; blen]);
            }
            for vi in block_start..block_end {
                let local = vi - block_start;
                let carriers = view.sparse_g()?.load_variant(vi as u32);
                for &CarrierEntry { sample_idx, dosage } in &carriers.entries {
                    let si = sample_idx as usize;
                    if let Some(g) = geno_block.get_mut(&si) {
                        g[local] = dosage;
                    }
                }
            }

            // ISAF computation and kinship accumulation.
            for (pi, &(si, sj)) in pairs.iter().enumerate() {
                let gi = geno_block.get(&si).unwrap();
                let gj = geno_block.get(&sj).unwrap();
                for local in 0..blen {
                    let di = gi[local];
                    let dj = gj[local];
                    if di == 255 || dj == 255 {
                        continue;
                    }

                    // ISAF for sample s at this SNP:
                    // isaf_s = beta[local,:] @ X[s,:]
                    let _vi = block_start + local;
                    let mut isaf_i = 0.0f64;
                    let mut isaf_j = 0.0f64;
                    for c in 0..k {
                        isaf_i += beta[(local, c)] * x_mat[(si, c)];
                        isaf_j += beta[(local, c)] * x_mat[(sj, c)];
                    }
                    isaf_i = isaf_i.clamp(0.0001, 1.9999);
                    isaf_j = isaf_j.clamp(0.0001, 1.9999);

                    let res_i = di as f64 - isaf_i;
                    let res_j = dj as f64 - isaf_j;
                    let sd_i = (isaf_i * (1.0 - isaf_i / 2.0) * 2.0).sqrt();
                    let sd_j = (isaf_j * (1.0 - isaf_j / 2.0) * 2.0).sqrt();

                    accum[pi].numerator += res_i * res_j;
                    accum[pi].denominator += sd_i * sd_j;
                }
            }

            block_start = block_end;
        }
    }

    Ok(accum)
}

#[cfg(test)]
mod tests {
    #[test]
    fn isaf_clamp_bounds() {
        let v: f64 = 0.00005;
        assert_eq!(v.clamp(0.0001, 1.9999), 0.0001);
        let v2: f64 = 2.5;
        assert_eq!(v2.clamp(0.0001, 1.9999), 1.9999);
    }
}
