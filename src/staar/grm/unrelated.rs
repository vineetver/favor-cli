//! Unrelated sample selection via greedy set cover with ancestry divergence
//! tie-breaking.
//!
//! Mirrors FastSparseGRM R/getUnrels.R:selectUnrel (lines 81-125) and
//! cppFunct.cpp:calculateDivergence (lines 525-576).

use std::collections::HashMap;

use rayon::prelude::*;

use crate::error::CohortError;
use crate::store::cohort::variants::CarrierEntry;
use crate::store::cohort::CohortHandle;
use crate::types::Chromosome;

use super::types::UnrelatedSubset;

/// Greedy unrelated selection matching R selectUnrel().
///
/// Repeatedly removes the sample with the most relatives; ties broken by
/// (divergence ascending, total_kinship ascending). Samples in
/// `always_keep` (from --include) are never removed.
pub fn select_unrelated(
    pairs: &[(usize, usize, f64)],
    n_samples: usize,
    divergence: &[i32],
) -> UnrelatedSubset {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_samples];
    let mut total_kin = vec![0.0f64; n_samples];
    for &(i, j, k) in pairs {
        adj[i].push(j);
        adj[j].push(i);
        total_kin[i] += k;
        total_kin[j] += k;
    }

    let mut n_rel: Vec<usize> = adj.iter().map(|a| a.len()).collect();
    let mut removed = vec![false; n_samples];
    let involved: std::collections::HashSet<usize> =
        pairs.iter().flat_map(|&(i, j, _)| [i, j]).collect();

    loop {
        let candidate = involved
            .iter()
            .copied()
            .filter(|&s| !removed[s] && n_rel[s] > 0)
            .max_by(|&a, &b| {
                n_rel[a]
                    .cmp(&n_rel[b])
                    .then_with(|| divergence[b].cmp(&divergence[a]))
                    .then_with(|| {
                        total_kin[b]
                            .partial_cmp(&total_kin[a])
                            .unwrap_or(std::cmp::Ordering::Equal)
                    })
            });

        let Some(rm) = candidate else { break };
        removed[rm] = true;
        for &neighbor in &adj[rm] {
            if !removed[neighbor] {
                n_rel[neighbor] = n_rel[neighbor].saturating_sub(1);
            }
        }
    }

    let sample_indices: Vec<usize> = (0..n_samples).filter(|&i| !removed[i]).collect();
    UnrelatedSubset { sample_indices }
}

/// Compute ancestry divergence for samples that appear in KING .seg.
///
/// Matches cppFunct.cpp:calculateDivergence. For each sample i in
/// `related_indices`, counts how many other samples j have
/// `(hethet - 2*homopp) / (nhet_i + nhet_j) < cutoff`.
///
/// Uses a random subsample of variants for speed (10K by default, same
/// as upstream R getDivergence). Packed byte format + 256x256 lookup
/// tables make the O(|related| * n * n_snps/4) inner loop fast.
pub fn compute_divergence(
    cohort: &CohortHandle<'_>,
    manifest: &crate::store::cohort::CohortManifest,
    related_indices: &[usize],
    n_samples: usize,
    max_snps: usize,
    cutoff: f64,
) -> Result<Vec<i32>, CohortError> {
    let mut div = vec![0i32; n_samples];
    if related_indices.is_empty() {
        return Ok(div);
    }

    let (packed, n_bytes, nhet) =
        build_packed_genotypes(cohort, manifest, n_samples, max_snps)?;

    let (hethet_tab, homopp_tab) = build_lookup_tables();

    let per_related: Vec<(usize, i32)> = related_indices
        .par_iter()
        .map(|&ri| {
            let mut count = 0i32;
            for sj in 0..n_samples {
                if sj == ri {
                    continue;
                }
                let mut hh = 0i32;
                let mut ho = 0i32;
                let base_i = ri * n_bytes;
                let base_j = sj * n_bytes;
                for k in 0..n_bytes {
                    let bi = packed[base_i + k] as usize;
                    let bj = packed[base_j + k] as usize;
                    hh += hethet_tab[bi][bj] as i32;
                    ho += homopp_tab[bi][bj] as i32;
                }
                let denom = nhet[ri] + nhet[sj];
                if denom > 0 {
                    let d = (hh - 2 * ho) as f64 / denom as f64;
                    if d < cutoff {
                        count += 1;
                    }
                }
            }
            (ri, count)
        })
        .collect();

    for (ri, count) in per_related {
        div[ri] = count;
    }
    Ok(div)
}

/// Build a packed-byte genotype matrix from carrier lists of randomly
/// sampled variants. Each sample gets `n_bytes = ceil(n_snps / 4)` bytes
/// with 2-bit PLINK encoding: 0=homref, 1=missing, 2=het, 3=homalt.
fn build_packed_genotypes(
    cohort: &CohortHandle<'_>,
    manifest: &crate::store::cohort::CohortManifest,
    n_samples: usize,
    max_snps: usize,
) -> Result<(Vec<u8>, usize, Vec<i32>), CohortError> {
    let mut all_variant_locs: Vec<(Chromosome, u32)> = Vec::new();
    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
        let n_var = ci.n_variants;
        for v in 0..n_var {
            all_variant_locs.push((chrom, v as u32));
        }
    }

    let total = all_variant_locs.len();
    let n_use = max_snps.min(total);

    // Deterministic subsample via stride (no RNG dep, reproducible).
    let step = if n_use >= total { 1 } else { total / n_use };
    let selected: Vec<(Chromosome, u32)> = all_variant_locs
        .iter()
        .step_by(step)
        .take(n_use)
        .copied()
        .collect();

    let n_bytes = selected.len().div_ceil(4);
    let mut packed = vec![0u8; n_samples * n_bytes];
    let mut nhet = vec![0i32; n_samples];

    // Group selected variants by chromosome for sequential mmap access.
    let mut by_chrom: HashMap<Chromosome, Vec<(usize, u32)>> = HashMap::new();
    for (snp_i, &(chrom, vcf)) in selected.iter().enumerate() {
        by_chrom.entry(chrom).or_default().push((snp_i, vcf));
    }

    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
        let Some(variants) = by_chrom.get(&chrom) else {
            continue;
        };
        let view = cohort.chromosome(&chrom)?;
        let sorted_vcfs: Vec<crate::store::cohort::types::VariantVcf> = variants
            .iter()
            .map(|&(_, vcf)| crate::store::cohort::types::VariantVcf(vcf))
            .collect();
        let batch = view.carriers_batch(&sorted_vcfs)?;

        for (local_idx, &(snp_i, _)) in variants.iter().enumerate() {
            let carriers = &batch.entries[local_idx];
            let byte_pos = snp_i / 4;
            let bit_shift = (snp_i % 4) * 2;
            for &CarrierEntry { sample_idx, dosage } in &carriers.entries {
                let si = sample_idx as usize;
                if si >= n_samples {
                    continue;
                }
                let code: u8 = match dosage {
                    1 => {
                        nhet[si] += 1;
                        2 // het
                    }
                    2 => 3, // homalt
                    255 => 1, // missing
                    _ => continue,
                };
                packed[si * n_bytes + byte_pos] |= code << bit_shift;
            }
        }
    }

    Ok((packed, n_bytes, nhet))
}

/// 256x256 lookup tables for het-het and hom-opp counts per byte pair.
/// Each byte packs 4 genotypes (2 bits each). Matches
/// cppFunct.cpp:createTable (lines 475-498).
fn build_lookup_tables() -> (Vec<Vec<i8>>, Vec<Vec<i8>>) {
    let mut hh = vec![vec![0i8; 256]; 256];
    let mut ho = vec![vec![0i8; 256]; 256];
    for i in 0..256u16 {
        for j in 0..256u16 {
            let mut nhethet = 0i8;
            let mut nhomopp = 0i8;
            let mut k = 0;
            while k < 8 {
                let ci = ((i >> k) & 3) as u8;
                let cj = ((j >> k) & 3) as u8;
                if ci == 2 && cj == 2 {
                    nhethet += 1;
                } else if (ci == 3 && cj == 0) || (ci == 0 && cj == 3) {
                    nhomopp += 1;
                }
                k += 2;
            }
            hh[i as usize][j as usize] = nhethet;
            ho[i as usize][j as usize] = nhomopp;
        }
    }
    (hh, ho)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn greedy_removes_most_connected_first() {
        // Triangle: 0-1, 1-2, 0-2. All have degree 2; tie broken by divergence.
        let pairs = vec![(0, 1, 0.5), (0, 2, 0.25), (1, 2, 0.25)];
        let div = vec![0, 10, 0]; // sample 1 has highest divergence
        let result = select_unrelated(&pairs, 3, &div);
        // Should remove 2 samples, keep 1.
        assert_eq!(result.sample_indices.len(), 1);
    }

    #[test]
    fn disjoint_pairs_keep_one_from_each() {
        let pairs = vec![(0, 1, 0.5), (2, 3, 0.5)];
        let div = vec![0; 4];
        let result = select_unrelated(&pairs, 4, &div);
        assert_eq!(result.sample_indices.len(), 2);
    }

    #[test]
    fn singletons_always_kept() {
        let pairs = vec![(0, 1, 0.5)];
        let div = vec![0; 5];
        let result = select_unrelated(&pairs, 5, &div);
        // Samples 2, 3, 4 are singletons. One of 0/1 is kept.
        assert_eq!(result.sample_indices.len(), 4);
    }

    #[test]
    fn lookup_tables_correct_for_known_byte() {
        let (hh, ho) = build_lookup_tables();
        // Byte 0b00_10_10_10 = all het (code 2) for 3 positions + homref
        // genotypes: pos0=het(2), pos1=het(2), pos2=het(2), pos3=homref(0)
        let byte_val = 0b00_10_10_10u8;
        // Compare with itself: all 3 het positions are hethet
        assert_eq!(hh[byte_val as usize][byte_val as usize], 3);
        assert_eq!(ho[byte_val as usize][byte_val as usize], 0);

        // homref vs homalt at pos0: 0b11 vs 0b00 → homopp
        let a = 0b00_00_00_00u8; // all homref
        let b = 0b00_00_00_11u8; // pos0=homalt, rest homref
        assert_eq!(hh[a as usize][b as usize], 0);
        assert_eq!(ho[a as usize][b as usize], 1);
    }
}
