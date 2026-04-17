//! `favor grm` subcommand: FastSparseGRM pipeline.

use std::path::PathBuf;

use serde_json::json;

use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::staar::grm::{cache, estimate, king, pca, unrelated};
use crate::staar::grm::types::GrmArtifact;
use crate::store::cohort::CohortId;
use crate::store::ids::CacheKey;

pub struct GrmArgs {
    pub cohort: String,
    pub king_seg: PathBuf,
    pub degree: u8,
    pub n_pcs: usize,
    pub block_size: usize,
    pub output: Option<PathBuf>,
}

pub fn run(
    engine: &Engine,
    args: GrmArgs,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    if !args.king_seg.exists() {
        return Err(CohortError::Input(format!(
            "KING .seg file not found: '{}'",
            args.king_seg.display()
        )));
    }
    if args.degree == 0 || args.degree > 10 {
        return Err(CohortError::Input(format!(
            "--degree must be 1..10, got {}",
            args.degree
        )));
    }

    let cohort_id = CohortId::new(args.cohort.trim().to_string());
    let cohort = engine.cohort(&cohort_id);
    let store_result = cohort.load()?;
    let manifest = &store_result.manifest;

    let fp = cache::fingerprint(
        &manifest.key,
        &args.king_seg,
        args.degree,
        args.n_pcs,
    )?;
    let cache_dir = args.output.clone().unwrap_or_else(|| {
        engine
            .store()
            .layout()
            .grm_cache_dir(&cohort_id, &CacheKey::new(&fp))
    });

    if cache::probe(&cache_dir) {
        out.status(&format!(
            "GRM cache hit at {}",
            cache_dir.display()
        ));
        out.result_json(&json!({
            "status": "cache_hit",
            "cache_dir": cache_dir.to_string_lossy(),
            "kinship": cache::grm_tsv_path(&cache_dir).to_string_lossy(),
            "pca": cache::pca_tsv_path(&cache_dir).to_string_lossy(),
        }));
        return Ok(());
    }

    if dry_run {
        out.result_json(&json!({
            "command": "grm",
            "cohort_id": cohort_id.as_str(),
            "king_seg": args.king_seg.to_string_lossy(),
            "degree": args.degree,
            "n_pcs": args.n_pcs,
            "n_samples": manifest.n_samples,
            "n_variants": manifest.n_variants,
            "output_dir": cache_dir.to_string_lossy(),
        }));
        return Ok(());
    }

    let sample_ids = store_result.geno.sample_names.clone();
    let n_samples = sample_ids.len();
    out.status(&format!(
        "GRM: {} samples, {} variants across {} chromosomes",
        n_samples, manifest.n_variants, manifest.chromosomes.len()
    ));

    // 1. Parse KING .seg.
    out.status("  Parsing KING .seg file...");
    let seg_entries = king::parse_king_seg(&args.king_seg, args.degree)?;
    out.status(&format!("    {} related pairs after degree-{} filter", seg_entries.len(), args.degree));

    // Map to cohort indices.
    let king_ids: Vec<String> = sample_ids
        .iter()
        .map(|s| format!("{}_{}", s, s))
        .collect();
    let (candidate_pairs, _id_map) = king::map_to_cohort_indices(&seg_entries, &king_ids);
    out.status(&format!("    {} pairs mapped to cohort", candidate_pairs.len()));

    if candidate_pairs.is_empty() {
        return Err(CohortError::Input(
            "No KING pairs mapped to cohort samples. Check that KING sample IDs \
             match the VCF sample IDs (FID_IID format)."
                .into(),
        ));
    }

    // 2. Compute divergence + select unrelated.
    out.status("  Computing ancestry divergence...");
    let related_indices: Vec<usize> = {
        let mut s: std::collections::HashSet<usize> = std::collections::HashSet::new();
        for &(i, j, _) in &candidate_pairs {
            s.insert(i);
            s.insert(j);
        }
        let mut v: Vec<usize> = s.into_iter().collect();
        v.sort_unstable();
        v
    };
    let divergence = unrelated::compute_divergence(
        &cohort, manifest, &related_indices, n_samples, 10_000, -0.025,
    )?;

    out.status("  Selecting unrelated samples...");
    let unrel = unrelated::select_unrelated(&candidate_pairs, n_samples, &divergence);
    out.status(&format!("    {} unrelated samples selected", unrel.sample_indices.len()));

    let unrelated_mask: Vec<bool> = (0..n_samples)
        .map(|i| unrel.sample_indices.contains(&i))
        .collect();

    // 3. Randomized PCA.
    out.status(&format!("  Randomized PCA ({} components)...", args.n_pcs));
    let all_mask = vec![true; n_samples];
    let pca_scores = pca::randomized_pca(
        &cohort,
        manifest,
        &unrelated_mask,
        &all_mask,
        args.n_pcs,
        10,
    )?;
    out.status(&format!(
        "    PCA complete: {} eigenvalues",
        pca_scores.eigenvalues.len()
    ));

    // 4. Estimate sparse GRM.
    out.status("  Estimating sparse GRM...");
    let grm = estimate::estimate_grm(
        &cohort,
        manifest,
        &pca_scores,
        &unrelated_mask,
        &candidate_pairs,
        n_samples,
        args.n_pcs,
        args.block_size,
        args.degree,
        out,
    )?;
    let n_off_diag = grm.triplets.iter().filter(|&&(i, j, _)| i < j).count();
    out.status(&format!("    {} off-diagonal kinship pairs", n_off_diag));

    // 5. Cache.
    let artifact = GrmArtifact {
        grm,
        pca: pca_scores,
        unrelated: unrel,
        sample_ids: sample_ids.clone(),
    };
    cache::save(&cache_dir, &artifact, &fp, args.degree, args.n_pcs)?;

    let kin_path = cache::grm_tsv_path(&cache_dir);
    let pca_path = cache::pca_tsv_path(&cache_dir);
    out.success(&format!(
        "GRM cached at {}\n  kinship: {}\n  pca:     {}",
        cache_dir.display(),
        kin_path.display(),
        pca_path.display(),
    ));
    out.result_json(&json!({
        "status": "ok",
        "cache_dir": cache_dir.to_string_lossy(),
        "kinship": kin_path.to_string_lossy(),
        "pca": pca_path.to_string_lossy(),
        "n_kinship_pairs": n_off_diag,
        "n_unrelated": artifact.unrelated.sample_indices.len(),
    }));
    Ok(())
}
