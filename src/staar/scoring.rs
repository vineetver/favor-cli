//! Per-chromosome score-test execution.
//!
//! Owns the per-gene and per-window scoring loops that used to live inside
//! `pipeline.rs::run_score_tests`. Dispatch on `ScoringMode` is explicit and
//! happens once per gene, not three times nested.
//!
//! Layout per chromosome:
//! 1. Resolve mask predicates from the requested mask categories
//! 2. Open `SparseG`, `VariantIndex`, and the cached `ChromScoreCache`
//! 3. Score every gene in parallel via `score_one_gene`
//! 4. Score sliding/SCANG windows via `score_one_window`

use std::path::Path;

use faer::Mat;
use rayon::prelude::*;

use crate::error::FavorError;
use crate::output::Output;
use crate::staar::carrier::sparse_score;
use crate::staar::carrier::{AnalysisVectors, CarrierList, VariantIndex, VariantIndexEntry};
use crate::staar::masks::{self, MaskGroup};
use crate::staar::pipeline::{ScoringContext, StaarConfig};
use crate::staar::score;
use crate::staar::score_cache::{self, ChromScoreCache, GeneKBlock};
use crate::staar::sparse_g::SparseG;
use crate::staar::store::StoreManifest;
use crate::staar::{self, GeneResult, MaskCategory, MaskType, ScoringMode};
use crate::types::{AnnotatedVariant, Chromosome};

/// Function-pointer mask predicate over an `AnnotatedVariant`.
type MaskPredicate = fn(&AnnotatedVariant) -> bool;

/// Per-mask result vectors. Order matches the predicate vector built by
/// `MaskPlan::build`.
pub type ResultSet = Vec<(MaskType, Vec<GeneResult>)>;

/// Compiled mask plan: gene predicates + window/scang flags + result-vector
/// indices for window outputs.
struct MaskPlan {
    /// Per-gene mask predicates, in result-vector order.
    gene_predicates: Vec<(MaskType, MaskPredicate)>,
    /// Result vectors, one per mask. Sized to `gene_predicates.len()` plus
    /// optional window/scang slots appended.
    results: ResultSet,
    /// Result-vector index for sliding-window output, if requested.
    window_slot: Option<usize>,
    /// Result-vector index for SCANG output, if requested.
    scang_slot: Option<usize>,
}

impl MaskPlan {
    fn build(categories: &[MaskCategory], out: &dyn Output) -> Self {
        let mut gene_predicates: Vec<(MaskType, MaskPredicate)> = Vec::new();
        let mut want_windows = false;
        let mut want_scang = false;

        for cat in categories {
            match cat {
                MaskCategory::Coding => {
                    for &(ref mt, pred) in masks::CODING_MASKS {
                        gene_predicates.push((mt.clone(), pred));
                    }
                }
                MaskCategory::Noncoding => {
                    for &(ref mt, pred) in masks::NONCODING_MASKS {
                        gene_predicates.push((mt.clone(), pred));
                    }
                }
                MaskCategory::SlidingWindow => want_windows = true,
                MaskCategory::Scang => want_scang = true,
                MaskCategory::Custom => out.warn("Custom BED: not yet implemented"),
            }
        }

        let mut results: ResultSet = gene_predicates
            .iter()
            .map(|(mt, _)| (mt.clone(), Vec::new()))
            .collect();

        let window_slot = want_windows.then(|| {
            results.push((MaskType::SlidingWindow, Vec::new()));
            results.len() - 1
        });
        let scang_slot = want_scang.then(|| {
            results.push((MaskType::Scang, Vec::new()));
            results.len() - 1
        });

        Self {
            gene_predicates,
            results,
            window_slot,
            scang_slot,
        }
    }

    fn has_gene_masks(&self) -> bool {
        !self.gene_predicates.is_empty()
    }

    fn has_window_work(&self) -> bool {
        self.window_slot.is_some() || self.scang_slot.is_some()
    }
}

/// Per-chromosome handles loaded once and shared across all genes.
struct ChromCtx<'a> {
    name: &'a str,
    chrom: Chromosome,
    sparse_g: SparseG,
    variant_index: VariantIndex,
    cache: ChromScoreCache,
}

impl<'a> ChromCtx<'a> {
    fn open(
        store_dir: &Path,
        name: &'a str,
        cache_dir: &Path,
    ) -> Result<Self, FavorError> {
        let chrom_dir = store_dir.join(format!("chromosome={name}"));
        let sparse_g = SparseG::open(&chrom_dir)?;
        let variant_index = VariantIndex::load(&chrom_dir)?;
        let cache = score_cache::load_chromosome(cache_dir, name, &variant_index)?;
        let chrom: Chromosome = name.parse().unwrap_or(Chromosome::Autosome(1));
        Ok(Self {
            name,
            chrom,
            sparse_g,
            variant_index,
            cache,
        })
    }
}

/// Public entry point: score every chromosome in `manifest` and return
/// per-mask gene/window result vectors plus individual p-values.
pub fn run_score_tests(
    store_dir: &Path,
    manifest: &StoreManifest,
    config: &StaarConfig,
    analysis: &AnalysisVectors,
    ctx: &ScoringContext,
    out: &dyn Output,
) -> Result<(ResultSet, Vec<(usize, f64)>), FavorError> {
    out.status("Running score tests (carrier-indexed sparse)...");

    let mut plan = MaskPlan::build(&config.mask_categories, out);

    for ci in &manifest.chromosomes {
        let chrom_ctx = ChromCtx::open(store_dir, &ci.name, &ctx.cache_dir)?;
        out.status(&format!(
            "  chr{}: {} variants, {} genes",
            chrom_ctx.name,
            chrom_ctx.variant_index.len(),
            chrom_ctx.variant_index.n_genes(),
        ));

        if plan.has_gene_masks() {
            score_chrom_genes(&chrom_ctx, &plan.gene_predicates, config, analysis, ctx, &mut plan.results, out);
        }

        if plan.has_window_work() {
            score_chrom_windows(&chrom_ctx, config, analysis, ctx, &mut plan, out);
        }
    }

    Ok((plan.results, Vec::new()))
}

// ---------------------------------------------------------------------------
// Per-gene scoring
// ---------------------------------------------------------------------------

/// All quantities a single gene needs for mask scoring, in a single
/// maf-pass-filtered local frame. `qualifying` indices into the per-mask loop
/// always index into vectors of length `gene_vcfs.len()`.
struct GeneInputs {
    gene_name: String,
    gene_vcfs: Vec<u32>,
    /// One full carrier per local index. Used by SPA + AI-STAAR scorers.
    /// `Standard` mode also needs this only when K is missing from cache.
    carriers: Option<Vec<CarrierList>>,
    /// Local-indexed full U/K. Always materialised; small per gene because
    /// large genes (m > MAX_K_VARIANTS) take the recompute path.
    u_full: Mat<f64>,
    k_full: Mat<f64>,
    start: u32,
    end: u32,
}

/// Score every gene on this chromosome in parallel and append results into
/// `results` in mask-vector order.
fn score_chrom_genes(
    chrom_ctx: &ChromCtx<'_>,
    mask_predicates: &[(MaskType, MaskPredicate)],
    config: &StaarConfig,
    analysis: &AnalysisVectors,
    ctx: &ScoringContext,
    results: &mut ResultSet,
    out: &dyn Output,
) {
    let gene_names: Vec<&String> = chrom_ctx.cache.gene_blocks.keys().collect();
    let mode = ctx.scoring_mode;
    let maf_cutoff = config.maf_cutoff;

    let per_gene: Vec<Vec<(usize, GeneResult)>> = gene_names
        .par_iter()
        .filter_map(|gene_name| {
            let block = chrom_ctx.cache.gene_blocks.get(*gene_name)?;
            let inputs = compile_gene_inputs(
                gene_name,
                block,
                chrom_ctx,
                analysis,
                maf_cutoff,
                mode,
            )?;
            score_gene_masks(&inputs, mask_predicates, chrom_ctx, analysis, ctx)
        })
        .collect();

    for gene_results in per_gene {
        for (mask_idx, result) in gene_results {
            results[mask_idx].1.push(result);
        }
    }

    for (idx, (mask_type, _)) in mask_predicates.iter().enumerate() {
        let n = results[idx].1.len();
        if n > 0 {
            out.status(&format!("    {}: {} groups", mask_type.file_stem(), n));
        }
    }
}

/// Build the per-gene local frame: maf filter, gene_vcfs, and full U/K.
///
/// For `block.has_k() == true` (small/medium genes), U/K come from the cache
/// and are sliced into local order.
///
/// For `block.has_k() == false` (large genes), carriers are loaded from
/// `SparseG` and U/K are recomputed via `score_gene_sparse`.
///
/// In SPA / AI-STAAR mode `carriers` is always loaded since the kernels need
/// raw genotypes; the cached U/K is unused but the work to slice is small.
fn compile_gene_inputs(
    gene_name: &str,
    block: &GeneKBlock,
    chrom_ctx: &ChromCtx<'_>,
    analysis: &AnalysisVectors,
    maf_cutoff: f64,
    mode: ScoringMode,
) -> Option<GeneInputs> {
    if block.m() < 2 {
        return None;
    }
    let all_entries = chrom_ctx.variant_index.all_entries();

    let maf_pass: Vec<usize> = (0..block.m())
        .filter(|&local| {
            let global = block.variant_offsets[local] as usize;
            all_entries[global].maf < maf_cutoff
        })
        .collect();
    if maf_pass.len() < 2 {
        return None;
    }

    let gene_vcfs: Vec<u32> = maf_pass
        .iter()
        .map(|&local| block.variant_offsets[local])
        .collect();

    let positions: Vec<u32> = gene_vcfs
        .iter()
        .map(|&v| all_entries[v as usize].position)
        .collect();
    let start = *positions.iter().min().unwrap();
    let end = *positions.iter().max().unwrap();

    let needs_carriers = !block.has_k() || mode != ScoringMode::Standard;
    let carriers = needs_carriers.then(|| chrom_ctx.sparse_g.load_variants(&gene_vcfs));

    if let Some(c) = &carriers {
        if c.len() < 2 {
            return None;
        }
    }

    let (u_full, k_full) = if block.has_k() {
        // Slice cached chrom-wide U and gene-local K_flat into a maf-filtered frame.
        let m_pass = gene_vcfs.len();
        let u_full = Mat::from_fn(m_pass, 1, |i, _| {
            chrom_ctx.cache.u_all[gene_vcfs[i] as usize]
        });
        let k_full = sparse_score::slice_sumstats_flat(&block.k_flat, block.m(), &maf_pass);
        (u_full, k_full)
    } else {
        // Large gene: K not cached, recompute from carriers.
        let carriers_ref = carriers.as_ref().expect("needs_carriers true when !has_k");
        sparse_score::score_gene_sparse(carriers_ref, analysis)
    };

    Some(GeneInputs {
        gene_name: gene_name.to_string(),
        gene_vcfs,
        carriers,
        u_full,
        k_full,
        start,
        end,
    })
}

/// Score every mask predicate against one gene's compiled inputs.
fn score_gene_masks(
    inputs: &GeneInputs,
    mask_predicates: &[(MaskType, MaskPredicate)],
    chrom_ctx: &ChromCtx<'_>,
    analysis: &AnalysisVectors,
    ctx: &ScoringContext,
) -> Option<Vec<(usize, GeneResult)>> {
    let all_entries = chrom_ctx.variant_index.all_entries();
    let n_vcf = analysis.n_vcf_total;

    let mut results = Vec::new();
    for (mask_idx, (_mt, predicate)) in mask_predicates.iter().enumerate() {
        let qualifying: Vec<usize> = (0..inputs.gene_vcfs.len())
            .filter(|&i| {
                let v = inputs.gene_vcfs[i] as usize;
                predicate(
                    &all_entries[v]
                        .to_annotated_variant_with_gene(chrom_ctx.chrom, &inputs.gene_name),
                )
            })
            .collect();
        if qualifying.len() < 2 {
            continue;
        }

        let mafs: Vec<f64> = qualifying
            .iter()
            .map(|&i| all_entries[inputs.gene_vcfs[i] as usize].maf)
            .collect();
        let ann_matrix = build_annotation_matrix(&inputs.gene_vcfs, &qualifying, all_entries);

        let staar = run_one_test(ctx, &qualifying, &mafs, &ann_matrix, inputs, analysis);

        let cmac: u32 = mafs
            .iter()
            .map(|&maf| (2.0 * maf * n_vcf as f64).round() as u32)
            .sum();

        results.push((
            mask_idx,
            GeneResult {
                ensembl_id: inputs.gene_name.clone(),
                gene_symbol: inputs.gene_name.clone(),
                chromosome: chrom_ctx.chrom,
                start: inputs.start,
                end: inputs.end,
                n_variants: qualifying.len() as u32,
                cumulative_mac: cmac,
                staar,
            },
        ));
    }

    if results.is_empty() {
        None
    } else {
        Some(results)
    }
}

/// Build an `[11 × m]` annotation weight matrix indexed by `qualifying[i]`.
fn build_annotation_matrix(
    gene_vcfs: &[u32],
    qualifying: &[usize],
    all_entries: &[VariantIndexEntry],
) -> Vec<Vec<f64>> {
    (0..11)
        .map(|ch| {
            qualifying
                .iter()
                .map(|&i| all_entries[gene_vcfs[i] as usize].weights.0[ch])
                .collect()
        })
        .collect()
}

/// Dispatch one mask scoring run on the chosen backend.
///
/// `ScoringMode::Standard` slices cached U/K. `Spa` and `AiStaar` need raw
/// carriers; both are guaranteed to have them because `compile_gene_inputs`
/// loaded carriers when `mode != Standard`.
fn run_one_test(
    ctx: &ScoringContext,
    qualifying: &[usize],
    mafs: &[f64],
    ann_matrix: &[Vec<f64>],
    inputs: &GeneInputs,
    analysis: &AnalysisVectors,
) -> score::StaarResult {
    match ctx.scoring_mode {
        ScoringMode::Standard => {
            let (u_sub, k_sub) = sparse_score::slice_sumstats(&inputs.u_full, &inputs.k_full, qualifying);
            score::run_staar_from_sumstats(&u_sub, &k_sub, ann_matrix, mafs, analysis.n_pheno)
        }
        ScoringMode::Spa => {
            let carriers = inputs.carriers.as_ref().expect("Spa requires carriers");
            let subset: Vec<CarrierList> = qualifying.iter().map(|&i| carriers[i].clone()).collect();
            sparse_score::run_staar_sparse(&subset, analysis, ann_matrix, mafs, true)
        }
        ScoringMode::AiStaar => {
            let ai = ctx.ancestry.as_ref().expect("AiStaar requires ancestry info");
            let carriers = inputs.carriers.as_ref().expect("AiStaar requires carriers");
            let subset: Vec<CarrierList> = qualifying.iter().map(|&i| carriers[i].clone()).collect();
            // `--ancestry-col --spa --binary-trait` composes both: SPA fires
            // inside the per-population kernel, AI-STAAR Cauchy-combines.
            staar::ancestry::run_ai_staar_gene(&subset, analysis, ai, ann_matrix, ctx.spa_in_ai_staar)
        }
    }
}

// ---------------------------------------------------------------------------
// Per-window scoring (sliding + SCANG)
// ---------------------------------------------------------------------------

/// Score sliding-window and SCANG groups for one chromosome and append into
/// the result vectors.
fn score_chrom_windows(
    chrom_ctx: &ChromCtx<'_>,
    config: &StaarConfig,
    analysis: &AnalysisVectors,
    ctx: &ScoringContext,
    plan: &mut MaskPlan,
    out: &dyn Output,
) {
    let maf_cutoff = config.maf_cutoff;

    let (chrom_variants, global_indices): (Vec<AnnotatedVariant>, Vec<usize>) = chrom_ctx
        .variant_index
        .all_entries()
        .iter()
        .enumerate()
        .filter(|(_, e)| e.maf < maf_cutoff)
        .map(|(gi, e)| (e.to_annotated_variant(chrom_ctx.chrom), gi))
        .unzip();
    let chrom_indices: Vec<usize> = (0..chrom_variants.len()).collect();
    let n_vcf = analysis.n_vcf_total;

    let score_window = |g: &MaskGroup| -> Option<GeneResult> {
        score_one_window(g, &chrom_variants, &global_indices, chrom_ctx, analysis, ctx, n_vcf)
    };

    if let Some(slot) = plan.window_slot {
        let groups = masks::build_sliding_windows(
            &chrom_variants,
            &chrom_indices,
            chrom_ctx.chrom,
            config.window_size,
            config.window_size / 2,
        );
        if !groups.is_empty() {
            let r: Vec<GeneResult> = groups.iter().filter_map(score_window).collect();
            if !r.is_empty() {
                out.status(&format!("    sliding_window: {} windows", r.len()));
                plan.results[slot].1.extend(r);
            }
        }
    }

    if let Some(slot) = plan.scang_slot {
        let scang_all = masks::build_scang_windows(
            &chrom_variants,
            &chrom_indices,
            chrom_ctx.chrom,
            &config.scang_params,
        );
        for (wsize, groups) in &scang_all {
            let r: Vec<GeneResult> = groups.iter().filter_map(score_window).collect();
            if !r.is_empty() {
                out.status(&format!("    scang L={wsize}: {} windows", r.len()));
                plan.results[slot].1.extend(r);
            }
        }
    }
}

fn score_one_window(
    group: &MaskGroup,
    chrom_variants: &[AnnotatedVariant],
    global_indices: &[usize],
    chrom_ctx: &ChromCtx<'_>,
    analysis: &AnalysisVectors,
    ctx: &ScoringContext,
    n_vcf: usize,
) -> Option<GeneResult> {
    let m = group.variant_indices.len();
    if m < 2 {
        return None;
    }

    let win_globals: Vec<usize> = group
        .variant_indices
        .iter()
        .map(|&ci| global_indices[ci])
        .collect();
    let mafs: Vec<f64> = group
        .variant_indices
        .iter()
        .map(|&ci| chrom_variants[ci].maf)
        .collect();
    let ann_matrix: Vec<Vec<f64>> = (0..11)
        .map(|ch| {
            group
                .variant_indices
                .iter()
                .map(|&ci| chrom_variants[ci].annotation.weights.0[ch])
                .collect()
        })
        .collect();

    let staar = match ctx.scoring_mode {
        ScoringMode::AiStaar => {
            let ai = ctx
                .ancestry
                .as_ref()
                .expect("AiStaar requires ancestry info");
            let win_vcfs: Vec<u32> = win_globals.iter().map(|&i| i as u32).collect();
            let carriers = chrom_ctx.sparse_g.load_variants(&win_vcfs);
            staar::ancestry::run_ai_staar_gene(&carriers, analysis, ai, &ann_matrix, false)
        }
        // SPA on windows is not currently wired (windows always go through
        // cached U/K). Treat Standard and Spa the same here.
        ScoringMode::Standard | ScoringMode::Spa => {
            let u_win = score_cache::slice_window_u(&chrom_ctx.cache, &win_globals);
            let k_win = score_cache::assemble_window_k(&chrom_ctx.cache, &win_globals);
            score::run_staar_from_sumstats(&u_win, &k_win, &ann_matrix, &mafs, analysis.n_pheno)
        }
    };

    let cmac: u32 = mafs
        .iter()
        .map(|&maf| (2.0 * maf * n_vcf as f64).round() as u32)
        .sum();

    Some(GeneResult {
        ensembl_id: group.name.clone(),
        gene_symbol: group.name.clone(),
        chromosome: group.chromosome,
        start: group.start,
        end: group.end,
        n_variants: m as u32,
        cumulative_mac: cmac,
        staar,
    })
}

