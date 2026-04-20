//! Per-chromosome score-test execution.
//!
//! Owns the per-gene and per-window scoring loops that used to live inside
//! `pipeline.rs::run_score_tests`. Dispatch on `ScoringMode` is explicit and
//! happens once per gene, not three times nested.
//!
//! Layout per chromosome:
//! 1. Resolve mask predicates from the requested mask categories
//! 2. Open a ChromosomeView and the cached ChromScoreCache via that view
//! 3. Score every gene in parallel via `score_one_gene`
//! 4. Score sliding/SCANG windows via `score_one_window`

use std::path::Path;

use faer::Mat;
use rayon::prelude::*;

use crate::error::CohortError;
use crate::output::Output;
use crate::staar::carrier::sparse_score::{self, carriers_to_dense_compact};
use crate::staar::carrier::AnalysisVectors;
use crate::staar::masks::{self, MaskGroup};
use crate::staar::multi::{self, MultiNull};
use crate::staar::score;
use crate::staar::{self, GeneResult, MaskCategory, MaskType, ScoringMode};
use crate::store::cache::score_cache::{self, ChromScoreCache, GeneKBlock};
use crate::store::cohort::types::{from_u32_slice, VariantVcf};
use crate::store::cohort::variants::{CarrierList, VariantIndexEntry};
use crate::store::cohort::{ChromosomeView, CohortHandle, CohortManifest};
use crate::types::{AnnotatedVariant, Chromosome};

use super::ancestry::AncestryInfo;
use super::masks::ScangParams;

/// Function-pointer mask predicate over an `AnnotatedVariant`.
type MaskPredicate = fn(&AnnotatedVariant) -> bool;

/// (window_size, scored_windows) batch for SCANG MC threshold computation.
type ScoredScangBatch = Vec<(u32, Vec<(GeneResult, Vec<usize>)>)>;

/// Per-mask result vectors. Order matches the predicate vector built by
/// `MaskPlan::build`.
pub type ResultSet = Vec<(MaskType, Vec<GeneResult>)>;

pub struct ScoringRequest<'a> {
    pub mask_categories: &'a [MaskCategory],
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub scang_params: &'a ScangParams,
    pub mode: ScoringMode,
    pub cache_dir: &'a Path,
    pub ancestry: Option<&'a AncestryInfo>,
    pub spa_active: bool,
    /// Emit per-variant score-test rows alongside mask/window results.
    /// Matches R STAARpipeline's `Individual_Analysis` output (single-trait
    /// gaussian no-kinship path only in this build).
    pub individual: bool,
    /// Minor-allele-count threshold for individual analysis. Matches R's
    /// `mac_cutoff` (default 20 in `Individual_Analysis.R`).
    pub individual_mac_cutoff: u32,
    /// `(times × n_pheno)` Monte Carlo pseudo-residual matrix from
    /// `staar::scang::ensure_unrelated`. Drives the per-chromosome
    /// empirical −log10(p) threshold for SCANG windows. `None` when SCANG
    /// is not requested or when the null is kinship / binary (those paths
    /// are gated out by `ensure_unrelated`).
    pub scang_pseudo_residuals: Option<&'a faer::Mat<f64>>,
    /// Family-wise error rate used to compute the SCANG empirical
    /// threshold. SCANG default from R/SCANG.r:30.
    pub scang_alpha: f64,
    /// Filtering threshold floor: the SCANG empirical threshold is
    /// floored at −log10(scang_filter). R/SCANG.r:32.
    pub scang_filter: f64,
}

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
    fn build(categories: &[MaskCategory]) -> Self {
        let mut gene_predicates: Vec<(MaskType, MaskPredicate)> = Vec::new();
        let mut want_windows = false;
        let mut want_scang = false;

        for cat in categories {
            match cat {
                MaskCategory::Coding => gene_predicates.extend_from_slice(masks::CODING_MASKS),
                MaskCategory::Noncoding => {
                    gene_predicates.extend_from_slice(masks::NONCODING_MASKS)
                }
                MaskCategory::SlidingWindow => want_windows = true,
                MaskCategory::Scang => want_scang = true,
            }
        }

        let mut results: ResultSet = gene_predicates
            .iter()
            .map(|(mt, _)| (*mt, Vec::new()))
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

/// Per-chromosome handles loaded once and shared across all genes. Holds
/// a `ChromosomeView` (which lazily mmap's `sparse_g.bin` and parses
/// `variants.parquet`) plus the per-phenotype score cache.
struct ChromCtx<'v> {
    name: String,
    chrom: Chromosome,
    view: ChromosomeView<'v>,
    cache: ChromScoreCache,
}

impl<'v> ChromCtx<'v> {
    fn open(
        cohort: &'v CohortHandle<'v>,
        chrom: Chromosome,
        cache_dir: &Path,
    ) -> Result<Self, CohortError> {
        let view = cohort.chromosome(&chrom)?;
        let cache = score_cache::load_chromosome(cache_dir, &view)?;
        Ok(Self {
            name: chrom.label(),
            chrom,
            view,
            cache,
        })
    }
}

/// Public entry point: score every chromosome in `manifest` and return
/// per-mask gene/window result vectors plus per-variant individual rows.
pub fn run_score_tests(
    cohort: &CohortHandle<'_>,
    manifest: &CohortManifest,
    request: &ScoringRequest<'_>,
    analysis: &AnalysisVectors,
    out: &dyn Output,
) -> Result<(ResultSet, Vec<crate::staar::output::IndividualRow>), CohortError> {
    out.status("Running score tests (carrier-indexed sparse)...");

    let mut plan = MaskPlan::build(request.mask_categories);
    let mut individual_rows: Vec<crate::staar::output::IndividualRow> = Vec::new();

    // `score_one_window` always takes the cached-U/K path; SPA on windows
    // isn't wired yet. Warn once so --spa users know their window p-values
    // aren't SPA-corrected.
    if request.mode == ScoringMode::Spa && plan.has_window_work() {
        out.warn("--spa: window scoring uses non-SPA p-values; SPA only applies to gene masks");
    }

    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
        let chrom_ctx = ChromCtx::open(cohort, chrom, request.cache_dir)?;
        let n_variants = chrom_ctx.view.index()?.len();
        let n_genes = chrom_ctx.view.index()?.n_genes();
        out.status(&format!(
            "  chr{}: {} variants, {} genes",
            chrom_ctx.name, n_variants, n_genes,
        ));

        if plan.has_gene_masks() {
            score_chrom_genes(
                &chrom_ctx,
                &plan.gene_predicates,
                request,
                analysis,
                &mut plan.results,
                out,
            )?;
        }

        if plan.has_window_work() {
            score_chrom_windows(&chrom_ctx, request, analysis, &mut plan, out);
        }

        if request.individual {
            score_chrom_individual(&chrom_ctx, request, analysis, &mut individual_rows)?;
        }
    }

    Ok((plan.results, individual_rows))
}

/// Per-variant score test for every variant whose MAC clears `mac_cutoff`.
/// Mirrors R STAARpipeline `Individual_Analysis` on the gaussian-unrelated
/// path. Writer re-sorts by position to match R's final `order(results[,2])`.
fn score_chrom_individual(
    chrom_ctx: &ChromCtx<'_>,
    request: &ScoringRequest<'_>,
    analysis: &AnalysisVectors,
    rows: &mut Vec<crate::staar::output::IndividualRow>,
) -> Result<(), CohortError> {
    use crate::store::cohort::types::VariantVcf;

    let index = chrom_ctx.view.index()?;
    let entries = index.all_entries();
    let n_pheno = analysis.n_pheno as f64;
    let mac_cutoff = request.individual_mac_cutoff;
    let chrom_label = chrom_ctx.chrom.label();

    // R `Individual_Analysis.R:324`: MAF > (mac_cutoff − 0.5) / (2 n).
    let maf_floor = (mac_cutoff as f64 - 0.5) / (2.0 * n_pheno);

    let mut passing: Vec<(u32, &crate::store::cohort::variants::VariantIndexEntry)> =
        Vec::with_capacity(entries.len() / 100);
    for (vcf, entry) in entries.iter().enumerate() {
        if entry.maf > maf_floor {
            passing.push((vcf as u32, entry));
        }
    }
    if passing.is_empty() {
        return Ok(());
    }

    // `all_entries` is vcf-ordered, so `passing` is already sorted as
    // `carriers_batch` requires for a single sequential mmap walk.
    let vcfs: Vec<VariantVcf> = passing.iter().map(|(v, _)| VariantVcf(*v)).collect();
    let batch = chrom_ctx.view.carriers_batch(&vcfs)?;

    // R AI_Individual_Analysis.R:52-55 returns NULL when use_SPA=TRUE, and
    // we only validate the gaussian-unrelated path. Gate the AI call on both.
    let ai_ctx = request
        .ancestry
        .filter(|_| !request.spa_active && analysis.kinship.is_none());

    let mut chrom_rows: Vec<crate::staar::output::IndividualRow> = passing
        .par_iter()
        .zip(batch.entries.par_iter())
        .map(|(&(_vcf, entry), carriers)| {
            let result = crate::staar::carrier::sparse_score::individual_score_test(
                carriers, analysis,
            );
            let pvalue_log10 = if result.pvalue > 0.0 {
                -result.pvalue.log10()
            } else {
                f64::INFINITY
            };
            let (ai_pvalue, pop_mafs) = match ai_ctx {
                Some(ancestry) => {
                    let (p, mafs) = crate::staar::ancestry::ai_individual_score_test(
                        carriers, analysis, ancestry,
                    );
                    (Some(p), mafs)
                }
                None => (None, Vec::new()),
            };
            crate::staar::output::IndividualRow {
                chromosome: chrom_label.clone(),
                position: entry.position,
                ref_allele: entry.ref_allele.clone(),
                alt_allele: entry.alt_allele.clone(),
                maf: entry.maf,
                alt_af: result.alt_af,
                n: analysis.n_pheno as u32,
                mac: result.mac,
                pvalue: result.pvalue,
                pvalue_log10,
                score: result.score,
                score_se: result.score_se,
                est: result.est,
                est_se: result.est_se,
                ai_pvalue,
                pop_mafs,
            }
        })
        .collect();

    rows.append(&mut chrom_rows);
    Ok(())
}

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

/// Score every gene on this chromosome in parallel and append results
/// into `results` in mask-vector order. Per-gene returns
/// `Result<Option<_>>` so legitimate skips (no qualifying variants) and
/// real failures don't both look like `None` from a `filter_map`.
fn score_chrom_genes(
    chrom_ctx: &ChromCtx<'_>,
    mask_predicates: &[(MaskType, MaskPredicate)],
    request: &ScoringRequest<'_>,
    analysis: &AnalysisVectors,
    results: &mut ResultSet,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let gene_names: Vec<&String> = chrom_ctx.cache.gene_blocks.keys().collect();
    let mode = request.mode;
    let maf_cutoff = request.maf_cutoff;

    let per_gene: Vec<Vec<(usize, GeneResult)>> = gene_names
        .par_iter()
        .map(
            |gene_name| -> Result<Option<Vec<(usize, GeneResult)>>, CohortError> {
                let Some(block) = chrom_ctx.cache.gene_blocks.get(*gene_name) else {
                    return Ok(None);
                };
                let Some(inputs) =
                    compile_gene_inputs(gene_name, block, chrom_ctx, analysis, maf_cutoff, mode)?
                else {
                    return Ok(None);
                };
                score_gene_masks(&inputs, mask_predicates, chrom_ctx, analysis, request)
            },
        )
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .flatten()
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
    Ok(())
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
) -> Result<Option<GeneInputs>, CohortError> {
    if block.m() < 2 {
        return Ok(None);
    }
    let variant_index = chrom_ctx.view.index()?;
    let all_entries = variant_index.all_entries();

    let maf_pass: Vec<usize> = (0..block.m())
        .filter(|&local| {
            let global = block.variant_offsets[local] as usize;
            all_entries[global].maf < maf_cutoff
        })
        .collect();
    if maf_pass.len() < 2 {
        return Ok(None);
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
    let carriers = if needs_carriers {
        let vcfs = crate::store::cohort::types::from_u32_slice(&gene_vcfs);
        Some(chrom_ctx.view.carriers_batch(vcfs)?.entries)
    } else {
        None
    };

    if let Some(c) = &carriers {
        if c.len() < 2 {
            return Ok(None);
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

    Ok(Some(GeneInputs {
        gene_name: gene_name.to_string(),
        gene_vcfs,
        carriers,
        u_full,
        k_full,
        start,
        end,
    }))
}

/// Score every mask predicate against one gene. `Ok(None)` is the "no
/// mask had ≥2 qualifying variants" skip; real errors propagate.
fn score_gene_masks(
    inputs: &GeneInputs,
    mask_predicates: &[(MaskType, MaskPredicate)],
    chrom_ctx: &ChromCtx<'_>,
    analysis: &AnalysisVectors,
    request: &ScoringRequest<'_>,
) -> Result<Option<Vec<(usize, GeneResult)>>, CohortError> {
    let all_entries = chrom_ctx.view.index()?.all_entries();
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

        let staar = run_one_test(request, &qualifying, &mafs, &ann_matrix, inputs, analysis);

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
                emthr: f64::NAN,
            },
        ));
    }

    if results.is_empty() {
        Ok(None)
    } else {
        Ok(Some(results))
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
    request: &ScoringRequest<'_>,
    qualifying: &[usize],
    mafs: &[f64],
    ann_matrix: &[Vec<f64>],
    inputs: &GeneInputs,
    analysis: &AnalysisVectors,
) -> score::StaarResult {
    match request.mode {
        ScoringMode::Standard => {
            let (u_sub, k_sub) =
                sparse_score::slice_sumstats(&inputs.u_full, &inputs.k_full, qualifying);
            score::run_staar_from_sumstats(&u_sub, &k_sub, ann_matrix, mafs, analysis.n_pheno)
        }
        ScoringMode::Spa => {
            let carriers = inputs.carriers.as_ref().expect("Spa requires carriers");
            let subset: Vec<CarrierList> =
                qualifying.iter().map(|&i| carriers[i].clone()).collect();
            sparse_score::run_staar_sparse(&subset, analysis, ann_matrix, mafs, true)
        }
        ScoringMode::AiStaar => {
            let ai = request.ancestry.expect("AiStaar requires ancestry info");
            let carriers = inputs.carriers.as_ref().expect("AiStaar requires carriers");
            let subset: Vec<CarrierList> =
                qualifying.iter().map(|&i| carriers[i].clone()).collect();
            // `--ancestry-col --spa --binary-trait` composes both: SPA fires
            // inside the per-population kernel, AI-STAAR Cauchy-combines.
            staar::ancestry::run_ai_staar_gene(
                &subset,
                analysis,
                ai,
                ann_matrix,
                request.spa_active,
            )
        }
    }
}

/// Score sliding-window and SCANG groups for one chromosome and append into
/// the result vectors.
fn score_chrom_windows(
    chrom_ctx: &ChromCtx<'_>,
    request: &ScoringRequest<'_>,
    analysis: &AnalysisVectors,
    plan: &mut MaskPlan,
    out: &dyn Output,
) {
    let maf_cutoff = request.maf_cutoff;

    let variant_index = match chrom_ctx.view.index() {
        Ok(i) => i,
        Err(e) => {
            out.warn(&format!(
                "    skipping windows on chr{}: {e}",
                chrom_ctx.name
            ));
            return;
        }
    };
    let (chrom_variants, global_indices): (Vec<AnnotatedVariant>, Vec<usize>) = variant_index
        .all_entries()
        .iter()
        .enumerate()
        .filter(|(_, e)| e.maf < maf_cutoff)
        .map(|(gi, e)| (e.to_annotated_variant(chrom_ctx.chrom), gi))
        .unzip();
    let chrom_indices: Vec<usize> = (0..chrom_variants.len()).collect();
    let n_vcf = analysis.n_vcf_total;

    let score_window = |g: &MaskGroup| -> Option<GeneResult> {
        score_one_window(
            g,
            &chrom_variants,
            &global_indices,
            chrom_ctx,
            analysis,
            request,
            n_vcf,
        )
    };

    if let Some(slot) = plan.window_slot {
        let groups = masks::build_sliding_windows(
            &chrom_variants,
            &chrom_indices,
            chrom_ctx.chrom,
            request.window_size,
            request.window_size / 2,
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
            request.scang_params,
        );
        let mut all_scored: ScoredScangBatch = scang_all
            .iter()
            .map(|(wsize, groups)| {
                let scored: Vec<(GeneResult, Vec<usize>)> = groups
                    .iter()
                    .filter_map(|g| {
                        let gr = score_window(g)?;
                        let win_globals: Vec<usize> = g
                            .variant_indices
                            .iter()
                            .map(|&ci| global_indices[ci])
                            .collect();
                        Some((gr, win_globals))
                    })
                    .collect();
                (*wsize, scored)
            })
            .collect();

        let emthr = request.scang_pseudo_residuals.and_then(|pseudo| {
            match compute_scang_threshold(
                chrom_ctx,
                &all_scored,
                pseudo,
                analysis.n_pheno,
                request.scang_alpha,
                request.scang_filter,
                out,
            ) {
                Ok(th) => Some(th),
                Err(e) => {
                    out.warn(&format!(
                        "    scang MC threshold on chr{}: {e}; falling back to filter floor",
                        chrom_ctx.name,
                    ));
                    None
                }
            }
        });
        let floor = -request.scang_filter.ln() / std::f64::consts::LN_10;
        let threshold = emthr.unwrap_or(floor);

        for (wsize, scored) in all_scored.iter_mut() {
            let kept: Vec<GeneResult> = scored
                .drain(..)
                .filter_map(|(mut gr, _)| {
                    gr.emthr = threshold;
                    let neg_log10 = if gr.staar.staar_o > 0.0 {
                        -gr.staar.staar_o.log10()
                    } else {
                        f64::INFINITY
                    };
                    if neg_log10 >= threshold {
                        Some(gr)
                    } else {
                        None
                    }
                })
                .collect();
            if !kept.is_empty() {
                out.status(&format!(
                    "    scang L={wsize}: {} windows (th0={:.3})",
                    kept.len(),
                    threshold,
                ));
                plan.results[slot].1.extend(kept);
            }
        }
    }
}

/// Compute the per-chromosome SCANG empirical −log10(p) threshold using
/// `pseudo_residuals` from `staar::scang::ensure_unrelated`.
///
/// Algorithm mirrors SCANG R/SCANG.r:297-328:
/// 1. Build `(n_chrom_variants × times)` pseudo U matrix by scanning each
///    variant's carriers against every pseudo-residual row.
/// 2. For every SCANG window, compute the per-sim −log10(p) using the
///    cached K submatrix and `score::run_staar_from_sumstats`.
/// 3. Track the max stat per simulation across all windows.
/// 4. Return `max(quantile(1-alpha, max_stats), -log10(filter))`.
#[allow(clippy::type_complexity)]
fn compute_scang_threshold(
    chrom_ctx: &ChromCtx<'_>,
    scored: &[(u32, Vec<(GeneResult, Vec<usize>)>)],
    pseudo: &faer::Mat<f64>,
    n_pheno: usize,
    alpha: f64,
    filter: f64,
    out: &dyn Output,
) -> Result<f64, CohortError> {
    use faer::Mat;

    let times = pseudo.nrows();
    let n_pheno_pseudo = pseudo.ncols();
    if n_pheno_pseudo != n_pheno {
        return Err(CohortError::Input(format!(
            "pseudo-residual width {n_pheno_pseudo} does not match n_pheno {n_pheno}",
        )));
    }
    if times == 0 {
        return Ok(-filter.ln() / std::f64::consts::LN_10);
    }

    let variant_index = chrom_ctx.view.index()?;
    let n_variants = variant_index.len();
    let all_vcfs: Vec<crate::store::cohort::types::VariantVcf> = (0..n_variants as u32)
        .map(crate::store::cohort::types::VariantVcf)
        .collect();
    let carriers = chrom_ctx.view.carriers_batch(&all_vcfs)?.entries;

    // u_sim[(variant_id, sim_j)] = G_j' z_j for this chromosome's variant
    // against row j of the pseudo-residuals. One scan per variant.
    let mut u_sim = Mat::<f64>::zeros(n_variants, times);
    for (gi, carrier) in carriers.iter().enumerate() {
        for entry in &carrier.entries {
            if entry.dosage == 255 {
                continue;
            }
            let pi = entry.sample_idx as usize;
            let d = entry.dosage as f64;
            for t in 0..times {
                u_sim[(gi, t)] += d * pseudo[(t, pi)];
            }
        }
    }

    let mut max_stat = vec![0.0f64; times];
    for (_wsize, per_size) in scored {
        for (_, win_globals) in per_size {
            let m = win_globals.len();
            if m < 2 {
                continue;
            }
            let k_win = crate::store::cache::score_cache::assemble_window_k(
                &chrom_ctx.cache,
                win_globals,
            );
            // Per-sim weighted burden stat: T_j = (1' u_win_j)^2 / (1' K_win 1).
            // Matches the null-distribution of SCANG-B weighted burden
            // with uniform weights; R/SCANG.r fires a full SKAT + Burden
            // + ACAT-O omnibus per sim. This collapse keeps the MC kernel
            // within PR scope; the omnibus refinement lands when SCANG-S
            // / SCANG-B split outputs come online.
            let mut one_k_one = 0.0;
            for i in 0..m {
                for j in 0..m {
                    one_k_one += k_win[(i, j)];
                }
            }
            if !(one_k_one.is_finite() && one_k_one > 0.0) {
                continue;
            }
            for t in 0..times {
                let mut sum_u = 0.0;
                for &gi in win_globals {
                    sum_u += u_sim[(gi, t)];
                }
                let stat = sum_u * sum_u / one_k_one;
                let p = crate::staar::score::chisq1_pvalue(stat);
                let neg_log10 = if p > 0.0 {
                    -p.log10()
                } else {
                    f64::INFINITY.min(1e308)
                };
                if neg_log10 > max_stat[t] {
                    max_stat[t] = neg_log10;
                }
            }
        }
    }

    let th = crate::staar::scang::chrom_threshold(&max_stat, alpha, filter);
    out.status(&format!(
        "    scang MC threshold on chr{} = {th:.3} ({} sims × windows)",
        chrom_ctx.name, times,
    ));
    Ok(th)
}

fn score_one_window(
    group: &MaskGroup,
    chrom_variants: &[AnnotatedVariant],
    global_indices: &[usize],
    chrom_ctx: &ChromCtx<'_>,
    analysis: &AnalysisVectors,
    request: &ScoringRequest<'_>,
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

    let staar = match request.mode {
        ScoringMode::AiStaar => {
            let ai = request.ancestry.expect("AiStaar requires ancestry info");
            let win_vcfs: Vec<crate::store::cohort::types::VariantVcf> = win_globals
                .iter()
                .map(|&i| crate::store::cohort::types::VariantVcf(i as u32))
                .collect();
            let carriers = match chrom_ctx.view.carriers_batch(&win_vcfs) {
                Ok(b) => b.entries,
                Err(_) => return None,
            };
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
        emthr: f64::NAN,
    })
}

// =====================================================================
// Joint multi-trait scoring path
// =====================================================================
//
// This section mirrors the single-trait loop above but never builds an
// `AnalysisVectors` or touches the on-disk score cache. The multi kernel
// in `staar::multi` wants a dense `(n_pheno × m)` genotype matrix per
// gene plus the joint null model; there is no per-variant U/K
// representation to cache on the multi path. Gene iteration, mask plan,
// annotation matrix construction, and window building are all reused
// from the single-trait code above.

/// Arguments for the joint multi-trait scoring run. Narrower than
/// `ScoringRequest` because the multi path has no `ScoringMode` axis —
/// SPA, AI-STAAR, and kinship-aware scoring are all rejected at config
/// build time (see `commands::staar::build_config`).
pub struct MultiScoringRequest<'a> {
    pub mask_categories: &'a [MaskCategory],
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub scang_params: &'a ScangParams,
}

/// Public entry point: score every chromosome in `manifest` against the
/// joint multi-trait null model and return per-mask gene/window result
/// vectors. Mirrors `run_score_tests` but never touches a score cache.
pub fn run_multi_score_tests(
    cohort: &CohortHandle<'_>,
    manifest: &CohortManifest,
    request: &MultiScoringRequest<'_>,
    null: &MultiNull,
    pheno_mask: &[bool],
    out: &dyn Output,
) -> Result<ResultSet, CohortError> {
    out.status("Running joint multi-trait score tests...");

    // Compact VCF → phenotype index map. Same compaction rule as
    // `AnalysisVectors::from_null_model`: samples with phenotype data get
    // consecutive indices starting at 0; the rest stay None so carriers
    // touching them are skipped.
    let vcf_to_pheno: Vec<Option<u32>> = {
        let mut map = Vec::with_capacity(pheno_mask.len());
        let mut next: u32 = 0;
        for &has in pheno_mask {
            if has {
                map.push(Some(next));
                next += 1;
            } else {
                map.push(None);
            }
        }
        map
    };
    debug_assert_eq!(vcf_to_pheno.iter().filter(|p| p.is_some()).count(), null.n_samples);
    let n_pheno = null.n_samples;

    let mut plan = MaskPlan::build(request.mask_categories);

    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
        let view = cohort.chromosome(&chrom)?;
        let n_variants = view.index()?.len();
        let n_genes = view.index()?.n_genes();
        out.status(&format!(
            "  chr{}: {} variants, {} genes",
            chrom.label(),
            n_variants,
            n_genes,
        ));

        if plan.has_gene_masks() {
            score_chrom_genes_multi(
                &view,
                &vcf_to_pheno,
                n_pheno,
                null,
                &plan.gene_predicates,
                request.maf_cutoff,
                chrom,
                &mut plan.results,
                out,
            )?;
        }

        if plan.has_window_work() {
            score_chrom_windows_multi(
                &view,
                &vcf_to_pheno,
                n_pheno,
                null,
                request,
                chrom,
                &mut plan,
                out,
            );
        }
    }

    Ok(plan.results)
}

/// Per-chromosome multi-trait gene loop. Walks each gene in the variant
/// index, maf-filters the gene's variants, loads carriers, materializes
/// one `(n_pheno × m)` genotype matrix per gene, then evaluates every
/// mask predicate against a column slice of that matrix.
#[allow(clippy::too_many_arguments)]
fn score_chrom_genes_multi(
    view: &ChromosomeView<'_>,
    vcf_to_pheno: &[Option<u32>],
    n_pheno: usize,
    null: &MultiNull,
    mask_predicates: &[(MaskType, MaskPredicate)],
    maf_cutoff: f64,
    chrom: Chromosome,
    results: &mut ResultSet,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let variant_index = view.index()?;
    let all_entries = variant_index.all_entries();
    let gene_names: Vec<String> = variant_index.gene_names().map(|s| s.to_string()).collect();
    let n_vcf = vcf_to_pheno.len();

    let per_gene: Vec<Vec<(usize, GeneResult)>> = gene_names
        .par_iter()
        .map(|gene_name| -> Result<Vec<(usize, GeneResult)>, CohortError> {
            let raw_vcfs: &[u32] = variant_index.gene_variant_vcfs(gene_name);
            if raw_vcfs.len() < 2 {
                return Ok(Vec::new());
            }

            // MAF filter → indices of variants that pass, plus their
            // global VCF positions and positions on the chromosome.
            let gene_vcfs: Vec<u32> = raw_vcfs
                .iter()
                .copied()
                .filter(|&v| all_entries[v as usize].maf < maf_cutoff)
                .collect();
            if gene_vcfs.len() < 2 {
                return Ok(Vec::new());
            }

            let positions: Vec<u32> = gene_vcfs
                .iter()
                .map(|&v| all_entries[v as usize].position)
                .collect();
            let start = *positions.iter().min().unwrap();
            let end = *positions.iter().max().unwrap();

            // `gene_variant_vcfs` is sorted by the membership loader, so
            // the maf filter preserves order and `carriers_batch` sees a
            // monotonic slice as required.
            debug_assert!(gene_vcfs.windows(2).all(|w| w[0] <= w[1]));
            let vcfs: &[VariantVcf] = from_u32_slice(&gene_vcfs);
            let carriers = view.carriers_batch(vcfs)?.entries;
            if carriers.len() < 2 {
                return Ok(Vec::new());
            }
            let g0 = carriers_to_dense_compact(&carriers, vcf_to_pheno, n_pheno);

            let mut gene_results = Vec::new();
            for (mask_idx, (_mt, predicate)) in mask_predicates.iter().enumerate() {
                let qualifying: Vec<usize> = (0..gene_vcfs.len())
                    .filter(|&i| {
                        let v = gene_vcfs[i] as usize;
                        predicate(
                            &all_entries[v]
                                .to_annotated_variant_with_gene(chrom, gene_name),
                        )
                    })
                    .collect();
                if qualifying.len() < 2 {
                    continue;
                }

                let g_mask = slice_g0_columns(&g0, &qualifying);
                let mafs: Vec<f64> = qualifying
                    .iter()
                    .map(|&i| all_entries[gene_vcfs[i] as usize].maf)
                    .collect();
                let ann_matrix = build_annotation_matrix(&gene_vcfs, &qualifying, all_entries);

                let staar = multi::run_multi_staar(&g_mask, null, &ann_matrix, &mafs);

                let cmac: u32 = mafs
                    .iter()
                    .map(|&maf| (2.0 * maf * n_vcf as f64).round() as u32)
                    .sum();

                gene_results.push((
                    mask_idx,
                    GeneResult {
                        ensembl_id: gene_name.clone(),
                        gene_symbol: gene_name.clone(),
                        chromosome: chrom,
                        start,
                        end,
                        n_variants: qualifying.len() as u32,
                        cumulative_mac: cmac,
                        staar,
                        emthr: f64::NAN,
                    },
                ));
            }
            Ok(gene_results)
        })
        .collect::<Result<Vec<_>, _>>()?;

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

    Ok(())
}

/// Per-chromosome multi-trait window loop. Mirrors `score_chrom_windows`
/// but builds a dense G₀ from carriers and calls the joint kernel.
#[allow(clippy::too_many_arguments)]
fn score_chrom_windows_multi(
    view: &ChromosomeView<'_>,
    vcf_to_pheno: &[Option<u32>],
    n_pheno: usize,
    null: &MultiNull,
    request: &MultiScoringRequest<'_>,
    chrom: Chromosome,
    plan: &mut MaskPlan,
    out: &dyn Output,
) {
    let variant_index = match view.index() {
        Ok(i) => i,
        Err(e) => {
            out.warn(&format!("    skipping windows on chr{}: {e}", chrom.label()));
            return;
        }
    };
    let (chrom_variants, global_indices): (Vec<AnnotatedVariant>, Vec<usize>) = variant_index
        .all_entries()
        .iter()
        .enumerate()
        .filter(|(_, e)| e.maf < request.maf_cutoff)
        .map(|(gi, e)| (e.to_annotated_variant(chrom), gi))
        .unzip();
    let chrom_indices: Vec<usize> = (0..chrom_variants.len()).collect();
    let n_vcf = vcf_to_pheno.len();

    let score_window = |group: &MaskGroup| -> Option<GeneResult> {
        score_one_window_multi(
            group,
            &chrom_variants,
            &global_indices,
            view,
            vcf_to_pheno,
            n_pheno,
            null,
            n_vcf,
        )
    };

    if let Some(slot) = plan.window_slot {
        let groups = masks::build_sliding_windows(
            &chrom_variants,
            &chrom_indices,
            chrom,
            request.window_size,
            request.window_size / 2,
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
            chrom,
            request.scang_params,
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

/// Score one multi-trait window. Returns `None` for windows with fewer
/// than two carriers after loading, to match the single-trait path's
/// skip semantics.
#[allow(clippy::too_many_arguments)]
fn score_one_window_multi(
    group: &MaskGroup,
    chrom_variants: &[AnnotatedVariant],
    global_indices: &[usize],
    view: &ChromosomeView<'_>,
    vcf_to_pheno: &[Option<u32>],
    n_pheno: usize,
    null: &MultiNull,
    n_vcf: usize,
) -> Option<GeneResult> {
    let m = group.variant_indices.len();
    if m < 2 {
        return None;
    }

    let win_globals: Vec<u32> = group
        .variant_indices
        .iter()
        .map(|&ci| global_indices[ci] as u32)
        .collect();
    // Windows emit variant_indices in position order, so globals are
    // monotonic — `carriers_batch`'s sorted-input contract is honoured.
    debug_assert!(win_globals.windows(2).all(|w| w[0] <= w[1]));
    let vcfs: &[VariantVcf] = from_u32_slice(&win_globals);
    let carriers = view.carriers_batch(vcfs).ok()?.entries;
    if carriers.len() < 2 {
        return None;
    }
    let g = carriers_to_dense_compact(&carriers, vcf_to_pheno, n_pheno);

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

    let staar = multi::run_multi_staar(&g, null, &ann_matrix, &mafs);

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
        emthr: f64::NAN,
    })
}

/// Slice a subset of columns out of a dense `(n_pheno, m)` genotype
/// matrix. The multi path builds one G₀ per gene and evaluates every
/// mask predicate against a column subset of that matrix, so this
/// small helper lives next to the loop that calls it.
fn slice_g0_columns(g: &Mat<f64>, columns: &[usize]) -> Mat<f64> {
    let n = g.nrows();
    let k = columns.len();
    Mat::from_fn(n, k, |i, j| g[(i, columns[j])])
}
