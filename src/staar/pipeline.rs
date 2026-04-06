//! STAAR analysis pipeline: ensure_store → fit_null_model → score_all → write_results.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use faer::Mat;
use rayon::prelude::*;

use crate::column::STAAR_WEIGHTS;
use crate::data::{AnnotatedSet, VariantSet, VariantSetKind};
use crate::error::FavorError;
use crate::ingest::ColumnContract;
use crate::output::Output;
use crate::staar::carrier::AnalysisVectors;
use crate::staar::carrier::sparse_score;
use crate::staar::masks::{self, MaskGroup, ScangParams};
use crate::staar::sparse_g::SparseG;
use crate::staar::model::{augment_covariates, load_known_loci, load_phenotype};
use crate::staar::output::{write_individual_results, write_results};
use crate::staar::score;
use crate::staar::score_cache;
use crate::staar::store::{self, GenoStoreResult, STAAR_ANNOTATION_COLUMNS};
use crate::staar::{self, GeneResult, MaskCategory, MaskType, TraitType};
use crate::resource::Resources;
use crate::types::{AnnotatedVariant, Chromosome};

// ---------------------------------------------------------------------------
// Config
// ---------------------------------------------------------------------------

pub struct StaarConfig {
    pub genotypes: PathBuf,
    pub phenotype: PathBuf,
    pub annotations: PathBuf,
    pub trait_names: Vec<String>,
    pub covariates: Vec<String>,
    pub mask_categories: Vec<MaskCategory>,
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub individual: bool,
    pub spa: bool,
    pub ancestry_col: Option<String>,
    /// AI-STAAR ensemble base test count B. Adds 2*(B+1) STAAR runs per gene.
    pub ai_base_tests: usize,
    /// AI-STAAR ensemble weight RNG seed.
    pub ai_seed: u64,
    pub scang_params: ScangParams,
    pub known_loci: Option<PathBuf>,
    pub emit_sumstats: bool,
    pub rebuild_store: bool,
    pub column_map: HashMap<String, String>,
    pub output_dir: PathBuf,
    pub store_dir: PathBuf,
}

impl StaarConfig {
    pub fn results_dir(&self) -> PathBuf {
        self.output_dir
            .join("results")
            .join(&self.trait_names[0])
    }

    pub fn sumstats_dir(&self) -> PathBuf {
        self.output_dir
            .join("sumstats")
            .join(&self.trait_names[0])
    }
}

pub type ResultSet = Vec<(MaskType, Vec<GeneResult>)>;

// ---------------------------------------------------------------------------
// Pipeline
// ---------------------------------------------------------------------------

pub struct StaarPipeline<'a> {
    config: StaarConfig,
    out: &'a dyn Output,
    res: Resources,
}

struct ScoringContext {
    analysis: AnalysisVectors,
    use_spa: bool,
    cache_dir: PathBuf,
    ancestry: Option<staar::ancestry::AncestryInfo>,
}

impl<'a> StaarPipeline<'a> {
    pub fn new(config: StaarConfig, out: &'a dyn Output) -> Result<Self, FavorError> {
        let res = store::setup_resources(out)?;
        Ok(Self { config, out, res })
    }

    pub fn run(self) -> Result<(), FavorError> {
        std::fs::create_dir_all(&self.config.output_dir)
            .map_err(|e| FavorError::Resource(format!(
                "Cannot create output directory '{}': {e}", self.config.output_dir.display()
            )))?;

        // Upfront validation: tier check (typed) + column check (string)
        if self.config.annotations.exists() {
            // Typed tier check: STAAR needs all 11 annotation weight columns
            if let Ok(annotated) = AnnotatedSet::open(&self.config.annotations) {
                let weight_cols: Vec<crate::column::Col> = STAAR_WEIGHTS.to_vec();
                annotated.supports(&weight_cols)?;
            }

            // String-based column check: validates raw annotation struct columns
            let ann_vs = VariantSet::open(&self.config.annotations)?;
            let contract = ColumnContract {
                command: "staar",
                required: STAAR_ANNOTATION_COLUMNS,
            };
            let missing = contract.check(ann_vs.columns());
            if !missing.is_empty() {
                let tier_hint = match ann_vs.kind() {
                    Some(VariantSetKind::Annotated { tier: crate::config::Tier::Base }) => {
                        " Your data was annotated with base tier. Re-run: `favor annotate --full`."
                    }
                    _ => " Re-run: `favor annotate --full`.",
                };
                return Err(FavorError::DataMissing(format!(
                    "Missing annotation columns in {}:\n{}\n\
                     STAAR requires favor-full annotations.{}",
                    self.config.annotations.display(),
                    ColumnContract::format_missing(&missing),
                    tier_hint,
                )));
            }
        }

        let store = self.ensure_store()?;
        // Load variants via VariantIndex (same parquet the scoring path uses)
        // then convert to AnnotatedVariant for sumstats/individual/result output.
        let variants = load_rare_variants(
            &store.store_dir, &store.manifest, self.config.maf_cutoff,
        )?;
        let n_rare = variants.len() as i64;
        let primary_trait = &self.config.trait_names[0];

        // Phenotype loading
        let geno_for_pheno = staar::genotype::GenotypeResult {
            sample_names: store.sample_names.clone(),
            output_dir: store
                .geno_output_dir
                .clone()
                .unwrap_or_else(|| store.store_dir.clone()),
        };
        let pheno = load_phenotype(
            &store.engine,
            &self.config.phenotype,
            &self.config.covariates,
            &geno_for_pheno,
            primary_trait,
            self.config.ancestry_col.as_deref(),
            self.config.ai_base_tests,
            self.config.ai_seed,
            &self.config.column_map,
            self.out,
        )?;
        let (y, mut x, trait_type, n, pheno_mask, ancestry) = (
            pheno.y,
            pheno.x,
            pheno.trait_type,
            pheno.n,
            pheno.pheno_mask,
            pheno.ancestry,
        );

        if let Some(ref loci_path) = self.config.known_loci {
            let x_cond = load_known_loci(&store.engine, &geno_for_pheno, loci_path, n, self.out)?;
            x = augment_covariates(&x, &x_cond);
            self.out.status(&format!(
                "  Conditional: {} known loci added as covariates",
                x_cond.ncols()
            ));
        }

        let use_spa = self.config.spa && trait_type == TraitType::Binary;
        if self.config.spa && trait_type == TraitType::Continuous {
            self.out
                .warn("--spa ignored: saddlepoint approximation only applies to binary traits");
        }
        if use_spa {
            self.out
                .status("  SPA enabled: saddlepoint approximation for Burden and ACAT-V");
        }

        // Step 2/4: Fit null model
        let null_model = self.fit_null_model(&y, &x, trait_type)?;

        // Sumstats export (early exit)
        if self.config.emit_sumstats {
            let sumstats_dir = self.config.sumstats_dir();
            std::fs::create_dir_all(&sumstats_dir)
                .map_err(|e| FavorError::Resource(format!(
                    "Cannot create sumstats directory '{}': {e}", sumstats_dir.display()
                )))?;
            let meta = staar::meta::StudyMeta {
                favor_meta_version: 1,
                trait_type: format!("{trait_type:?}"),
                trait_name: self.config.trait_names[0].clone(),
                n_samples: n,
                sigma2: null_model.sigma2,
                maf_cutoff: self.config.maf_cutoff,
                covariates: self.config.covariates.clone(),
                segment_size: 500_000,
            };
            let analysis = AnalysisVectors::from_null_model(&null_model, &pheno_mask);
            return staar::meta::emit_sumstats(
                &store.store_dir,
                &analysis,
                &variants,
                &sumstats_dir,
                &meta,
                self.out,
            );
        }

        // Layer 2: Build or load per-phenotype score cache (U/K for all genes)
        let analysis = AnalysisVectors::from_null_model(&null_model, &pheno_mask);
        let sc_dir = ensure_score_cache(
            &store.store_dir,
            &store.manifest,
            &analysis,
            &store.manifest.key,
            primary_trait,
            &self.config.covariates,
            self.config.known_loci.as_deref(),
            self.out,
        )?;

        let ctx = ScoringContext {
            analysis,
            use_spa,
            cache_dir: sc_dir,
            ancestry,
        };

        // Layer 3: Score tests from cache
        let (results, individual_pvals) = self.score_all(&store, &ctx)?;

        // Step 4/4: Write results
        let results_dir = self.config.results_dir();
        std::fs::create_dir_all(&results_dir)
            .map_err(|e| FavorError::Resource(format!(
                "Cannot create results directory '{}': {e}", results_dir.display()
            )))?;

        if self.config.individual && !individual_pvals.is_empty() {
            write_individual_results(
                &individual_pvals,
                &variants,
                &results_dir,
                self.out,
            )?;
        }

        write_results(
            &results,
            &self.config.trait_names,
            self.config.maf_cutoff,
            &results_dir,
            &null_model,
            trait_type,
            n,
            n_rare,
            self.out,
        )?;

        Ok(())
    }

    // -----------------------------------------------------------------------
    // Pipeline steps
    // -----------------------------------------------------------------------

    pub fn ensure_store(&self) -> Result<GenoStoreResult, FavorError> {
        let geno_staging_dir = self.config.output_dir.join(".geno_staging");
        store::build_or_load_store(
            &self.config.genotypes,
            &self.config.annotations,
            &self.config.store_dir,
            &geno_staging_dir,
            self.config.rebuild_store,
            &self.res,
            self.out,
        )
    }

    pub fn fit_null_model(
        &self,
        y: &Mat<f64>,
        x: &Mat<f64>,
        trait_type: TraitType,
    ) -> Result<staar::model::NullModel, FavorError> {
        self.out.status("Step 2/4: Fitting null model...");
        let null_model = match trait_type {
            TraitType::Continuous => staar::model::fit_glm(y, x),
            TraitType::Binary => staar::model::fit_logistic(y, x, 25),
        };
        self.out
            .status(&format!("  sigma2 = {:.4}", null_model.sigma2));
        Ok(null_model)
    }

}

// ---------------------------------------------------------------------------
// Layer 2: Per-phenotype score cache
// ---------------------------------------------------------------------------

#[allow(clippy::too_many_arguments)]
fn ensure_score_cache(
    store_dir: &Path,
    manifest: &store::StoreManifest,
    analysis: &AnalysisVectors,
    store_key: &str,
    trait_name: &str,
    covariates: &[String],
    known_loci: Option<&Path>,
    out: &dyn Output,
) -> Result<PathBuf, FavorError> {
    let key = score_cache::cache_key(store_key, trait_name, covariates, known_loci);
    let dir = score_cache::cache_dir(store_dir, &key);

    if score_cache::probe(store_dir, manifest, &key) {
        out.status("  Score cache: hit (reusing cached U/K)");
        return Ok(dir);
    }

    out.status("Step 2b/4: Building score cache (all U/K, no MAF filter)...");
    std::fs::create_dir_all(&dir)
        .map_err(|e| FavorError::Resource(format!(
            "Cannot create score cache directory '{}': {e}", dir.display()
        )))?;

    for ci in &manifest.chromosomes {
        let chrom_dir = store_dir.join(format!("chromosome={}", ci.name));
        let sg = SparseG::open(&chrom_dir)?;
        let vi = VariantIndex::load(&chrom_dir)?;
        score_cache::build_chromosome(&sg, &vi, analysis, &dir, &ci.name, out)?;
    }

    Ok(dir)
}

// ---------------------------------------------------------------------------
// Pipeline (continued)
// ---------------------------------------------------------------------------

impl<'a> StaarPipeline<'a> {
    fn score_all(
        &self,
        store: &GenoStoreResult,
        ctx: &ScoringContext,
    ) -> Result<(ResultSet, Vec<(usize, f64)>), FavorError> {
        run_score_tests(
            &store.store_dir,
            &store.manifest,
            &self.config,
            &ctx.analysis,
            ctx,
            self.out,
        )
    }
}

// ---------------------------------------------------------------------------
// Variant loading — single path through VariantIndex
// ---------------------------------------------------------------------------

use crate::staar::carrier::VariantIndex;

/// Load rare variants from the store via VariantIndex.
fn load_rare_variants(
    store_dir: &Path,
    manifest: &store::StoreManifest,
    maf_cutoff: f64,
) -> Result<Vec<AnnotatedVariant>, FavorError> {
    let mut all = Vec::with_capacity(manifest.n_variants);
    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().unwrap_or(Chromosome::Autosome(1));
        let chrom_dir = store_dir.join(format!("chromosome={}", ci.name));
        let index = VariantIndex::load(&chrom_dir)?;
        for entry in index.all_entries() {
            if entry.maf < maf_cutoff {
                all.push(entry.to_annotated_variant(chrom));
            }
        }
    }
    Ok(all)
}

// ---------------------------------------------------------------------------
// Scoring — sparse G + aligned variant vectors
// ---------------------------------------------------------------------------

type MaskPredicate = fn(&AnnotatedVariant) -> bool;

fn run_score_tests(
    store_dir: &Path,
    manifest: &store::StoreManifest,
    config: &StaarConfig,
    analysis: &AnalysisVectors,
    ctx: &ScoringContext,
    out: &dyn Output,
) -> Result<(ResultSet, Vec<(usize, f64)>), FavorError> {
    out.status("Step 3/4: Running score tests (carrier-indexed sparse)...");

    let mut mask_predicates: Vec<(MaskType, MaskPredicate)> = Vec::new();
    let mut has_windows = false;
    let mut has_scang = false;
    for cat in &config.mask_categories {
        match cat {
            MaskCategory::Coding => {
                for &(ref mt, pred) in masks::CODING_MASKS {
                    mask_predicates.push((mt.clone(), pred));
                }
            }
            MaskCategory::Noncoding => {
                for &(ref mt, pred) in masks::NONCODING_MASKS {
                    mask_predicates.push((mt.clone(), pred));
                }
            }
            MaskCategory::SlidingWindow => has_windows = true,
            MaskCategory::Scang => has_scang = true,
            MaskCategory::Custom => out.warn("Custom BED: not yet implemented"),
        }
    }

    let mut all_results: ResultSet = mask_predicates
        .iter()
        .map(|(mt, _)| (mt.clone(), Vec::new()))
        .collect();

    let window_idx = if has_windows {
        all_results.push((MaskType::SlidingWindow, Vec::new()));
        Some(all_results.len() - 1)
    } else {
        None
    };

    let scang_idx = if has_scang {
        all_results.push((MaskType::Scang, Vec::new()));
        Some(all_results.len() - 1)
    } else {
        None
    };

    let individual_pvals: Vec<(usize, f64)> = Vec::new();
    let maf_cutoff = config.maf_cutoff;

    for ci in &manifest.chromosomes {
        let chrom_name = &ci.name;
        let chrom_parsed: Chromosome = chrom_name.parse().unwrap_or(Chromosome::Autosome(1));
        let chrom_dir = store_dir.join(format!("chromosome={chrom_name}"));

        let sparse_g = SparseG::open(&chrom_dir)?;
        let variant_index = VariantIndex::load(&chrom_dir)?;

        out.status(&format!(
            "  chr{chrom_name}: {} variants, {} genes",
            variant_index.len(),
            variant_index.n_genes(),
        ));

        // Layer 3: Gene masks — score from cached U/K
        if !mask_predicates.is_empty() {
            let cache = score_cache::load_chromosome(&ctx.cache_dir, chrom_name, &variant_index)?;
            let all_entries = variant_index.all_entries();
            let use_spa = ctx.use_spa;
            let ancestry = ctx.ancestry.as_ref();
            let n_vcf = analysis.n_vcf_total;

            let gene_names: Vec<&String> = cache.gene_blocks.keys().collect();
            let per_gene_results: Vec<Vec<(usize, GeneResult)>> = gene_names
                .par_iter()
                .filter_map(|gene_name| {
                    let block = cache.gene_blocks.get(*gene_name)?;
                    if block.m() < 2 { return None; }

                    // MAF filter
                    let maf_pass: Vec<usize> = (0..block.m())
                        .filter(|&local| {
                            let global = block.variant_offsets[local] as usize;
                            all_entries[global].maf < maf_cutoff
                        })
                        .collect();
                    if maf_pass.len() < 2 { return None; }

                    let start = maf_pass.iter()
                        .map(|&l| all_entries[block.variant_offsets[l] as usize].position)
                        .min().unwrap();
                    let end = maf_pass.iter()
                        .map(|&l| all_entries[block.variant_offsets[l] as usize].position)
                        .max().unwrap();

                    // Large-gene fallback: K not cached, load from SparseG
                    if !block.has_k() {
                        let gene_vcfs: Vec<u32> = maf_pass.iter()
                            .map(|&l| block.variant_offsets[l])
                            .collect();
                        let carriers = sparse_g.load_variants(&gene_vcfs);
                        if carriers.len() < 2 { return None; }
                        let (full_u, full_k) = sparse_score::score_gene_sparse(&carriers, analysis);

                        let mut results = Vec::new();
                        for (mask_idx, (_mt, predicate)) in mask_predicates.iter().enumerate() {
                            let qualifying: Vec<usize> = (0..gene_vcfs.len())
                                .filter(|&i| {
                                    let v = gene_vcfs[i] as usize;
                                    predicate(&all_entries[v].to_annotated_variant_with_gene(chrom_parsed, gene_name))
                                })
                                .collect();
                            if qualifying.len() < 2 { continue; }

                            let mafs: Vec<f64> = qualifying.iter().map(|&i| all_entries[gene_vcfs[i] as usize].maf).collect();
                            let ann_matrix: Vec<Vec<f64>> = (0..11)
                                .map(|ch| qualifying.iter().map(|&i| all_entries[gene_vcfs[i] as usize].weights.0[ch]).collect())
                                .collect();

                            let sr = if let Some(ai) = ancestry {
                                let subset: Vec<_> = qualifying.iter().map(|&i| carriers[i].clone()).collect();
                                staar::ancestry::run_ai_staar_gene(&subset, analysis, ai, &ann_matrix, use_spa)
                            } else if use_spa {
                                let subset: Vec<_> = qualifying.iter().map(|&i| carriers[i].clone()).collect();
                                sparse_score::run_staar_sparse(&subset, analysis, &ann_matrix, &mafs, true)
                            } else {
                                let (u_sub, k_sub) = sparse_score::slice_sumstats(&full_u, &full_k, &qualifying);
                                score::run_staar_from_sumstats(&u_sub, &k_sub, &ann_matrix, &mafs, analysis.n_pheno)
                            };

                            let cmac: u32 = mafs.iter().map(|&maf| (2.0 * maf * n_vcf as f64).round() as u32).sum();
                            results.push((mask_idx, GeneResult {
                                ensembl_id: gene_name.to_string(),
                                gene_symbol: gene_name.to_string(),
                                chromosome: chrom_parsed,
                                start, end,
                                n_variants: qualifying.len() as u32,
                                cumulative_mac: cmac,
                                staar: sr,
                            }));
                        }
                        return if results.is_empty() { None } else { Some(results) };
                    }

                    // Standard cached path: slice U/K for each mask
                    let mut results = Vec::new();
                    for (mask_idx, (_mt, predicate)) in mask_predicates.iter().enumerate() {
                        let qualifying: Vec<usize> = maf_pass.iter()
                            .filter(|&&local| {
                                let global = block.variant_offsets[local] as usize;
                                predicate(&all_entries[global].to_annotated_variant_with_gene(chrom_parsed, gene_name))
                            })
                            .copied()
                            .collect();
                        if qualifying.len() < 2 { continue; }

                        let mafs: Vec<f64> = qualifying.iter()
                            .map(|&l| all_entries[block.variant_offsets[l] as usize].maf)
                            .collect();
                        let ann_matrix: Vec<Vec<f64>> = (0..11)
                            .map(|ch| qualifying.iter()
                                .map(|&l| all_entries[block.variant_offsets[l] as usize].weights.0[ch])
                                .collect())
                            .collect();

                        let sr = if let Some(ai) = ancestry {
                            let mask_vcfs: Vec<u32> = qualifying.iter()
                                .map(|&l| block.variant_offsets[l])
                                .collect();
                            let subset = sparse_g.load_variants(&mask_vcfs);
                            staar::ancestry::run_ai_staar_gene(&subset, analysis, ai, &ann_matrix, use_spa)
                        } else if use_spa {
                            // SPA: load carriers from SparseG for this mask
                            let mask_vcfs: Vec<u32> = qualifying.iter()
                                .map(|&l| block.variant_offsets[l])
                                .collect();
                            let subset = sparse_g.load_variants(&mask_vcfs);
                            sparse_score::run_staar_sparse(&subset, analysis, &ann_matrix, &mafs, true)
                        } else {
                            let global_indices: Vec<usize> = qualifying.iter()
                                .map(|&l| block.variant_offsets[l] as usize)
                                .collect();
                            let u_sub = score_cache::slice_window_u(&cache, &global_indices);
                            let k_sub = sparse_score::slice_sumstats_flat(
                                &block.k_flat, block.m(), &qualifying,
                            );
                            score::run_staar_from_sumstats(
                                &u_sub, &k_sub, &ann_matrix, &mafs, analysis.n_pheno,
                            )
                        };

                        let cmac: u32 = qualifying.iter()
                            .map(|&l| {
                                let maf = all_entries[block.variant_offsets[l] as usize].maf;
                                (2.0 * maf * n_vcf as f64).round() as u32
                            })
                            .sum();

                        results.push((mask_idx, GeneResult {
                            ensembl_id: gene_name.to_string(),
                            gene_symbol: gene_name.to_string(),
                            chromosome: chrom_parsed,
                            start, end,
                            n_variants: qualifying.len() as u32,
                            cumulative_mac: cmac,
                            staar: sr,
                        }));
                    }

                    if results.is_empty() { None } else { Some(results) }
                })
                .collect();

            for gene_results in per_gene_results {
                for (mask_idx, result) in gene_results {
                    all_results[mask_idx].1.push(result);
                }
            }

            for (idx, (mask_type, _)) in mask_predicates.iter().enumerate() {
                let n = all_results[idx].1.len();
                if n > 0 {
                    out.status(&format!("    {}: {} groups", mask_type.file_stem(), n));
                }
            }
        }

        // Sliding windows and SCANG
        if window_idx.is_some() || scang_idx.is_some() {
            let (chrom_variants, global_indices): (Vec<AnnotatedVariant>, Vec<usize>) =
                variant_index.all_entries().iter().enumerate()
                    .filter(|(_, e)| e.maf < maf_cutoff)
                    .map(|(gi, e)| (e.to_annotated_variant(chrom_parsed), gi))
                    .unzip();
            let chrom_indices: Vec<usize> = (0..chrom_variants.len()).collect();

            let window_cache = score_cache::load_chromosome(
                &ctx.cache_dir, chrom_name, &variant_index,
            )?;
            let n_vcf = analysis.n_vcf_total;
            let ancestry = ctx.ancestry.as_ref();

            let score_window = |g: &MaskGroup| -> Option<GeneResult> {
                let m = g.variant_indices.len();
                if m < 2 { return None; }

                let win_globals: Vec<usize> = g.variant_indices.iter()
                    .map(|&ci| global_indices[ci])
                    .collect();

                let mafs: Vec<f64> = g.variant_indices.iter()
                    .map(|&ci| chrom_variants[ci].maf)
                    .collect();
                let ann_matrix: Vec<Vec<f64>> = (0..11)
                    .map(|ch| g.variant_indices.iter()
                        .map(|&ci| chrom_variants[ci].annotation.weights.0[ch])
                        .collect())
                    .collect();

                let sr = if let Some(ai) = ancestry {
                    let win_vcfs: Vec<u32> = win_globals.iter().map(|&i| i as u32).collect();
                    let carriers = sparse_g.load_variants(&win_vcfs);
                    staar::ancestry::run_ai_staar_gene(&carriers, analysis, ai, &ann_matrix, ctx.use_spa)
                } else {
                    let u_win = score_cache::slice_window_u(&window_cache, &win_globals);
                    let k_win = score_cache::assemble_window_k(&window_cache, &win_globals);
                    score::run_staar_from_sumstats(
                        &u_win, &k_win, &ann_matrix, &mafs, analysis.n_pheno,
                    )
                };
                let cmac: u32 = mafs.iter()
                    .map(|&maf| (2.0 * maf * n_vcf as f64).round() as u32)
                    .sum();

                Some(GeneResult {
                    ensembl_id: g.name.clone(),
                    gene_symbol: g.name.clone(),
                    chromosome: g.chromosome,
                    start: g.start,
                    end: g.end,
                    n_variants: m as u32,
                    cumulative_mac: cmac,
                    staar: sr,
                })
            };

            if let Some(wi) = window_idx {
                let window_groups = masks::build_sliding_windows(
                    &chrom_variants, &chrom_indices, chrom_parsed,
                    config.window_size, config.window_size / 2,
                );
                if !window_groups.is_empty() {
                    let r: Vec<GeneResult> = window_groups.iter()
                        .filter_map(score_window).collect();
                    if !r.is_empty() {
                        out.status(&format!("    sliding_window: {} windows", r.len()));
                        all_results[wi].1.extend(r);
                    }
                }
            }

            if let Some(si) = scang_idx {
                let scang_all = masks::build_scang_windows(
                    &chrom_variants, &chrom_indices, chrom_parsed,
                    &config.scang_params,
                );
                for (wsize, groups) in &scang_all {
                    let r: Vec<GeneResult> = groups.iter()
                        .filter_map(score_window).collect();
                    if !r.is_empty() {
                        out.status(&format!("    scang L={wsize}: {} windows", r.len()));
                        all_results[si].1.extend(r);
                    }
                }
            }
        }
    }

    Ok((all_results, individual_pvals))
}

