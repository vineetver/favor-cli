//! STAAR analysis pipeline orchestration.
//!
//! Owns `StaarConfig`, `StaarPipeline`, the per-run manifest, and the typed
//! stage functions. Per-gene and per-window scoring live in `scoring.rs`.
//!
//! Stage layout (each is a small function on `StaarPipeline`):
//! ```text
//! validate → ensure_store → load_phenotype → fit_null_model
//!   → ensure_score_cache → run_scoring → write_results
//! ```
//! After each stage the manifest is rewritten atomically so an interrupted
//! run can be inspected to see exactly how far it got.

use std::collections::HashMap;
use std::path::PathBuf;

use faer::Mat;

use crate::column::STAAR_WEIGHTS;
use crate::error::CohortError;
use crate::ingest::ColumnContract;
use crate::output::{bail_if_cancelled, Output};
use crate::resource::Resources;
use crate::staar::carrier::AnalysisVectors;
use crate::staar::masks::ScangParams;
use crate::staar::model::{augment_covariates, load_known_loci, load_phenotype, NullModel};
use crate::staar::output::{write_individual_results, write_results};
use crate::staar::run_manifest::{
    self, ArtifactKind, CacheDecision, CacheOutcome, ConfigHashInputs, ResumeDecision,
    ResumeSummary, RunManifest, Stage,
};
use crate::staar::scoring::{self, ResultSet};
use crate::staar::{
    self, ancestry::AncestryInfo, MaskCategory, NullModelKind, RunMode, ScoringMode, TraitType,
};
use crate::store::cache::score_cache;
use crate::store::cohort::{
    self, CohortHandle, CohortId, CohortManifest, GenoStoreResult, STAAR_ANNOTATION_COLUMNS,
};
use crate::store::list::{AnnotatedSet, VariantSet, VariantSetKind};
use crate::store::Store;
use crate::types::{AnnotatedVariant, Chromosome};

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
    /// AI-STAAR ensemble base test count B.
    pub ai_base_tests: usize,
    /// AI-STAAR ensemble RNG seed.
    pub ai_seed: u64,
    pub scang_params: ScangParams,
    /// Kinship matrix files for mixed-model analysis. Empty = no kinship.
    pub kinship: Vec<PathBuf>,
    /// Phenotype column for heteroscedastic residual variance partition.
    pub kinship_groups: Option<String>,
    pub known_loci: Option<PathBuf>,
    pub emit_sumstats: bool,
    pub rebuild_store: bool,
    pub column_map: HashMap<String, String>,
    pub output_dir: PathBuf,
    pub cohort_id: CohortId,
}

impl StaarConfig {
    pub fn results_dir(&self) -> PathBuf {
        self.output_dir.join("results").join(&self.trait_names[0])
    }

    pub fn sumstats_dir(&self) -> PathBuf {
        self.output_dir.join("sumstats").join(&self.trait_names[0])
    }

    /// What this run actually does end-to-end.
    pub fn run_mode(&self) -> RunMode {
        if self.emit_sumstats {
            RunMode::EmitSumstats
        } else {
            RunMode::Analyze
        }
    }

    /// Choose the scoring backend from the trait type and CLI flags.
    ///
    /// Trait type is required because SPA only applies to binary traits;
    /// for continuous traits with `--spa` we drop back to Standard and warn
    /// at the call site.
    pub fn scoring_mode(&self, trait_type: TraitType) -> ScoringMode {
        if self.ancestry_col.is_some() {
            ScoringMode::AiStaar
        } else if self.spa && trait_type == TraitType::Binary {
            ScoringMode::Spa
        } else {
            ScoringMode::Standard
        }
    }

    /// True if this run uses any kinship-aware path (REML or PQL).
    /// Single source of truth for the predicate, used by `null_model_kind`
    /// and the SPA-vs-kinship warning in `warn_mode_combinations`.
    pub fn has_kinship(&self) -> bool {
        !self.kinship.is_empty() || self.kinship_groups.is_some()
    }

    /// Pick the null-model fitter from trait type + kinship presence.
    pub fn null_model_kind(&self, trait_type: TraitType) -> NullModelKind {
        match (trait_type, self.has_kinship()) {
            (TraitType::Continuous, false) => NullModelKind::Glm,
            (TraitType::Binary, false) => NullModelKind::Logistic,
            (TraitType::Continuous, true) => NullModelKind::KinshipReml,
            (TraitType::Binary, true) => NullModelKind::KinshipPql,
        }
    }
}

impl StaarConfig {
    /// Borrowed view fed into `run_manifest::compute_config_hash`.
    pub fn hash_inputs(&self) -> ConfigHashInputs<'_> {
        ConfigHashInputs {
            genotypes: &self.genotypes,
            annotations: &self.annotations,
            phenotype: &self.phenotype,
            trait_names: &self.trait_names,
            covariates: &self.covariates,
            maf_cutoff: self.maf_cutoff,
            run_mode: self.run_mode(),
            kinship: &self.kinship,
            kinship_groups: self.kinship_groups.as_deref(),
        }
    }

    pub fn config_hash(&self) -> String {
        run_manifest::compute_config_hash(&self.hash_inputs())
    }

    pub fn new_run_manifest(&self) -> RunManifest {
        RunManifest::new(
            self.run_mode(),
            self.trait_names[0].clone(),
            self.config_hash(),
        )
    }
}

/// Output of `stage_run_scoring`. `individual_pvals` is reserved for the
/// `--individual` per-variant path; today the producer in
/// `scoring::run_score_tests` returns `Vec::new()`. Field stays so a
/// future producer doesn't have to ripple a new tuple shape up the stack.
pub struct ScoringOutput {
    pub results: ResultSet,
    pub individual_pvals: Vec<(usize, f64)>,
}

/// Per-run scoring context: fixed once after the null model is fit.
pub struct ScoringContext {
    pub analysis: AnalysisVectors,
    pub scoring_mode: ScoringMode,
    pub cache_dir: PathBuf,
    pub ancestry: Option<AncestryInfo>,
    /// True iff `--spa` was requested AND the trait is binary. Recorded once
    /// at construction so the `--ancestry-col --spa --binary-trait` triple
    /// keeps SPA inside the per-population AI-STAAR kernel without
    /// recomputing the predicate at every gene.
    pub spa_active: bool,
}

/// Output of `load_phenotype` stage.
struct PhenoStageOut {
    y: Mat<f64>,
    x: Mat<f64>,
    trait_type: TraitType,
    n: usize,
    pheno_mask: Vec<bool>,
    ancestry: Option<AncestryInfo>,
    /// Compact sample list aligned to `pheno_mask` (only samples with both
    /// phenotype and genotype). Used by kinship loaders.
    compact_samples: Vec<String>,
    /// Genotype-result view rebuilt once and shared with kinship loading.
    /// Lives here so `stage_fit_null_model` doesn't have to reach back into
    /// `GenoStoreResult` and reconstruct one of these by hand.
    genotype_result: staar::genotype::GenotypeResult,
}

pub struct StaarPipeline<'a> {
    config: StaarConfig,
    store: Store,
    out: &'a dyn Output,
    res: Resources,
    manifest: RunManifest,
    resume: ResumeDecision,
}

impl<'a> StaarPipeline<'a> {
    pub fn new(config: StaarConfig, store: Store, out: &'a dyn Output) -> Result<Self, CohortError> {
        let res = cohort::setup_resources(out)?;
        let manifest = config.new_run_manifest();
        Ok(Self {
            config,
            store,
            out,
            res,
            manifest,
            resume: ResumeDecision::Fresh,
        })
    }

    /// Stage runner. Each stage gates on the previous one returning Ok and
    /// rewrites `run.json` after every transition.
    pub fn run(mut self) -> Result<(), CohortError> {
        std::fs::create_dir_all(&self.config.output_dir).map_err(|e| {
            CohortError::Resource(format!(
                "Cannot create output directory '{}': {e}",
                self.config.output_dir.display()
            ))
        })?;

        let prior = RunManifest::probe(&self.config.output_dir);
        self.resume = run_manifest::plan_resume(prior, &self.manifest.config_hash);
        self.manifest.outputs.resume_summary = Some(ResumeSummary::from(&self.resume));
        match &self.resume {
            ResumeDecision::Fresh => {}
            ResumeDecision::Discarded(reason) => {
                self.out
                    .warn(&format!("  prior run.json discarded: {reason}"));
            }
            ResumeDecision::Resume { skippable } if skippable.is_empty() => {
                self.out
                    .status("  prior run.json found (same config); no stages eligible to skip");
            }
            ResumeDecision::Resume { skippable } => {
                let labels: Vec<&str> = skippable.iter().map(|s| s.label()).collect();
                self.out.status(&format!(
                    "  resuming: skipping {} (durable on disk)",
                    labels.join(", ")
                ));
            }
        }

        let result = self.run_stages();
        if result.is_err() {
            // Keep the manifest as a forensic record on failure.
            let _ = self.manifest.write(&self.config.output_dir);
        }
        result
    }

    fn run_stages(&mut self) -> Result<(), CohortError> {
        self.stage(Stage::Validate, |p| p.stage_validate())?;

        let store = self.stage(Stage::EnsureStore, |p| p.stage_ensure_store())?;
        self.manifest.outputs.store_dir = Some(store.store_dir.clone());
        self.manifest.write(&self.config.output_dir)?;

        let pheno = self.stage(Stage::LoadPhenotype, |p| p.stage_load_phenotype(&store))?;

        let null_model = self.stage(Stage::FitNullModel, |p| {
            p.stage_fit_null_model(&pheno, &store)
        })?;

        // Sumstats early-exit run mode.
        if self.config.run_mode() == RunMode::EmitSumstats {
            self.stage(Stage::EmitSumstats, |p| {
                p.stage_emit_sumstats(&store, &null_model, &pheno)
            })?;
            self.manifest.write(&self.config.output_dir)?;
            return Ok(());
        }

        let analysis = AnalysisVectors::from_null_model(&null_model, &pheno.pheno_mask)?;

        let cache_dir = self.stage(Stage::EnsureScoreCache, |p| {
            p.stage_ensure_score_cache(&store, &analysis)
        })?;
        self.manifest.outputs.score_cache_dir = Some(cache_dir.clone());
        self.manifest.write(&self.config.output_dir)?;

        let ctx = ScoringContext {
            analysis,
            scoring_mode: self.config.scoring_mode(pheno.trait_type),
            cache_dir,
            ancestry: pheno.ancestry.clone(),
            spa_active: self.config.spa && pheno.trait_type == TraitType::Binary,
        };
        self.warn_mode_combinations(pheno.trait_type, ctx.scoring_mode);

        let scoring =
            self.stage(Stage::RunScoring, |p| p.stage_run_scoring(&store, &ctx))?;

        // load_rare_variants walks every chromosome and reads each VariantIndex
        // through the cohort handle. Hold the result in a local so the
        // WriteResults stage doesn't trigger a second walk.
        let cohort = self.store.cohort(&self.config.cohort_id);
        let variants = load_rare_variants(&cohort, &store.manifest, self.config.maf_cutoff)?;
        let n_rare = variants.len() as i64;

        self.stage(Stage::WriteResults, |p| {
            p.stage_write_results(
                &scoring,
                &variants,
                &null_model,
                pheno.trait_type,
                pheno.n,
                n_rare,
            )
        })?;
        self.manifest.outputs.results_dir = Some(self.config.results_dir());
        self.manifest.write(&self.config.output_dir)?;

        Ok(())
    }

    /// Run one stage with begin/complete/fail bookkeeping and a manifest
    /// rewrite at every transition.
    ///
    /// Manifest write semantics:
    ///   - `begin` write is best-effort (an in-progress marker isn't worth
    ///     killing a long-running stage over).
    ///   - `complete` write is required: if the operator can't trust the
    ///     manifest, resume becomes impossible. Surfaced as a `Resource`
    ///     error to the caller.
    ///   - `fail` write is best-effort: we already have a real error to
    ///     return; a secondary fsync failure shouldn't mask it.
    fn stage<T, F>(&mut self, stage: Stage, body: F) -> Result<T, CohortError>
    where
        F: FnOnce(&mut Self) -> Result<T, CohortError>,
    {
        bail_if_cancelled(self.out)?;
        self.manifest.begin(stage);
        let _ = self.manifest.write(&self.config.output_dir);
        match body(self) {
            Ok(v) => {
                self.manifest.complete(stage);
                self.manifest.write(&self.config.output_dir)?;
                bail_if_cancelled(self.out)?;
                Ok(v)
            }
            Err(e) => {
                self.manifest.fail(stage);
                let _ = self.manifest.write(&self.config.output_dir);
                Err(e)
            }
        }
    }

    fn stage_validate(&mut self) -> Result<(), CohortError> {
        // `build_config` already verified the annotation path exists; if it
        // disappeared between then and now we want a hard error, not silent
        // skip — surface the missing-file message rather than no-op.
        if !self.config.annotations.exists() {
            return Err(CohortError::Input(format!(
                "Annotations no longer at '{}'. Re-run `cohort annotate` or check the path.",
                self.config.annotations.display()
            )));
        }

        // Typed tier check: STAAR needs all 11 annotation weight columns.
        if let Ok(annotated) = AnnotatedSet::open(&self.config.annotations) {
            let weight_cols: Vec<crate::column::Col> = STAAR_WEIGHTS.to_vec();
            annotated.supports(&weight_cols)?;
        }

        // Raw annotation column contract.
        let ann_vs = VariantSet::open(&self.config.annotations)?;
        let contract = ColumnContract {
            command: "staar",
            required: STAAR_ANNOTATION_COLUMNS,
        };
        let missing = contract.check(ann_vs.columns());
        if missing.is_empty() {
            return Ok(());
        }
        let tier_hint = match ann_vs.kind() {
            Some(VariantSetKind::Annotated {
                tier: crate::config::Tier::Base,
            }) => " Your data was annotated with base tier. Re-run: `cohort annotate --full`.",
            _ => " Re-run: `cohort annotate --full`.",
        };
        Err(CohortError::DataMissing(format!(
            "Missing annotation columns in {}:\n{}\n\
             STAAR requires FAVOR full-tier annotations.{}",
            self.config.annotations.display(),
            ColumnContract::format_missing(&missing),
            tier_hint,
        )))
    }

    fn stage_ensure_store(&mut self) -> Result<GenoStoreResult, CohortError> {
        // One probe — the recorded cache decision and the path actually
        // taken in `build_or_load` are sourced from the same StoreProbe
        // so they cannot drift.
        let cohort = self.store.cohort(&self.config.cohort_id);
        let probe = if self.config.rebuild_store {
            self.manifest.outputs.cache_decisions.push(CacheDecision {
                artifact: ArtifactKind::GenotypeStore,
                outcome: CacheOutcome::Rebuilt,
                reason: "rebuild requested by --rebuild-store".into(),
            });
            // Drop any score-cache + lookup index that referenced this
            // cohort under its current id; the rebuild changes the
            // content fingerprint and would orphan them otherwise.
            if let Err(e) = self.store.cache().prune_cohort(&self.config.cohort_id) {
                self.out
                    .warn(&format!("  prune cache before rebuild: {e}"));
            }
            cohort::StoreProbe {
                store_dir: cohort.dir().to_path_buf(),
                manifest: None,
                miss_reason: None,
            }
        } else {
            let probe = cohort.probe(&self.config.genotypes, &self.config.annotations);
            let decision = if probe.manifest.is_some() {
                CacheDecision {
                    artifact: ArtifactKind::GenotypeStore,
                    outcome: CacheOutcome::Hit,
                    reason: "content fingerprint matched".into(),
                }
            } else {
                let reason = probe
                    .miss_reason
                    .as_ref()
                    .map(cohort::describe_miss)
                    .unwrap_or_else(|| "unknown miss".into());
                CacheDecision {
                    artifact: ArtifactKind::GenotypeStore,
                    outcome: CacheOutcome::Miss,
                    reason,
                }
            };
            self.manifest.outputs.cache_decisions.push(decision);
            probe
        };

        let geno_staging_dir = self.config.output_dir.join(".geno_staging");
        cohort.build_or_load(
            &self.config.genotypes,
            &self.config.annotations,
            &geno_staging_dir,
            self.config.rebuild_store,
            probe,
            &self.res,
            self.out,
        )
    }

    fn stage_load_phenotype(
        &mut self,
        store: &GenoStoreResult,
    ) -> Result<PhenoStageOut, CohortError> {
        let genotype_result = store.to_genotype_result();
        let primary_trait = &self.config.trait_names[0];

        let pheno = load_phenotype(
            &store.engine,
            &self.config.phenotype,
            &self.config.covariates,
            &genotype_result,
            primary_trait,
            self.config.ancestry_col.as_deref(),
            self.config.ai_base_tests,
            self.config.ai_seed,
            &self.config.column_map,
            self.out,
        )?;
        let mut x = pheno.x;

        if let Some(ref loci_path) = self.config.known_loci {
            let x_cond =
                load_known_loci(&store.engine, &genotype_result, loci_path, pheno.n, self.out)?;
            x = augment_covariates(&x, &x_cond);
            self.out.status(&format!(
                "  Conditional: {} known loci added as covariates",
                x_cond.ncols()
            ));
        }

        let compact_samples: Vec<String> = store
            .sample_names
            .iter()
            .zip(pheno.pheno_mask.iter())
            .filter_map(|(name, &has)| if has { Some(name.clone()) } else { None })
            .collect();

        Ok(PhenoStageOut {
            y: pheno.y,
            x,
            trait_type: pheno.trait_type,
            n: pheno.n,
            pheno_mask: pheno.pheno_mask,
            ancestry: pheno.ancestry,
            compact_samples,
            genotype_result,
        })
    }

    fn stage_fit_null_model(
        &mut self,
        pheno: &PhenoStageOut,
        store: &GenoStoreResult,
    ) -> Result<NullModel, CohortError> {
        self.out.status("Fitting null model...");

        match self.config.null_model_kind(pheno.trait_type) {
            NullModelKind::Glm => {
                let nm = staar::model::fit_glm(&pheno.y, &pheno.x);
                self.out.status(&format!("  sigma2 = {:.4}", nm.sigma2));
                Ok(nm)
            }
            NullModelKind::Logistic => {
                let nm = staar::model::fit_logistic(&pheno.y, &pheno.x, 25);
                self.out.status(&format!("  sigma2 = {:.4}", nm.sigma2));
                Ok(nm)
            }
            kind @ (NullModelKind::KinshipReml | NullModelKind::KinshipPql) => {
                self.fit_kinship_null_model(kind, pheno, store)
            }
        }
    }

    /// Kinship-aware null model: load kinship matrices, build the group
    /// partition (or fall back to a single group), fit REML or PQL, and
    /// translate the resulting `KinshipState` into a `NullModel`.
    fn fit_kinship_null_model(
        &mut self,
        kind: NullModelKind,
        pheno: &PhenoStageOut,
        store: &GenoStoreResult,
    ) -> Result<NullModel, CohortError> {
        let kinships = staar::kinship::load_kinship(
            &self.config.kinship,
            &pheno.compact_samples,
            self.out,
        )?;
        let groups = match self.config.kinship_groups.as_deref() {
            Some(col) => staar::kinship::load_groups(
                &store.engine,
                &self.config.phenotype,
                col,
                &pheno.genotype_result,
                &pheno.pheno_mask,
                &self.config.column_map,
                self.out,
            )?,
            None => staar::kinship::GroupPartition::single(pheno.n),
        };

        let state = match kind {
            NullModelKind::KinshipReml => {
                staar::kinship::fit_reml(&pheno.y, &pheno.x, &kinships, &groups, None)?
            }
            NullModelKind::KinshipPql => staar::kinship::fit_pql_glmm(
                &pheno.y,
                &pheno.x,
                &kinships,
                &groups,
                self.out,
            )?,
            _ => unreachable!("fit_kinship_null_model called with non-kinship kind"),
        };
        self.out.status(&format!(
            "  AI-REML converged in {} iterations ({} boundary refits), h² = {:?}",
            state.n_iter, state.outer_refits, state.h2,
        ));

        let n = pheno.y.nrows();
        let mut residuals = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            residuals[(i, 0)] = state.p_y[(i, 0)];
        }
        Ok(NullModel {
            residuals,
            x_matrix: pheno.x.clone(),
            xtx_inv: state.cov.clone(),
            sigma2: 1.0,
            n_samples: n,
            fitted_values: None,
            working_weights: None,
            kinship: Some(state),
        })
    }

    fn stage_emit_sumstats(
        &mut self,
        store: &GenoStoreResult,
        null_model: &NullModel,
        pheno: &PhenoStageOut,
    ) -> Result<(), CohortError> {
        let sumstats_dir = self.config.sumstats_dir();
        std::fs::create_dir_all(&sumstats_dir).map_err(|e| {
            CohortError::Resource(format!(
                "Cannot create sumstats directory '{}': {e}",
                sumstats_dir.display()
            ))
        })?;
        let cohort = self.store.cohort(&self.config.cohort_id);
        let variants = load_rare_variants(&cohort, &store.manifest, self.config.maf_cutoff)?;
        let analysis = AnalysisVectors::from_null_model(null_model, &pheno.pheno_mask)?;
        let meta = staar::meta::StudyMeta {
            cohort_meta_version: 1,
            trait_type: format!("{:?}", pheno.trait_type),
            trait_name: self.config.trait_names[0].clone(),
            n_samples: pheno.n,
            sigma2: null_model.sigma2,
            maf_cutoff: self.config.maf_cutoff,
            covariates: self.config.covariates.clone(),
            segment_size: 500_000,
        };
        staar::meta::emit_sumstats(
            &self.store,
            &self.config.cohort_id,
            &analysis,
            &variants,
            &sumstats_dir,
            &meta,
            self.out,
        )?;
        self.manifest.outputs.sumstats_dir = Some(sumstats_dir);
        Ok(())
    }

    fn stage_ensure_score_cache(
        &mut self,
        store: &GenoStoreResult,
        analysis: &AnalysisVectors,
    ) -> Result<PathBuf, CohortError> {
        let key = score_cache::cache_key(
            &store.manifest.key,
            &self.config.trait_names[0],
            &self.config.covariates,
            self.config.known_loci.as_deref(),
            &self.config.kinship,
            self.config.kinship_groups.as_deref(),
        );
        let cache_key = crate::store::ids::CacheKey::new(&key);
        let dir = self
            .store
            .cache()
            .score_cache_dir(&self.config.cohort_id, &cache_key);

        match score_cache::probe(&dir, &store.manifest) {
            None => {
                self.out.status("  Score cache: hit (reusing cached U/K)");
                self.manifest.outputs.cache_decisions.push(CacheDecision {
                    artifact: ArtifactKind::ScoreCache,
                    outcome: CacheOutcome::Hit,
                    reason: format!("cache_key={key} matched"),
                });
                return Ok(dir);
            }
            Some(miss) => {
                self.manifest.outputs.cache_decisions.push(CacheDecision {
                    artifact: ArtifactKind::ScoreCache,
                    outcome: CacheOutcome::Miss,
                    reason: format!("cache_key={key}: {}", miss.describe()),
                });
            }
        }

        self.out
            .status("Building score cache (all U/K, no MAF filter)...");
        std::fs::create_dir_all(&dir).map_err(|e| {
            CohortError::Resource(format!("Cannot create score cache directory '{}': {e}", dir.display()))
        })?;

        let cohort = self.store.cohort(&self.config.cohort_id);
        for ci in &store.manifest.chromosomes {
            let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
            let view = cohort.chromosome(&chrom)?;
            score_cache::build_chromosome(
                view.sparse_g()?,
                view.index()?,
                analysis,
                &dir,
                &ci.name,
                self.out,
            )?;
        }

        Ok(dir)
    }

    fn stage_run_scoring(
        &mut self,
        store: &GenoStoreResult,
        ctx: &ScoringContext,
    ) -> Result<ScoringOutput, CohortError> {
        let (results, individual_pvals) = scoring::run_score_tests(
            &self.store,
            &self.config.cohort_id,
            &store.manifest,
            &self.config,
            &ctx.analysis,
            ctx,
            self.out,
        )?;
        Ok(ScoringOutput {
            results,
            individual_pvals,
        })
    }

    fn stage_write_results(
        &mut self,
        scoring: &ScoringOutput,
        variants: &[AnnotatedVariant],
        null_model: &NullModel,
        trait_type: TraitType,
        n: usize,
        n_rare: i64,
    ) -> Result<(), CohortError> {
        let results_dir = self.config.results_dir();
        std::fs::create_dir_all(&results_dir).map_err(|e| {
            CohortError::Resource(format!(
                "Cannot create results directory '{}': {e}",
                results_dir.display()
            ))
        })?;

        if self.config.individual && !scoring.individual_pvals.is_empty() {
            write_individual_results(&scoring.individual_pvals, variants, &results_dir, self.out)?;
        }

        write_results(
            &scoring.results,
            &self.config.trait_names,
            self.config.maf_cutoff,
            &results_dir,
            null_model,
            trait_type,
            n,
            n_rare,
            self.out,
        )
    }

    fn warn_mode_combinations(&self, trait_type: TraitType, mode: ScoringMode) {
        if self.config.spa && trait_type == TraitType::Continuous {
            self.out
                .warn("--spa ignored: saddlepoint approximation only applies to binary traits");
        }
        if mode == ScoringMode::Spa {
            self.out
                .status("  SPA enabled: saddlepoint approximation for Burden and ACAT-V");
        }
        if self.config.has_kinship() && mode == ScoringMode::Spa {
            self.out.warn(
                "--spa with --kinship: SPA is applied at the score-test layer; the kinship-aware path provides exact variance via the GLMM projection.",
            );
        }
    }
}

/// Walk the manifest's chromosomes and collect rare variants for output
/// writing (sumstats, individual, results headers).
fn load_rare_variants(
    cohort: &CohortHandle<'_>,
    manifest: &CohortManifest,
    maf_cutoff: f64,
) -> Result<Vec<AnnotatedVariant>, CohortError> {
    let mut all = Vec::with_capacity(manifest.n_variants);
    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().unwrap_or(Chromosome::Autosome(1));
        let view = cohort.chromosome(&chrom)?;
        for entry in view.index()?.all_entries() {
            if entry.maf < maf_cutoff {
                all.push(entry.to_annotated_variant(chrom));
            }
        }
    }
    Ok(all)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_config() -> StaarConfig {
        StaarConfig {
            genotypes: PathBuf::from("/tmp/g.vcf.gz"),
            phenotype: PathBuf::from("/tmp/p.tsv"),
            annotations: PathBuf::from("/tmp/a.annotated"),
            trait_names: vec!["BMI".into()],
            covariates: vec!["age".into(), "sex".into()],
            mask_categories: vec![MaskCategory::Coding],
            maf_cutoff: 0.01,
            window_size: 2000,
            individual: false,
            spa: false,
            ancestry_col: None,
            ai_base_tests: 5,
            ai_seed: 7590,
            scang_params: ScangParams {
                lmin: 40,
                lmax: 300,
                step: 10,
            },
            kinship: Vec::new(),
            kinship_groups: None,
            known_loci: None,
            emit_sumstats: false,
            rebuild_store: false,
            column_map: HashMap::new(),
            output_dir: PathBuf::from("/tmp/out"),
            cohort_id: CohortId::new("dummy"),
        }
    }

    #[test]
    fn run_mode_default_is_analyze() {
        let c = dummy_config();
        assert_eq!(c.run_mode(), RunMode::Analyze);
    }

    #[test]
    fn run_mode_emit_sumstats() {
        let mut c = dummy_config();
        c.emit_sumstats = true;
        assert_eq!(c.run_mode(), RunMode::EmitSumstats);
    }

    #[test]
    fn scoring_mode_standard_when_no_special_flags() {
        let c = dummy_config();
        assert_eq!(c.scoring_mode(TraitType::Continuous), ScoringMode::Standard);
        assert_eq!(c.scoring_mode(TraitType::Binary), ScoringMode::Standard);
    }

    #[test]
    fn scoring_mode_spa_only_for_binary() {
        let mut c = dummy_config();
        c.spa = true;
        // Continuous: SPA flag is ignored at the mode-selection layer.
        assert_eq!(c.scoring_mode(TraitType::Continuous), ScoringMode::Standard);
        assert_eq!(c.scoring_mode(TraitType::Binary), ScoringMode::Spa);
    }

    #[test]
    fn scoring_mode_ai_staar_overrides_spa() {
        let mut c = dummy_config();
        c.spa = true;
        c.ancestry_col = Some("super_population".into());
        assert_eq!(c.scoring_mode(TraitType::Binary), ScoringMode::AiStaar);
    }

    #[test]
    fn null_model_kind_dispatch() {
        let mut c = dummy_config();
        assert_eq!(c.null_model_kind(TraitType::Continuous), NullModelKind::Glm);
        assert_eq!(c.null_model_kind(TraitType::Binary), NullModelKind::Logistic);
        c.kinship.push(PathBuf::from("/no/such/k.tsv"));
        assert_eq!(c.null_model_kind(TraitType::Continuous), NullModelKind::KinshipReml);
        assert_eq!(c.null_model_kind(TraitType::Binary), NullModelKind::KinshipPql);
    }

    #[test]
    fn config_hash_round_trips_through_staar_config() {
        // Integration check that StaarConfig wires its fields into
        // ConfigHashInputs the way callers expect: changing any one of
        // them must change the digest. The exhaustive sort/byte-level
        // assertions live in run_manifest::tests::compute_config_hash_*.
        let c1 = dummy_config();
        let h1 = c1.config_hash();
        let mut c2 = dummy_config();
        c2.maf_cutoff = 0.05;
        assert_ne!(h1, c2.config_hash());
        let mut c3 = dummy_config();
        c3.emit_sumstats = true;
        assert_ne!(h1, c3.config_hash());
        let mut c4 = dummy_config();
        c4.covariates = vec!["sex".into(), "age".into()];
        // Covariate order doesn't matter (sorted before hashing).
        assert_eq!(h1, c4.config_hash());
    }

    #[test]
    fn has_kinship_predicate() {
        let mut c = dummy_config();
        assert!(!c.has_kinship());
        c.kinship.push(PathBuf::from("/k.tsv"));
        assert!(c.has_kinship());

        let mut c2 = dummy_config();
        c2.kinship_groups = Some("ancestry".into());
        assert!(c2.has_kinship());
    }

    #[test]
    fn new_run_manifest_uses_config_hash() {
        let c = dummy_config();
        let m = c.new_run_manifest();
        assert_eq!(m.run_mode, c.run_mode());
        assert_eq!(m.trait_name, c.trait_names[0]);
        assert_eq!(m.config_hash, c.config_hash());
    }

    #[test]
    fn manifest_resume_cycle_round_trips_through_tempdir() {
        // Stand-in for the prologue of `StaarPipeline::run` — a prior run
        // marks the durable stages complete, the next process probes the
        // tempdir and the resume planner is asked to skip them. Catches
        // refactors that break the StaarConfig → RunManifest → plan_resume
        // wiring without needing a real VCF.
        let dir = tempfile::tempdir().unwrap();
        let mut config = dummy_config();
        config.output_dir = dir.path().to_path_buf();

        let mut prior = config.new_run_manifest();
        prior.begin(Stage::Validate);
        prior.complete(Stage::Validate);
        prior.begin(Stage::EnsureStore);
        prior.complete(Stage::EnsureStore);
        prior.begin(Stage::EnsureScoreCache);
        prior.complete(Stage::EnsureScoreCache);
        prior.outputs.cache_decisions.push(CacheDecision {
            artifact: ArtifactKind::GenotypeStore,
            outcome: CacheOutcome::Hit,
            reason: "content fingerprint matched".into(),
        });
        prior.write(&config.output_dir).unwrap();

        let probed = RunManifest::probe(&config.output_dir).expect("run.json should round-trip");
        assert_eq!(probed.config_hash, config.config_hash());

        let decision = run_manifest::plan_resume(Some(probed), &config.config_hash());
        let skippable = decision.skippable().to_vec();
        assert!(skippable.contains(&Stage::EnsureStore));
        assert!(skippable.contains(&Stage::EnsureScoreCache));
        assert!(!skippable.contains(&Stage::FitNullModel));
    }

    #[test]
    fn manifest_resume_discards_when_config_changes() {
        // Other half of the contract: editing a hashed field invalidates
        // the prior run.
        let dir = tempfile::tempdir().unwrap();
        let mut config = dummy_config();
        config.output_dir = dir.path().to_path_buf();

        let prior = config.new_run_manifest();
        prior.write(&config.output_dir).unwrap();

        let mut config2 = dummy_config();
        config2.output_dir = dir.path().to_path_buf();
        config2.maf_cutoff = 0.05;

        let probed = RunManifest::probe(&config2.output_dir).unwrap();
        let decision = run_manifest::plan_resume(Some(probed), &config2.config_hash());
        match decision {
            ResumeDecision::Discarded(reason) => {
                assert!(reason.contains("config hash"));
            }
            other => panic!("expected Discarded, got {other:?}"),
        }
    }
}
