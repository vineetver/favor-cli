//! STAAR analysis pipeline orchestration.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use faer::Mat;

use crate::column::STAAR_PHRED_CHANNELS;
use crate::error::CohortError;
use crate::ingest::ColumnContract;
use crate::output::{bail_if_cancelled, Output};
use crate::runtime::Engine;
use crate::staar::carrier::AnalysisVectors;
use crate::staar::masks::ScangParams;
use crate::staar::model::{augment_covariates, load_known_loci, load_phenotype, NullModel};
use crate::staar::multi::MultiNullContinuous;
use crate::staar::output::{write_individual_results, write_results, NullMeta};
use crate::staar::run_manifest::{
    self, ArtifactKind, CacheDecision, CacheOutcome, ConfigHashInputs, ResumeDecision,
    ResumeSummary, RunManifest, Stage,
};
use crate::staar::scoring::{self, MultiScoringRequest, ResultSet, ScoringRequest};
use crate::staar::{
    self, ancestry::AncestryInfo, MaskCategory, NullModelKind, RunMode, ScoringMode, TraitType,
};
use crate::store::cache::{null_model_cache, score_cache};
use crate::store::cohort::{
    self, CohortHandle, CohortId, CohortManifest, GenoStoreResult, STAAR_ANNOTATION_COLUMNS,
};
use crate::store::list::VariantSet;
use crate::types::{AnnotatedVariant, Chromosome};

pub enum CohortSource {
    Existing,
    Fresh {
        genotypes: Vec<PathBuf>,
        annotations: PathBuf,
    },
}

pub struct StaarConfig {
    pub cohort_source: CohortSource,
    pub phenotype: PathBuf,
    pub trait_names: Vec<String>,
    pub covariates: Vec<String>,
    pub mask_categories: Vec<MaskCategory>,
    pub maf_cutoff: f64,
    pub window_size: u32,
    pub individual: bool,
    pub spa: bool,
    pub ancestry_col: Option<String>,
    pub ai_base_tests: usize,
    pub ai_seed: u64,
    pub scang_params: ScangParams,
    pub kinship: Vec<PathBuf>,
    pub kinship_groups: Option<String>,
    pub known_loci: Option<PathBuf>,
    pub null_model_path: Option<PathBuf>,
    pub run_mode: RunMode,
    pub rebuild_store: bool,
    pub column_map: HashMap<String, String>,
    pub output_dir: PathBuf,
    pub cohort_id: CohortId,
}

impl StaarConfig {
    pub fn genotypes(&self) -> Option<&[PathBuf]> {
        match &self.cohort_source {
            CohortSource::Fresh { genotypes, .. } => Some(genotypes.as_slice()),
            CohortSource::Existing => None,
        }
    }

    pub fn annotations(&self) -> Option<&Path> {
        match &self.cohort_source {
            CohortSource::Fresh { annotations, .. } => Some(annotations.as_path()),
            CohortSource::Existing => None,
        }
    }
}

impl StaarConfig {
    /// Directory where the per-mask parquet results and `staar.meta.json`
    /// land. Single-trait runs shelve results under the trait name so
    /// parallel trait runs in the same output dir do not stomp each
    /// other; multi-trait runs use the literal label `"multi"` because
    /// the trait list is already embedded in the metadata and any
    /// `t1+t2+...` join scheme produces filesystem-illegal or
    /// unreadably long paths.
    pub fn results_dir(&self) -> PathBuf {
        let label = if self.trait_names.len() > 1 {
            "multi".to_string()
        } else {
            self.trait_names[0].clone()
        };
        self.output_dir.join("results").join(label)
    }

    pub fn sumstats_dir(&self) -> PathBuf {
        self.output_dir.join("sumstats").join(&self.trait_names[0])
    }

    pub fn scoring_mode(&self, trait_type: TraitType) -> ScoringMode {
        if self.ancestry_col.is_some() {
            ScoringMode::AiStaar
        } else if self.spa && trait_type == TraitType::Binary {
            ScoringMode::Spa
        } else {
            ScoringMode::Standard
        }
    }

    pub fn has_kinship(&self) -> bool {
        !self.kinship.is_empty() || self.kinship_groups.is_some()
    }

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
    pub fn hash_inputs<'a>(
        &'a self,
        existing_cohort_manifest: Option<&'a Path>,
    ) -> ConfigHashInputs<'a> {
        let (genotypes, annotations): (Vec<PathBuf>, &Path) = match &self.cohort_source {
            CohortSource::Fresh {
                genotypes,
                annotations,
            } => (genotypes.clone(), annotations),
            CohortSource::Existing => {
                let p = existing_cohort_manifest
                    .expect("existing cohort config hash requires a resolved manifest path");
                (vec![p.to_path_buf()], p)
            }
        };
        ConfigHashInputs {
            genotypes,
            annotations,
            phenotype: &self.phenotype,
            trait_names: &self.trait_names,
            covariates: &self.covariates,
            maf_cutoff: self.maf_cutoff,
            run_mode: self.run_mode,
            kinship: &self.kinship,
            kinship_groups: self.kinship_groups.as_deref(),
        }
    }

    pub fn config_hash(&self, existing_cohort_manifest: Option<&Path>) -> String {
        run_manifest::compute_config_hash(&self.hash_inputs(existing_cohort_manifest))
    }

    pub fn new_run_manifest(&self, existing_cohort_manifest: Option<&Path>) -> RunManifest {
        // Joint multi-trait runs encode every trait in the manifest label
        // so `run.json` on disk records the exact trait list, while the
        // results directory gets the short `"multi"` alias.
        let label = if self.trait_names.len() > 1 {
            format!("multi:{}", self.trait_names.join("+"))
        } else {
            self.trait_names[0].clone()
        };
        RunManifest::new(
            self.run_mode,
            label,
            self.config_hash(existing_cohort_manifest),
        )
    }
}

pub struct ScoringOutput {
    pub results: ResultSet,
    pub individual: Vec<crate::staar::output::IndividualRow>,
}

struct PhenoStageOut {
    y: Mat<f64>,
    x: Mat<f64>,
    trait_type: TraitType,
    n: usize,
    pheno_mask: Vec<bool>,
    ancestry: Option<AncestryInfo>,
    compact_samples: Vec<String>,
}

pub struct StaarPipeline<'a> {
    config: StaarConfig,
    engine: &'a Engine,
    out: &'a dyn Output,
    manifest: RunManifest,
    resume: ResumeDecision,
}

impl<'a> StaarPipeline<'a> {
    pub fn new(
        config: StaarConfig,
        engine: &'a Engine,
        out: &'a dyn Output,
    ) -> Result<Self, CohortError> {
        let res = engine.resources();
        out.status(&format!(
            "STAAR: {} memory, {} threads ({})",
            res.memory_human(),
            res.threads,
            res.environment()
        ));
        let existing_manifest = existing_cohort_manifest_path(engine, &config);
        let manifest = config.new_run_manifest(existing_manifest.as_deref());
        Ok(Self {
            config,
            engine,
            out,
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

        // Validate / EnsureStore / LoadPhenotype are shared between the
        // single-trait and joint multi-trait paths. After phenotype load
        // the two paths diverge: single-trait threads an `AnalysisVectors`
        // through the score cache and the cached-U/K kernels; multi-trait
        // fits a joint null and walks genes directly against the sparse
        // genotype store with no cache. Branching at the type level here
        // keeps the single-trait stages below from ever seeing a
        // `(n × k)` Y.
        match self.config.run_mode {
            RunMode::Analyze | RunMode::EmitSumstats => self.run_single_trait(&store, &pheno),
            RunMode::MultiTrait => self.run_multi_trait(&store, &pheno),
        }
    }

    fn run_single_trait(
        &mut self,
        store: &GenoStoreResult,
        pheno: &PhenoStageOut,
    ) -> Result<(), CohortError> {
        let null_model = self.stage(Stage::FitNullModel, |p| {
            p.stage_fit_null_model(pheno, store)
        })?;

        // Sumstats early-exit run mode.
        if self.config.run_mode == RunMode::EmitSumstats {
            self.stage(Stage::EmitSumstats, |p| {
                p.stage_emit_sumstats(store, &null_model, pheno)
            })?;
            self.manifest.write(&self.config.output_dir)?;
            return Ok(());
        }

        let analysis = AnalysisVectors::from_null_model(&null_model, &pheno.pheno_mask)?;

        let cache_dir = self.stage(Stage::EnsureScoreCache, |p| {
            p.stage_ensure_score_cache(store, &analysis)
        })?;
        self.manifest.outputs.score_cache_dir = Some(cache_dir.clone());
        self.manifest.write(&self.config.output_dir)?;

        let scoring_mode = self.config.scoring_mode(pheno.trait_type);
        self.warn_mode_combinations(pheno.trait_type, scoring_mode);

        let scoring = self.stage(Stage::RunScoring, |p| {
            p.stage_run_scoring(
                store,
                &analysis,
                scoring_mode,
                &cache_dir,
                pheno.ancestry.as_ref(),
                pheno.trait_type,
            )
        })?;

        // load_rare_variants walks every chromosome and reads each VariantIndex
        // through the cohort handle. Hold the result in a local so the
        // WriteResults stage doesn't trigger a second walk.
        let cohort = self.cohort();
        let n_rare = load_rare_variants(&cohort, &store.manifest, self.config.maf_cutoff)?
            .len() as i64;

        self.stage(Stage::WriteResults, |p| {
            p.stage_write_results(
                &scoring,
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

    fn run_multi_trait(
        &mut self,
        store: &GenoStoreResult,
        pheno: &PhenoStageOut,
    ) -> Result<(), CohortError> {
        // Multi-trait joint STAAR is unrelated continuous only; the
        // binary reject fires inside `load_phenotype` via the trait-type
        // detection, so by the time we reach this method `pheno.y` is an
        // `(n, k)` continuous matrix with `k >= 2`. Fit the joint null
        // directly from the multi module.
        let null = self.stage(Stage::FitNullModel, |p| p.stage_fit_multi_null(pheno))?;

        let scoring = self.stage(Stage::RunScoring, |p| {
            p.stage_run_multi_scoring(store, &null, pheno)
        })?;

        let cohort = self.cohort();
        let variants = load_rare_variants(&cohort, &store.manifest, self.config.maf_cutoff)?;
        let n_rare = variants.len() as i64;

        self.stage(Stage::WriteResults, |p| {
            p.stage_write_multi_results(&scoring, &variants, &null, pheno.n, n_rare)
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
        // For the `Existing` cohort path the operator is loading a cohort
        // by id; the original annotated set may have been deleted. The
        // tier check happened at `favoringest` time, baked into the
        // cohort store via the rebuild fingerprint. Just probe the
        // manifest exists; the load stage gives a typed error otherwise.
        let annotations = match &self.config.cohort_source {
            CohortSource::Fresh { annotations, .. } => annotations,
            CohortSource::Existing => {
                let manifest_path = self.cohort_manifest_path();
                if !manifest_path.exists() {
                    return Err(CohortError::DataMissing(format!(
                        "Cohort '{}' not found at {}. Run `favoringest <vcf> \
                         --annotations <set> --cohort-id {}` first.",
                        self.config.cohort_id.as_str(),
                        manifest_path
                            .parent()
                            .map(|p| p.display().to_string())
                            .unwrap_or_default(),
                        self.config.cohort_id.as_str(),
                    )));
                }
                return Ok(());
            }
        };

        // `build_config` already verified the annotation path exists; if it
        // disappeared between then and now we want a hard error, not silent
        // skip — surface the missing-file message rather than no-op.
        if !annotations.exists() {
            return Err(CohortError::Input(format!(
                "Annotations no longer at '{}'. Re-run `favorannotate` or check the path.",
                annotations.display()
            )));
        }

        // Typed tier check: STAAR needs all 11 annotation weight columns.
        // `supports` errors if the set is not annotated, or if any weight
        // column is missing for the tier it was annotated against. The
        // `require_staar_weight_catalog` call below additionally confirms
        // the parquet actually has each column by name — catches partial
        // writes and out-of-band regeneration the tier metadata misses.
        let ann_vs = VariantSet::open(annotations)?;
        let weight_cols: Vec<crate::column::Col> = STAAR_PHRED_CHANNELS.to_vec();
        ann_vs.supports(&weight_cols)?;
        ann_vs.require_staar_weight_catalog()?;
        ann_vs.require_structural_annotation_catalog()?;

        // Raw annotation column contract.
        let contract = ColumnContract {
            command: "staar",
            required: STAAR_ANNOTATION_COLUMNS,
        };
        let missing = contract.check(ann_vs.columns());
        if missing.is_empty() {
            return Ok(());
        }
        let tier_hint =
            " Re-run: `favor annotate` to produce a complete annotated set.";
        Err(CohortError::DataMissing(format!(
            "Missing annotation columns in {}:\n{}\n\
             STAAR needs all 11 annotation weight channels.{}",
            annotations.display(),
            ColumnContract::format_missing(&missing),
            tier_hint,
        )))
    }

    fn stage_ensure_store(&mut self) -> Result<GenoStoreResult, CohortError> {
        // Existing cohort path: trust the manifest, no rebuild, no probe.
        // The operator already built this cohort via `favoringest`.
        if matches!(self.config.cohort_source, CohortSource::Existing) {
            self.manifest.outputs.cache_decisions.push(CacheDecision {
                artifact: ArtifactKind::GenotypeStore,
                outcome: CacheOutcome::Hit,
                reason: "loaded by cohort id".into(),
            });
            return self.cohort().load();
        }

        let (genotypes, annotations) = match &self.config.cohort_source {
            CohortSource::Fresh {
                genotypes,
                annotations,
            } => (genotypes.clone(), annotations.clone()),
            CohortSource::Existing => unreachable!("handled above"),
        };

        // One probe — the recorded cache decision and the path actually
        // taken in `build_or_load` are sourced from the same StoreProbe
        // so they cannot drift.
        let probe = if self.config.rebuild_store {
            self.manifest.outputs.cache_decisions.push(CacheDecision {
                artifact: ArtifactKind::GenotypeStore,
                outcome: CacheOutcome::Rebuilt,
                reason: "rebuild requested by --rebuild-store".into(),
            });
            // Drop any score-cache + lookup index that referenced this
            // cohort under its current id; the rebuild changes the
            // content fingerprint and would orphan them otherwise.
            if let Err(e) = self
                .engine
                .store()
                .cache()
                .prune_cohort(&self.config.cohort_id)
            {
                self.out.warn(&format!("  prune cache before rebuild: {e}"));
            }
            cohort::StoreProbe {
                store_dir: self.cohort().dir().to_path_buf(),
                manifest: None,
                miss_reason: None,
            }
        } else {
            let probe = self.cohort().probe(&genotypes, &annotations);
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
        self.cohort().build_or_load(
            crate::store::cohort::CohortSources {
                genotypes: &genotypes,
                annotations: &annotations,
            },
            crate::store::cohort::BuildOpts {
                staging_dir: &geno_staging_dir,
                rebuild: self.config.rebuild_store,
                probe,
            },
            self.engine,
            self.out,
        )
    }

    fn stage_load_phenotype(
        &mut self,
        store: &GenoStoreResult,
    ) -> Result<PhenoStageOut, CohortError> {
        let pheno = load_phenotype(
            self.engine.df(),
            &self.config.phenotype,
            &self.config.covariates,
            &store.geno,
            &self.config.trait_names,
            self.config.ancestry_col.as_deref(),
            self.config.ai_base_tests,
            self.config.ai_seed,
            &self.config.column_map,
            self.out,
        )?;
        let mut x = pheno.x;

        if let Some(ref loci_path) = self.config.known_loci {
            let x_cond = load_known_loci(self.engine.df(), &store.geno, loci_path, pheno.n, self.out)?;
            x = augment_covariates(&x, &x_cond);
            self.out.status(&format!(
                "  Conditional: {} known loci added as covariates",
                x_cond.ncols()
            ));
        }

        let compact_samples: Vec<String> = store
            .geno
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
        })
    }

    fn stage_fit_null_model(
        &mut self,
        pheno: &PhenoStageOut,
        store: &GenoStoreResult,
    ) -> Result<NullModel, CohortError> {
        if let Some(path) = self.config.null_model_path.as_ref() {
            if self.config.has_kinship() {
                return Err(CohortError::Input(
                    "--null-model is not yet compatible with --kinship; drop --kinship to \
                     use the imported null, or remove --null-model to refit with kinship."
                        .into(),
                ));
            }
            self.out.status(&format!(
                "Loading null model from {}...",
                path.display(),
            ));
            let nm = null_model_cache::load_from_file(path)?;
            if nm.n_samples != pheno.n {
                return Err(CohortError::Input(format!(
                    "imported null model has n_samples={}, phenotype has n={}",
                    nm.n_samples, pheno.n,
                )));
            }
            self.out.status(&format!("  sigma2 = {:.4}", nm.sigma2));
            return Ok(nm);
        }

        let has_kinship = self.config.has_kinship();

        // Cache lookup (skip for kinship models — KinshipState contains
        // sparse Cholesky factors that are not serializable in v1).
        let cache_dir = if !has_kinship {
            let key_str = null_model_cache::cache_key(
                &store.manifest.key,
                &self.config.trait_names[0],
                &self.config.covariates,
                self.config.known_loci.as_deref(),
                &self.config.kinship,
                self.config.kinship_groups.as_deref(),
            );
            let dir = self
                .engine
                .store()
                .cache()
                .null_model_cache_dir(
                    &self.config.cohort_id,
                    &crate::store::ids::CacheKey::new(&key_str),
                );
            if null_model_cache::probe(&dir) {
                self.out.status("Null model cache: hit");
                self.manifest.outputs.cache_decisions.push(CacheDecision {
                    artifact: ArtifactKind::GenotypeStore,
                    outcome: CacheOutcome::Hit,
                    reason: "null model cache hit".into(),
                });
                let nm = null_model_cache::load(&dir)?;
                self.out.status(&format!("  sigma2 = {:.4}", nm.sigma2));
                return Ok(nm);
            }
            Some(dir)
        } else {
            self.out
                .status("  Null model cache: skipped (kinship models not cached in v1)");
            None
        };

        self.out.status("Fitting null model...");

        let nm = match self.config.null_model_kind(pheno.trait_type) {
            NullModelKind::Glm => staar::model::fit_glm(&pheno.y, &pheno.x),
            NullModelKind::Logistic => staar::model::fit_logistic(&pheno.y, &pheno.x, 25),
            kind @ (NullModelKind::KinshipReml | NullModelKind::KinshipPql) => {
                self.fit_kinship_null_model(kind, pheno, store)?
            }
        };
        self.out.status(&format!("  sigma2 = {:.4}", nm.sigma2));

        if let Some(ref dir) = cache_dir {
            if let Err(e) = null_model_cache::save(dir, &nm) {
                self.out
                    .warn(&format!("  Null model cache: save failed ({e})"));
            } else {
                self.out.status("  Null model: fitted and cached");
            }
        }

        Ok(nm)
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
        let kinships =
            staar::kinship::load_kinship(&self.config.kinship, &pheno.compact_samples, self.out)?;
        let groups = match self.config.kinship_groups.as_deref() {
            Some(col) => staar::kinship::load_groups(
                self.engine.df(),
                &self.config.phenotype,
                col,
                &store.geno,
                &pheno.pheno_mask,
                &self.config.column_map,
                self.out,
            )?,
            None => staar::kinship::GroupPartition::single(pheno.n),
        };

        let budget = self.engine.resources().kinship_budget_bytes;
        let state = match kind {
            NullModelKind::KinshipReml => {
                staar::kinship::fit_reml(&pheno.y, &pheno.x, &kinships, &groups, None, budget)?
            }
            NullModelKind::KinshipPql => staar::kinship::fit_pql_glmm(
                &pheno.y, &pheno.x, &kinships, &groups, self.out, budget,
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

    /// Fit the joint multi-trait null model. Unrelated continuous only;
    /// the binary and kinship combinations are rejected at config build
    /// time, and the binary-first-trait reject fires inside
    /// `load_phenotype` before this stage ever runs.
    fn stage_fit_multi_null(
        &mut self,
        pheno: &PhenoStageOut,
    ) -> Result<MultiNullContinuous, CohortError> {
        self.out.status("Fitting joint multi-trait null model...");
        let null = MultiNullContinuous::fit(&pheno.y, &pheno.x);
        self.out.status(&format!(
            "  k = {} traits, n = {} samples, Σ_res fitted",
            null.n_pheno, null.n_samples,
        ));
        Ok(null)
    }

    fn stage_run_multi_scoring(
        &mut self,
        store: &GenoStoreResult,
        null: &MultiNullContinuous,
        pheno: &PhenoStageOut,
    ) -> Result<ScoringOutput, CohortError> {
        let request = MultiScoringRequest {
            mask_categories: &self.config.mask_categories,
            maf_cutoff: self.config.maf_cutoff,
            window_size: self.config.window_size,
            scang_params: &self.config.scang_params,
        };
        let results = scoring::run_multi_score_tests(
            &self.cohort(),
            &store.manifest,
            &request,
            null,
            &pheno.pheno_mask,
            self.out,
        )?;
        Ok(ScoringOutput {
            results,
            individual: Vec::new(),
        })
    }

    fn stage_write_multi_results(
        &mut self,
        scoring: &ScoringOutput,
        _variants: &[AnnotatedVariant],
        null: &MultiNullContinuous,
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

        // Multi-trait path does not produce `--individual` per-variant
        // p-values; the joint kernel works on dense G₀ rather than the
        // per-variant score cache the single-trait individual-test
        // producer feeds off.
        write_results(
            &scoring.results,
            &self.config.trait_names,
            self.config.maf_cutoff,
            &results_dir,
            NullMeta::Multi(null),
            TraitType::Continuous,
            n,
            n_rare,
            self.config.known_loci.is_some(),
            self.out,
        )
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
        let cohort = self.cohort();
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
            self.engine,
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
            .engine
            .store()
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
            CohortError::Resource(format!(
                "Cannot create score cache directory '{}': {e}",
                dir.display()
            ))
        })?;

        let cohort = self.cohort();
        for ci in &store.manifest.chromosomes {
            let chrom: Chromosome = ci.name.parse().map_err(|e: String| CohortError::Input(e))?;
            let view = cohort.chromosome(&chrom)?;
            score_cache::build_chromosome(
                &view,
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
        analysis: &AnalysisVectors,
        scoring_mode: ScoringMode,
        cache_dir: &Path,
        ancestry: Option<&AncestryInfo>,
        trait_type: TraitType,
    ) -> Result<ScoringOutput, CohortError> {
        let request = ScoringRequest {
            mask_categories: &self.config.mask_categories,
            maf_cutoff: self.config.maf_cutoff,
            window_size: self.config.window_size,
            scang_params: &self.config.scang_params,
            mode: scoring_mode,
            cache_dir,
            ancestry,
            spa_active: self.config.spa && trait_type == TraitType::Binary,
            // Single-trait gaussian no-kinship path only in this build.
            // Binary / kinship / multi land in follow-ups.
            individual: self.config.individual
                && trait_type == TraitType::Continuous
                && !self.config.has_kinship(),
            individual_mac_cutoff: 20,
        };
        let (results, individual) = scoring::run_score_tests(
            &self.cohort(),
            &store.manifest,
            &request,
            analysis,
            self.out,
        )?;
        Ok(ScoringOutput {
            results,
            individual,
        })
    }

    fn stage_write_results(
        &mut self,
        scoring: &ScoringOutput,
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

        let is_conditional = self.config.known_loci.is_some();

        if self.config.individual && !scoring.individual.is_empty() {
            write_individual_results(
                &scoring.individual,
                &results_dir,
                is_conditional,
                self.out,
            )?;
        }

        write_results(
            &scoring.results,
            &self.config.trait_names,
            self.config.maf_cutoff,
            &results_dir,
            NullMeta::Single {
                sigma2: null_model.sigma2,
            },
            trait_type,
            n,
            n_rare,
            is_conditional,
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
        if self.config.known_loci.is_some() && mode == ScoringMode::Spa {
            self.out.warn(
                "--spa with --known-loci: the known-loci adjustment is absorbed into the fitted null, but the SPA saddlepoint reads raw G; gene-mask p-values from this run are SPA-adjusted but not conditional. A separate --cond-spa path is tracked in STAARpipeline as Gene_Centric_*_cond_spa.",
            );
        }
    }

    fn cohort(&self) -> CohortHandle<'_> {
        self.engine.cohort(&self.config.cohort_id)
    }

    fn cohort_manifest_path(&self) -> PathBuf {
        self.cohort().dir().join("manifest.json")
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

fn existing_cohort_manifest_path(engine: &Engine, config: &StaarConfig) -> Option<PathBuf> {
    match config.cohort_source {
        CohortSource::Existing => {
            Some(engine.cohort(&config.cohort_id).dir().join("manifest.json"))
        }
        CohortSource::Fresh { .. } => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_config() -> StaarConfig {
        StaarConfig {
            cohort_source: CohortSource::Fresh {
                genotypes: vec![PathBuf::from("/tmp/g.vcf.gz")],
                annotations: PathBuf::from("/tmp/a.annotated"),
            },
            phenotype: PathBuf::from("/tmp/p.tsv"),
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
            null_model_path: None,
            run_mode: RunMode::Analyze,
            rebuild_store: false,
            column_map: HashMap::new(),
            output_dir: PathBuf::from("/tmp/out"),
            cohort_id: CohortId::new("dummy"),
        }
    }

    #[test]
    fn run_mode_default_is_analyze() {
        let c = dummy_config();
        assert_eq!(c.run_mode, RunMode::Analyze);
    }

    #[test]
    fn run_mode_emit_sumstats() {
        let mut c = dummy_config();
        c.run_mode = RunMode::EmitSumstats;
        assert_eq!(c.run_mode, RunMode::EmitSumstats);
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
        assert_eq!(
            c.null_model_kind(TraitType::Binary),
            NullModelKind::Logistic
        );
        c.kinship.push(PathBuf::from("/no/such/k.tsv"));
        assert_eq!(
            c.null_model_kind(TraitType::Continuous),
            NullModelKind::KinshipReml
        );
        assert_eq!(
            c.null_model_kind(TraitType::Binary),
            NullModelKind::KinshipPql
        );
    }

    #[test]
    fn config_hash_round_trips_through_staar_config() {
        // Integration check that StaarConfig wires its fields into
        // ConfigHashInputs the way callers expect: changing any one of
        // them must change the digest. The exhaustive sort/byte-level
        // assertions live in run_manifest::tests::compute_config_hash_*.
        let c1 = dummy_config();
        let h1 = c1.config_hash(None);
        let mut c2 = dummy_config();
        c2.maf_cutoff = 0.05;
        assert_ne!(h1, c2.config_hash(None));
        let mut c3 = dummy_config();
        c3.run_mode = RunMode::EmitSumstats;
        assert_ne!(h1, c3.config_hash(None));
        let mut c4 = dummy_config();
        c4.covariates = vec!["sex".into(), "age".into()];
        // Covariate order doesn't matter (sorted before hashing).
        assert_eq!(h1, c4.config_hash(None));
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
        let m = c.new_run_manifest(None);
        assert_eq!(m.run_mode, c.run_mode);
        assert_eq!(m.trait_name, c.trait_names[0]);
        assert_eq!(m.config_hash, c.config_hash(None));
    }

    /// Single-trait runs keep the existing per-trait results directory
    /// layout so two parallel runs in the same output dir do not stomp
    /// each other.
    #[test]
    fn results_dir_single_trait_is_per_trait() {
        let c = dummy_config();
        let dir = c.results_dir();
        assert!(
            dir.ends_with("results/BMI"),
            "expected per-trait directory, got {}",
            dir.display()
        );
    }

    /// Multi-trait runs shelve results under the stable `"multi"` label;
    /// the exact trait list lives in `staar.meta.json` so the directory
    /// name stays short and filesystem-legal.
    #[test]
    fn results_dir_multi_trait_is_multi_label() {
        let mut c = dummy_config();
        c.trait_names = vec!["BMI".into(), "HEIGHT".into(), "LDL".into()];
        c.run_mode = RunMode::MultiTrait;
        let dir = c.results_dir();
        assert!(
            dir.ends_with("results/multi"),
            "expected results/multi, got {}",
            dir.display()
        );
    }

    /// The manifest label for a multi run embeds every trait so
    /// forensic `run.json` reads still identify the exact inputs.
    #[test]
    fn new_run_manifest_label_for_multi_trait() {
        let mut c = dummy_config();
        c.trait_names = vec!["BMI".into(), "HEIGHT".into()];
        c.run_mode = RunMode::MultiTrait;
        let m = c.new_run_manifest(None);
        assert_eq!(m.run_mode, RunMode::MultiTrait);
        assert_eq!(m.trait_name, "multi:BMI+HEIGHT");
    }

    /// Config-hash invariance: adding a second trait must change the
    /// hash so the resume planner does not reuse a single-trait run.
    #[test]
    fn config_hash_changes_when_second_trait_added() {
        let c1 = dummy_config();
        let mut c2 = dummy_config();
        c2.trait_names.push("HEIGHT".into());
        c2.run_mode = RunMode::MultiTrait;
        assert_ne!(c1.config_hash(None), c2.config_hash(None));
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

        let mut prior = config.new_run_manifest(None);
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
        assert_eq!(probed.config_hash, config.config_hash(None));

        let decision = run_manifest::plan_resume(Some(probed), &config.config_hash(None));
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

        let prior = config.new_run_manifest(None);
        prior.write(&config.output_dir).unwrap();

        let mut config2 = dummy_config();
        config2.output_dir = dir.path().to_path_buf();
        config2.maf_cutoff = 0.05;

        let probed = RunManifest::probe(&config2.output_dir).unwrap();
        let decision = run_manifest::plan_resume(Some(probed), &config2.config_hash(None));
        match decision {
            ResumeDecision::Discarded(reason) => {
                assert!(reason.contains("config hash"));
            }
            other => panic!("expected Discarded, got {other:?}"),
        }
    }
}
