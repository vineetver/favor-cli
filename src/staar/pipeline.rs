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
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use faer::Mat;
use serde::{Deserialize, Serialize};

use crate::column::STAAR_WEIGHTS;
use crate::data::{AnnotatedSet, VariantSet, VariantSetKind};
use crate::error::FavorError;
use crate::ingest::ColumnContract;
use crate::output::Output;
use crate::resource::Resources;
use crate::staar::carrier::AnalysisVectors;
use crate::staar::masks::ScangParams;
use crate::staar::model::{augment_covariates, load_known_loci, load_phenotype, NullModel};
use crate::staar::output::{write_individual_results, write_results};
use crate::staar::score_cache;
use crate::staar::scoring::{self, ResultSet};
use crate::staar::store::{self, GenoStoreResult, StoreManifest, STAAR_ANNOTATION_COLUMNS};
use crate::staar::{
    self, ancestry::AncestryInfo, MaskCategory, NullModelKind, RunMode, ScoringMode, TraitType,
};
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
    pub store_dir: PathBuf,
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

    /// Pick the null-model fitter from trait type + kinship presence.
    pub fn null_model_kind(&self, trait_type: TraitType) -> NullModelKind {
        let has_kinship = !self.kinship.is_empty() || self.kinship_groups.is_some();
        match (trait_type, has_kinship) {
            (TraitType::Continuous, false) => NullModelKind::Glm,
            (TraitType::Binary, false) => NullModelKind::Logistic,
            (TraitType::Continuous, true) => NullModelKind::KinshipReml,
            (TraitType::Binary, true) => NullModelKind::KinshipPql,
        }
    }
}

// ---------------------------------------------------------------------------
// Run manifest (artifact lifecycle)
// ---------------------------------------------------------------------------

/// One stage in a `favor staar` run. Order matches execution order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Stage {
    Validate,
    EnsureStore,
    LoadPhenotype,
    FitNullModel,
    EmitSumstats,
    EnsureScoreCache,
    RunScoring,
    WriteResults,
}

impl Stage {
    pub fn label(self) -> &'static str {
        match self {
            Stage::Validate => "validate",
            Stage::EnsureStore => "ensure_store",
            Stage::LoadPhenotype => "load_phenotype",
            Stage::FitNullModel => "fit_null_model",
            Stage::EmitSumstats => "emit_sumstats",
            Stage::EnsureScoreCache => "ensure_score_cache",
            Stage::RunScoring => "run_scoring",
            Stage::WriteResults => "write_results",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StageStatus {
    InProgress,
    Completed,
    Failed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageRecord {
    pub stage: Stage,
    pub status: StageStatus,
    pub started_unix: u64,
    pub completed_unix: Option<u64>,
}

/// Top-level run manifest written into `{output_dir}/run.json`.
///
/// Captures enough state for an operator to ask:
///   - which run mode this is
///   - which stages are complete vs in-progress vs failed
///   - which artifacts were produced
///
/// Resumability hook: probe this file at start of `run()` to log which
/// stages were complete in the last attempt. A future `--resume` can short
/// circuit completed stages off this manifest.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    pub favor_version: String,
    pub run_mode: RunMode,
    pub trait_name: String,
    pub config_hash: String,
    pub started_unix: u64,
    pub stages: Vec<StageRecord>,
    pub outputs: RunOutputs,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct RunOutputs {
    pub store_dir: Option<PathBuf>,
    pub score_cache_dir: Option<PathBuf>,
    pub results_dir: Option<PathBuf>,
    pub sumstats_dir: Option<PathBuf>,
    /// Per-stage cache hit/miss decisions in execution order. Operators can
    /// inspect `run.json` to ask "why did STAAR rebuild the store?" or
    /// "did the score cache hit?" without grepping stdout logs.
    #[serde(default)]
    pub cache_decisions: Vec<CacheDecision>,
}

/// One cache hit/miss decision the pipeline made.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheDecision {
    pub artifact: ArtifactKind,
    pub outcome: CacheOutcome,
    /// Free-text reason — typically the output of `store::describe_miss`
    /// or `"cache_key=<hash> matched"` for hits.
    pub reason: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ArtifactKind {
    GenotypeStore,
    ScoreCache,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CacheOutcome {
    Hit,
    Miss,
    Rebuilt,
}

impl RunManifest {
    fn new(config: &StaarConfig) -> Self {
        Self {
            favor_version: env!("CARGO_PKG_VERSION").to_string(),
            run_mode: config.run_mode(),
            trait_name: config.trait_names[0].clone(),
            config_hash: config_hash(config),
            started_unix: now_unix(),
            stages: Vec::new(),
            outputs: RunOutputs::default(),
        }
    }

    fn begin(&mut self, stage: Stage) {
        self.stages.retain(|r| r.stage != stage);
        self.stages.push(StageRecord {
            stage,
            status: StageStatus::InProgress,
            started_unix: now_unix(),
            completed_unix: None,
        });
    }

    fn complete(&mut self, stage: Stage) {
        if let Some(rec) = self.stages.iter_mut().rev().find(|r| r.stage == stage) {
            rec.status = StageStatus::Completed;
            rec.completed_unix = Some(now_unix());
        }
    }

    fn fail(&mut self, stage: Stage) {
        if let Some(rec) = self.stages.iter_mut().rev().find(|r| r.stage == stage) {
            rec.status = StageStatus::Failed;
            rec.completed_unix = Some(now_unix());
        }
    }

    /// Atomic write to `{output_dir}/run.json` (tmp + fsync + rename + dir
    /// fsync). Reuses `store::write_atomic` so manifest, score cache, and
    /// run.json all share the same durability story.
    pub fn write(&self, output_dir: &Path) -> Result<(), FavorError> {
        let path = output_dir.join("run.json");
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| FavorError::Resource(format!("Serialize run manifest: {e}")))?;
        store::write_atomic(&path, json.as_bytes())
    }

    /// Load an existing manifest if one is present and parses cleanly.
    pub fn probe(output_dir: &Path) -> Option<Self> {
        let path = output_dir.join("run.json");
        let s = std::fs::read_to_string(path).ok()?;
        serde_json::from_str(&s).ok()
    }

    /// Stages that previously completed under this same config.
    pub fn completed_stages(&self) -> Vec<Stage> {
        self.stages
            .iter()
            .filter(|r| r.status == StageStatus::Completed)
            .map(|r| r.stage)
            .collect()
    }
}

fn now_unix() -> u64 {
    SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Stable hash of the inputs that determine which results a run produces.
/// Mirrors `score_cache::cache_key` but covers all fields a stage decision
/// depends on, so a config change invalidates `run.json` cleanly.
fn config_hash(config: &StaarConfig) -> String {
    use sha2::{Digest, Sha256};
    let mut h = Sha256::new();
    h.update(config.genotypes.to_string_lossy().as_bytes());
    h.update(b"|");
    h.update(config.annotations.to_string_lossy().as_bytes());
    h.update(b"|");
    h.update(config.phenotype.to_string_lossy().as_bytes());
    h.update(b"|");
    for t in &config.trait_names {
        h.update(t.as_bytes());
        h.update(b",");
    }
    h.update(b"|");
    let mut covs = config.covariates.clone();
    covs.sort();
    for c in &covs {
        h.update(c.as_bytes());
        h.update(b",");
    }
    h.update(b"|");
    h.update(config.maf_cutoff.to_le_bytes());
    h.update(format!("{:?}", config.run_mode()).as_bytes());
    h.update(b"|");
    for k in &config.kinship {
        h.update(k.to_string_lossy().as_bytes());
        h.update(b",");
    }
    if let Some(g) = &config.kinship_groups {
        h.update(g.as_bytes());
    }
    format!("{:x}", h.finalize())
}

// ---------------------------------------------------------------------------
// Pipeline state passed between stages
// ---------------------------------------------------------------------------

pub type Results = ResultSet;

/// Per-run scoring context: fixed once after the null model is fit.
pub struct ScoringContext {
    pub analysis: AnalysisVectors,
    pub scoring_mode: ScoringMode,
    pub cache_dir: PathBuf,
    pub ancestry: Option<AncestryInfo>,
    /// True only when SPA was requested AND the trait is binary. Threaded
    /// into the AI-STAAR backend so the `--ancestry-col --spa` combination
    /// keeps using SPA inside the per-population kernel.
    pub spa_in_ai_staar: bool,
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

// ---------------------------------------------------------------------------
// Pipeline
// ---------------------------------------------------------------------------

pub struct StaarPipeline<'a> {
    config: StaarConfig,
    out: &'a dyn Output,
    res: Resources,
    manifest: RunManifest,
}

impl<'a> StaarPipeline<'a> {
    pub fn new(config: StaarConfig, out: &'a dyn Output) -> Result<Self, FavorError> {
        let res = store::setup_resources(out)?;
        let manifest = RunManifest::new(&config);
        Ok(Self {
            config,
            out,
            res,
            manifest,
        })
    }

    /// Stage runner. Each stage gates on the previous one returning Ok and
    /// rewrites `run.json` after every transition.
    pub fn run(mut self) -> Result<(), FavorError> {
        std::fs::create_dir_all(&self.config.output_dir).map_err(|e| {
            FavorError::Resource(format!(
                "Cannot create output directory '{}': {e}",
                self.config.output_dir.display()
            ))
        })?;

        if let Some(prior) = RunManifest::probe(&self.config.output_dir) {
            if prior.config_hash == self.manifest.config_hash {
                let done = prior.completed_stages();
                if !done.is_empty() {
                    let labels: Vec<&str> = done.iter().map(|s| s.label()).collect();
                    self.out.status(&format!(
                        "  Found prior run.json (same config). Previously completed: {}",
                        labels.join(", ")
                    ));
                }
            } else {
                self.out
                    .warn("  run.json present from a different config — overwriting");
            }
        }

        let result = self.run_stages();
        if result.is_err() {
            // Keep the manifest as a forensic record on failure.
            let _ = self.manifest.write(&self.config.output_dir);
        }
        result
    }

    fn run_stages(&mut self) -> Result<(), FavorError> {
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
            spa_in_ai_staar: self.config.spa && pheno.trait_type == TraitType::Binary,
        };
        self.warn_mode_combinations(pheno.trait_type, ctx.scoring_mode);

        let (results, individual_pvals) =
            self.stage(Stage::RunScoring, |p| p.stage_run_scoring(&store, &ctx))?;

        let variants = load_rare_variants(&store.store_dir, &store.manifest, self.config.maf_cutoff)?;
        let n_rare = variants.len() as i64;

        self.stage(Stage::WriteResults, |p| {
            p.stage_write_results(
                &results,
                &individual_pvals,
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
    fn stage<T, F>(&mut self, stage: Stage, body: F) -> Result<T, FavorError>
    where
        F: FnOnce(&mut Self) -> Result<T, FavorError>,
    {
        self.manifest.begin(stage);
        // Ignore write errors here — the run is already in a directory we
        // just created. A failure to persist a transient in-progress marker
        // shouldn't kill the run.
        let _ = self.manifest.write(&self.config.output_dir);
        match body(self) {
            Ok(v) => {
                self.manifest.complete(stage);
                let _ = self.manifest.write(&self.config.output_dir);
                Ok(v)
            }
            Err(e) => {
                self.manifest.fail(stage);
                let _ = self.manifest.write(&self.config.output_dir);
                Err(e)
            }
        }
    }

    // -----------------------------------------------------------------------
    // Stages
    // -----------------------------------------------------------------------

    fn stage_validate(&mut self) -> Result<(), FavorError> {
        if !self.config.annotations.exists() {
            return Ok(());
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
            }) => " Your data was annotated with base tier. Re-run: `favor annotate --full`.",
            _ => " Re-run: `favor annotate --full`.",
        };
        Err(FavorError::DataMissing(format!(
            "Missing annotation columns in {}:\n{}\n\
             STAAR requires favor-full annotations.{}",
            self.config.annotations.display(),
            ColumnContract::format_missing(&missing),
            tier_hint,
        )))
    }

    fn stage_ensure_store(&mut self) -> Result<GenoStoreResult, FavorError> {
        // Probe first so we can record the cache decision in run.json
        // regardless of whether the store is hit or rebuilt below.
        let probe_result = if self.config.rebuild_store {
            self.manifest.outputs.cache_decisions.push(CacheDecision {
                artifact: ArtifactKind::GenotypeStore,
                outcome: CacheOutcome::Rebuilt,
                reason: "rebuild requested by --rebuild-store".into(),
            });
            None
        } else {
            let probe = store::probe(
                &self.config.store_dir,
                &self.config.genotypes,
                &self.config.annotations,
            );
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
                    .map(store::describe_miss)
                    .unwrap_or_else(|| "unknown miss".into());
                CacheDecision {
                    artifact: ArtifactKind::GenotypeStore,
                    outcome: CacheOutcome::Miss,
                    reason,
                }
            };
            self.manifest.outputs.cache_decisions.push(decision);
            probe.manifest.is_some().then_some(probe)
        };

        let geno_staging_dir = self.config.output_dir.join(".geno_staging");
        // We re-probe inside `build_or_load_store`; that's a small redundant
        // file read but keeps the build_or_load_store API one-call-and-done
        // for non-pipeline callers (notably the dry-run estimator).
        let _ = probe_result;
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

    fn stage_load_phenotype(
        &mut self,
        store: &GenoStoreResult,
    ) -> Result<PhenoStageOut, FavorError> {
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
    ) -> Result<NullModel, FavorError> {
        self.out.status("Fitting null model...");

        let kind = self.config.null_model_kind(pheno.trait_type);
        match kind {
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
            NullModelKind::KinshipReml | NullModelKind::KinshipPql => {
                let kinships = staar::kinship::load_kinship(
                    &self.config.kinship,
                    &pheno.compact_samples,
                    self.out,
                )?;
                let groups = if let Some(ref col) = self.config.kinship_groups {
                    staar::kinship::load_groups(
                        &store.engine,
                        &self.config.phenotype,
                        col,
                        &pheno.genotype_result,
                        &pheno.pheno_mask,
                        &self.config.column_map,
                        self.out,
                    )?
                } else {
                    staar::kinship::GroupPartition::single(pheno.n)
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
                    _ => unreachable!(),
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
        }
    }

    fn stage_emit_sumstats(
        &mut self,
        store: &GenoStoreResult,
        null_model: &NullModel,
        pheno: &PhenoStageOut,
    ) -> Result<(), FavorError> {
        let sumstats_dir = self.config.sumstats_dir();
        std::fs::create_dir_all(&sumstats_dir).map_err(|e| {
            FavorError::Resource(format!(
                "Cannot create sumstats directory '{}': {e}",
                sumstats_dir.display()
            ))
        })?;
        let variants = load_rare_variants(&store.store_dir, &store.manifest, self.config.maf_cutoff)?;
        let analysis = AnalysisVectors::from_null_model(null_model, &pheno.pheno_mask)?;
        let meta = staar::meta::StudyMeta {
            favor_meta_version: 1,
            trait_type: format!("{:?}", pheno.trait_type),
            trait_name: self.config.trait_names[0].clone(),
            n_samples: pheno.n,
            sigma2: null_model.sigma2,
            maf_cutoff: self.config.maf_cutoff,
            covariates: self.config.covariates.clone(),
            segment_size: 500_000,
        };
        staar::meta::emit_sumstats(
            &store.store_dir,
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
    ) -> Result<PathBuf, FavorError> {
        let key = score_cache::cache_key(
            &store.manifest.key,
            &self.config.trait_names[0],
            &self.config.covariates,
            self.config.known_loci.as_deref(),
            &self.config.kinship,
            self.config.kinship_groups.as_deref(),
        );
        let dir = score_cache::cache_dir(&store.store_dir, &key);

        if score_cache::probe(&store.store_dir, &store.manifest, &key) {
            self.out.status("  Score cache: hit (reusing cached U/K)");
            self.manifest.outputs.cache_decisions.push(CacheDecision {
                artifact: ArtifactKind::ScoreCache,
                outcome: CacheOutcome::Hit,
                reason: format!("cache_key={key} matched"),
            });
            return Ok(dir);
        }

        self.manifest.outputs.cache_decisions.push(CacheDecision {
            artifact: ArtifactKind::ScoreCache,
            outcome: CacheOutcome::Miss,
            reason: format!("no valid scores.bin for cache_key={key}"),
        });

        self.out
            .status("Building score cache (all U/K, no MAF filter)...");
        std::fs::create_dir_all(&dir).map_err(|e| {
            FavorError::Resource(format!("Cannot create score cache directory '{}': {e}", dir.display()))
        })?;

        for ci in &store.manifest.chromosomes {
            let chrom_dir = store.store_dir.join(format!("chromosome={}", ci.name));
            let sg = crate::staar::sparse_g::SparseG::open(&chrom_dir)?;
            let vi = crate::staar::carrier::VariantIndex::load(&chrom_dir)?;
            score_cache::build_chromosome(&sg, &vi, analysis, &dir, &ci.name, self.out)?;
        }

        Ok(dir)
    }

    fn stage_run_scoring(
        &mut self,
        store: &GenoStoreResult,
        ctx: &ScoringContext,
    ) -> Result<(Results, Vec<(usize, f64)>), FavorError> {
        scoring::run_score_tests(
            &store.store_dir,
            &store.manifest,
            &self.config,
            &ctx.analysis,
            ctx,
            self.out,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn stage_write_results(
        &mut self,
        results: &Results,
        individual_pvals: &[(usize, f64)],
        variants: &[AnnotatedVariant],
        null_model: &NullModel,
        trait_type: TraitType,
        n: usize,
        n_rare: i64,
    ) -> Result<(), FavorError> {
        let results_dir = self.config.results_dir();
        std::fs::create_dir_all(&results_dir).map_err(|e| {
            FavorError::Resource(format!(
                "Cannot create results directory '{}': {e}",
                results_dir.display()
            ))
        })?;

        if self.config.individual && !individual_pvals.is_empty() {
            write_individual_results(individual_pvals, variants, &results_dir, self.out)?;
        }

        write_results(
            results,
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
        let has_kinship = !self.config.kinship.is_empty() || self.config.kinship_groups.is_some();
        if has_kinship && mode == ScoringMode::Spa {
            self.out.warn(
                "--spa with --kinship: SPA is applied at the score-test layer; the kinship-aware path provides exact variance via the GLMM projection.",
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Walk the manifest's chromosomes and collect rare variants for output
/// writing (sumstats, individual, results headers).
fn load_rare_variants(
    store_dir: &Path,
    manifest: &StoreManifest,
    maf_cutoff: f64,
) -> Result<Vec<AnnotatedVariant>, FavorError> {
    let mut all = Vec::with_capacity(manifest.n_variants);
    for ci in &manifest.chromosomes {
        let chrom: Chromosome = ci.name.parse().unwrap_or(Chromosome::Autosome(1));
        let chrom_dir = store_dir.join(format!("chromosome={}", ci.name));
        let index = crate::staar::carrier::VariantIndex::load(&chrom_dir)?;
        for entry in index.all_entries() {
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
            store_dir: PathBuf::from("/tmp/out/store"),
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
    fn cache_decisions_round_trip_through_run_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let mut m = RunManifest::new(&dummy_config());
        m.outputs.cache_decisions.push(CacheDecision {
            artifact: ArtifactKind::GenotypeStore,
            outcome: CacheOutcome::Miss,
            reason: "no manifest.json on disk".into(),
        });
        m.outputs.cache_decisions.push(CacheDecision {
            artifact: ArtifactKind::ScoreCache,
            outcome: CacheOutcome::Hit,
            reason: "cache_key=abc matched".into(),
        });
        m.write(dir.path()).unwrap();

        let loaded = RunManifest::probe(dir.path()).unwrap();
        assert_eq!(loaded.outputs.cache_decisions.len(), 2);
        assert_eq!(loaded.outputs.cache_decisions[0].outcome, CacheOutcome::Miss);
        assert_eq!(loaded.outputs.cache_decisions[0].artifact, ArtifactKind::GenotypeStore);
        assert_eq!(loaded.outputs.cache_decisions[1].outcome, CacheOutcome::Hit);
        assert_eq!(loaded.outputs.cache_decisions[1].artifact, ArtifactKind::ScoreCache);
    }

    #[test]
    fn run_manifest_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let c = dummy_config();
        let mut m = RunManifest::new(&c);
        m.begin(Stage::Validate);
        m.complete(Stage::Validate);
        m.begin(Stage::EnsureStore);
        m.complete(Stage::EnsureStore);
        m.outputs.store_dir = Some(PathBuf::from("/tmp/store"));
        m.write(dir.path()).unwrap();

        let loaded = RunManifest::probe(dir.path()).unwrap();
        assert_eq!(loaded.config_hash, m.config_hash);
        assert_eq!(loaded.run_mode, RunMode::Analyze);
        assert_eq!(loaded.completed_stages(), vec![Stage::Validate, Stage::EnsureStore]);
        assert_eq!(loaded.outputs.store_dir, Some(PathBuf::from("/tmp/store")));
    }

    #[test]
    fn run_manifest_fail_marks_stage() {
        let mut m = RunManifest::new(&dummy_config());
        m.begin(Stage::FitNullModel);
        m.fail(Stage::FitNullModel);
        let rec = m.stages.last().unwrap();
        assert_eq!(rec.status, StageStatus::Failed);
        assert!(rec.completed_unix.is_some());
    }

    #[test]
    fn config_hash_changes_with_inputs() {
        let c1 = dummy_config();
        let h1 = config_hash(&c1);
        let mut c2 = dummy_config();
        c2.maf_cutoff = 0.05;
        assert_ne!(h1, config_hash(&c2));
        let mut c3 = dummy_config();
        c3.emit_sumstats = true;
        assert_ne!(h1, config_hash(&c3));
        let mut c4 = dummy_config();
        c4.covariates = vec!["sex".into(), "age".into()];
        // Covariate order doesn't matter (sorted before hashing).
        assert_eq!(h1, config_hash(&c4));
    }
}
