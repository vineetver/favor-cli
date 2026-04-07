//! `run.json` shape and resume planner for `cohort staar`.
//!
//! Owns the `Stage` enum, the on-disk manifest, the cache-decision log,
//! and `plan_resume`. Disk writes go through `store::write_atomic` so a
//! crash leaves either no `run.json` or the prior version.

use std::path::{Path, PathBuf};
use std::time::SystemTime;

use serde::{Deserialize, Serialize};

use crate::error::CohortError;
use crate::staar::store;
use crate::staar::RunMode;

/// One stage in a `cohort staar` run. Order matches execution order.
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
    /// All stages in execution order. The `stage_labels_match_serde` and
    /// `stage_all_matches_run_order` tests use this to catch a new stage
    /// added without updating both the serde rename and the human label.
    #[allow(dead_code)]
    pub const ALL: [Stage; 8] = [
        Stage::Validate,
        Stage::EnsureStore,
        Stage::LoadPhenotype,
        Stage::FitNullModel,
        Stage::EmitSumstats,
        Stage::EnsureScoreCache,
        Stage::RunScoring,
        Stage::WriteResults,
    ];

    /// Label used in run.json and operator-facing logs. Kept in sync with
    /// the serde rename via `stage_labels_match_serde`.
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

/// On-disk schema version for `run.json`. Bump on any field add/remove
/// in `RunManifest`, `StageRecord`, or `RunOutputs`.
pub const RUN_MANIFEST_VERSION: u32 = 1;

fn default_manifest_version() -> u32 {
    // Pre-versioned manifests on disk are treated as v1 for compatibility.
    1
}

/// `{output_dir}/run.json` — run mode, per-stage status, output paths,
/// cache decisions. `StaarPipeline::run` probes this at startup and
/// skips stages that are both `Completed` here *and* backed by durable
/// disk artifacts. In-memory stages (phenotype, null model, scoring)
/// always rerun.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunManifest {
    /// On-disk schema version. See `RUN_MANIFEST_VERSION`.
    #[serde(default = "default_manifest_version")]
    pub schema_version: u32,
    pub cohort_version: String,
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
    /// Per-stage cache hit/miss log in execution order so operators don't
    /// have to grep stdout to find out why a stage rebuilt.
    #[serde(default)]
    pub cache_decisions: Vec<CacheDecision>,
    /// Resume planner result for this invocation. Set once at the top of
    /// `run()`; `None` on a fresh manifest.
    #[serde(default)]
    pub resume_summary: Option<ResumeSummary>,
}

/// Serializable rendering of `ResumeDecision` for `run.json`.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "outcome")]
pub enum ResumeSummary {
    Fresh,
    Discarded { reason: String },
    Resume { skipped_stages: Vec<Stage> },
}

impl From<&ResumeDecision> for ResumeSummary {
    fn from(decision: &ResumeDecision) -> Self {
        match decision {
            ResumeDecision::Fresh => ResumeSummary::Fresh,
            ResumeDecision::Discarded(reason) => ResumeSummary::Discarded {
                reason: reason.clone(),
            },
            ResumeDecision::Resume { skippable } => ResumeSummary::Resume {
                skipped_stages: skippable.clone(),
            },
        }
    }
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
    /// `config_hash` is computed by the caller so this module doesn't
    /// have to know about `StaarConfig`.
    pub fn new(run_mode: RunMode, trait_name: String, config_hash: String) -> Self {
        Self {
            schema_version: RUN_MANIFEST_VERSION,
            cohort_version: env!("CARGO_PKG_VERSION").to_string(),
            run_mode,
            trait_name,
            config_hash,
            started_unix: now_unix(),
            stages: Vec::new(),
            outputs: RunOutputs::default(),
        }
    }

    /// Most recent stage marked `Failed`, if any. Used by tests today; a
    /// future `--resume` flag will read it.
    #[allow(dead_code)]
    pub fn last_failed_stage(&self) -> Option<Stage> {
        self.stages
            .iter()
            .rev()
            .find(|r| r.status == StageStatus::Failed)
            .map(|r| r.stage)
    }

    /// True if this build can read the on-disk `schema_version`.
    pub fn schema_compatible(&self) -> bool {
        self.schema_version == RUN_MANIFEST_VERSION
    }

    pub fn begin(&mut self, stage: Stage) {
        self.stages.retain(|r| r.stage != stage);
        self.stages.push(StageRecord {
            stage,
            status: StageStatus::InProgress,
            started_unix: now_unix(),
            completed_unix: None,
        });
    }

    pub fn complete(&mut self, stage: Stage) {
        if let Some(rec) = self.stages.iter_mut().rev().find(|r| r.stage == stage) {
            rec.status = StageStatus::Completed;
            rec.completed_unix = Some(now_unix());
        }
    }

    pub fn fail(&mut self, stage: Stage) {
        if let Some(rec) = self.stages.iter_mut().rev().find(|r| r.stage == stage) {
            rec.status = StageStatus::Failed;
            rec.completed_unix = Some(now_unix());
        }
    }

    /// Atomic write to `{output_dir}/run.json` via `store::write_atomic`.
    pub fn write(&self, output_dir: &Path) -> Result<(), CohortError> {
        let path = output_dir.join("run.json");
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| CohortError::Resource(format!("Serialize run manifest: {e}")))?;
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

pub fn now_unix() -> u64 {
    SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Borrowed view of the `StaarConfig` fields that determine results.
/// Lives here so the hashing logic doesn't have to import `StaarConfig`.
pub struct ConfigHashInputs<'a> {
    pub genotypes: &'a Path,
    pub annotations: &'a Path,
    pub phenotype: &'a Path,
    pub trait_names: &'a [String],
    pub covariates: &'a [String],
    pub maf_cutoff: f64,
    pub run_mode: RunMode,
    pub kinship: &'a [PathBuf],
    pub kinship_groups: Option<&'a str>,
}

/// Hash of the inputs that determine results. Wider than
/// `score_cache::cache_key` because the resume planner needs *every*
/// field a stage decision can depend on.
pub fn compute_config_hash(inputs: &ConfigHashInputs<'_>) -> String {
    use sha2::{Digest, Sha256};
    let mut h = Sha256::new();
    h.update(inputs.genotypes.to_string_lossy().as_bytes());
    h.update(b"|");
    h.update(inputs.annotations.to_string_lossy().as_bytes());
    h.update(b"|");
    h.update(inputs.phenotype.to_string_lossy().as_bytes());
    h.update(b"|");
    for t in inputs.trait_names {
        h.update(t.as_bytes());
        h.update(b",");
    }
    h.update(b"|");
    let mut covs = inputs.covariates.to_vec();
    covs.sort();
    for c in &covs {
        h.update(c.as_bytes());
        h.update(b",");
    }
    h.update(b"|");
    h.update(inputs.maf_cutoff.to_le_bytes());
    h.update(format!("{:?}", inputs.run_mode).as_bytes());
    h.update(b"|");
    for k in inputs.kinship {
        h.update(k.to_string_lossy().as_bytes());
        h.update(b",");
    }
    if let Some(g) = inputs.kinship_groups {
        h.update(g.as_bytes());
    }
    format!("{:x}", h.finalize())
}

/// What the resume planner decided to do with a prior `run.json`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ResumeDecision {
    /// No prior manifest on disk.
    Fresh,
    /// Prior manifest existed but was rejected (different config, schema
    /// mismatch, parse error). Carries a one-line reason for the operator.
    Discarded(String),
    /// Prior manifest under the same config; the listed stages can be
    /// short-circuited (their state is durable on disk and doesn't need to
    /// be recomputed in memory). Stages not in `skippable` will rerun.
    Resume {
        skippable: Vec<Stage>,
    },
}

impl ResumeDecision {
    /// Skippable stages on this run. The actual short-circuit happens
    /// inside `EnsureStore` / `EnsureScoreCache` via content fingerprints;
    /// this accessor is for tests and tooling that want the decision
    /// without re-deriving it.
    #[allow(dead_code)]
    pub fn skippable(&self) -> &[Stage] {
        match self {
            ResumeDecision::Resume { skippable } => skippable,
            _ => &[],
        }
    }
}

/// Decide what to do with a prior manifest: discard, resume with
/// short-circuited stages, or treat as fresh. Pure so it's testable
/// without spinning up a `StaarPipeline`.
pub fn plan_resume(prior: Option<RunManifest>, current_hash: &str) -> ResumeDecision {
    let Some(prior) = prior else {
        return ResumeDecision::Fresh;
    };
    if !prior.schema_compatible() {
        return ResumeDecision::Discarded(format!(
            "schema v{} on disk, build expects v{RUN_MANIFEST_VERSION}",
            prior.schema_version
        ));
    }
    if prior.config_hash != current_hash {
        return ResumeDecision::Discarded("config hash differs from prior run".into());
    }
    // Only durable-on-disk stages are skippable. Phenotype/null model/
    // scoring are in-memory and always rerun. Validate stays in the rerun
    // set so column contracts are re-checked under the current tier.
    let durable: &[Stage] = &[Stage::EnsureStore, Stage::EnsureScoreCache];
    let skippable: Vec<Stage> = prior
        .completed_stages()
        .into_iter()
        .filter(|s| durable.contains(s))
        .collect();
    ResumeDecision::Resume { skippable }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dummy_manifest() -> RunManifest {
        RunManifest::new(RunMode::Analyze, "BMI".into(), "abc".into())
    }

    #[test]
    fn run_manifest_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let mut m = dummy_manifest();
        m.begin(Stage::Validate);
        m.complete(Stage::Validate);
        m.begin(Stage::EnsureStore);
        m.complete(Stage::EnsureStore);
        m.outputs.store_dir = Some(PathBuf::from("/tmp/store"));
        m.write(dir.path()).unwrap();

        let loaded = RunManifest::probe(dir.path()).unwrap();
        assert_eq!(loaded.config_hash, m.config_hash);
        assert_eq!(loaded.run_mode, RunMode::Analyze);
        assert_eq!(
            loaded.completed_stages(),
            vec![Stage::Validate, Stage::EnsureStore]
        );
        assert_eq!(loaded.outputs.store_dir, Some(PathBuf::from("/tmp/store")));
    }

    #[test]
    fn run_manifest_fail_marks_stage() {
        let mut m = dummy_manifest();
        m.begin(Stage::FitNullModel);
        m.fail(Stage::FitNullModel);
        let rec = m.stages.last().unwrap();
        assert_eq!(rec.status, StageStatus::Failed);
        assert!(rec.completed_unix.is_some());
    }

    #[test]
    fn cache_decisions_round_trip_through_run_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let mut m = dummy_manifest();
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
        assert_eq!(
            loaded.outputs.cache_decisions[0].outcome,
            CacheOutcome::Miss
        );
        assert_eq!(
            loaded.outputs.cache_decisions[0].artifact,
            ArtifactKind::GenotypeStore
        );
        assert_eq!(
            loaded.outputs.cache_decisions[1].outcome,
            CacheOutcome::Hit
        );
        assert_eq!(
            loaded.outputs.cache_decisions[1].artifact,
            ArtifactKind::ScoreCache
        );
    }

    #[test]
    fn stage_labels_match_serde() {
        // The on-disk run.json key (serde rename_all = snake_case) and the
        // human-facing log label MUST agree. This test catches drift if a
        // future stage is added without updating both sites.
        for stage in Stage::ALL {
            let serialized = serde_json::to_string(&stage).unwrap();
            let serialized = serialized.trim_matches('"').to_string();
            assert_eq!(
                serialized,
                stage.label(),
                "Stage::label() drifted from serde rename for {:?}",
                stage
            );
        }
    }

    #[test]
    fn stage_all_matches_run_order() {
        // Stage::ALL is the canonical execution order. Any stage missing
        // from this constant will cause the resume planner to silently
        // ignore it.
        assert_eq!(Stage::ALL.len(), 8);
        assert_eq!(Stage::ALL[0], Stage::Validate);
        assert_eq!(Stage::ALL[Stage::ALL.len() - 1], Stage::WriteResults);
    }

    #[test]
    fn run_manifest_last_failed_stage() {
        let mut m = dummy_manifest();
        assert_eq!(m.last_failed_stage(), None);

        m.begin(Stage::Validate);
        m.complete(Stage::Validate);
        m.begin(Stage::EnsureStore);
        m.fail(Stage::EnsureStore);
        assert_eq!(m.last_failed_stage(), Some(Stage::EnsureStore));

        // A subsequent failure replaces the answer.
        m.begin(Stage::FitNullModel);
        m.fail(Stage::FitNullModel);
        assert_eq!(m.last_failed_stage(), Some(Stage::FitNullModel));
    }

    #[test]
    fn run_manifest_schema_version_set() {
        let m = dummy_manifest();
        assert_eq!(m.schema_version, RUN_MANIFEST_VERSION);
        assert!(m.schema_compatible());
    }

    #[test]
    fn plan_resume_fresh_when_no_prior() {
        let decision = plan_resume(None, "abc");
        assert_eq!(decision, ResumeDecision::Fresh);
        assert!(decision.skippable().is_empty());
    }

    #[test]
    fn plan_resume_discards_on_config_change() {
        let mut prior = dummy_manifest();
        prior.config_hash = "old-hash".into();
        let decision = plan_resume(Some(prior), "new-hash");
        match decision {
            ResumeDecision::Discarded(reason) => assert!(reason.contains("config hash")),
            other => panic!("expected Discarded, got {other:?}"),
        }
    }

    #[test]
    fn plan_resume_discards_on_schema_mismatch() {
        let mut prior = dummy_manifest();
        prior.schema_version = RUN_MANIFEST_VERSION + 1;
        let hash = prior.config_hash.clone();
        let decision = plan_resume(Some(prior), &hash);
        match decision {
            ResumeDecision::Discarded(reason) => assert!(reason.contains("schema")),
            other => panic!("expected Discarded, got {other:?}"),
        }
    }

    #[test]
    fn plan_resume_only_skips_durable_stages() {
        let mut prior = dummy_manifest();
        for stage in Stage::ALL {
            prior.begin(stage);
            prior.complete(stage);
        }
        let hash = prior.config_hash.clone();
        let decision = plan_resume(Some(prior), &hash);
        let skipped = decision.skippable().to_vec();
        assert!(skipped.contains(&Stage::EnsureStore));
        assert!(skipped.contains(&Stage::EnsureScoreCache));
        assert!(!skipped.contains(&Stage::Validate));
        assert!(!skipped.contains(&Stage::LoadPhenotype));
        assert!(!skipped.contains(&Stage::FitNullModel));
        assert!(!skipped.contains(&Stage::RunScoring));
        assert!(!skipped.contains(&Stage::WriteResults));
    }

    #[test]
    fn resume_summary_serializes_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let mut m = dummy_manifest();
        m.outputs.resume_summary = Some(ResumeSummary::Resume {
            skipped_stages: vec![Stage::EnsureStore, Stage::EnsureScoreCache],
        });
        m.write(dir.path()).unwrap();
        let loaded = RunManifest::probe(dir.path()).unwrap();
        match loaded.outputs.resume_summary {
            Some(ResumeSummary::Resume { skipped_stages }) => {
                assert_eq!(
                    skipped_stages,
                    vec![Stage::EnsureStore, Stage::EnsureScoreCache]
                );
            }
            other => panic!("expected Resume summary, got {other:?}"),
        }
    }

    #[test]
    fn run_manifest_legacy_version_default() {
        // Earlier builds wrote run.json without a schema_version field.
        // The serde default keeps those readable as v1.
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("run.json");
        let legacy = serde_json::json!({
            "cohort_version": "0.0.0",
            "run_mode": "analyze",
            "trait_name": "BMI",
            "config_hash": "abc",
            "started_unix": 0,
            "stages": [],
            "outputs": {}
        });
        std::fs::write(&path, serde_json::to_string(&legacy).unwrap()).unwrap();
        let loaded = RunManifest::probe(dir.path()).unwrap();
        assert_eq!(loaded.schema_version, 1);
        assert!(loaded.schema_compatible());
    }

    #[test]
    fn compute_config_hash_changes_with_inputs() {
        let trait_names = vec!["BMI".to_string()];
        let covs = vec!["age".to_string(), "sex".to_string()];
        let kinship: Vec<PathBuf> = Vec::new();
        let inputs = ConfigHashInputs {
            genotypes: Path::new("/tmp/g.vcf.gz"),
            annotations: Path::new("/tmp/a"),
            phenotype: Path::new("/tmp/p.tsv"),
            trait_names: &trait_names,
            covariates: &covs,
            maf_cutoff: 0.01,
            run_mode: RunMode::Analyze,
            kinship: &kinship,
            kinship_groups: None,
        };
        let h1 = compute_config_hash(&inputs);

        let inputs2 = ConfigHashInputs {
            maf_cutoff: 0.05,
            ..ConfigHashInputs {
                genotypes: Path::new("/tmp/g.vcf.gz"),
                annotations: Path::new("/tmp/a"),
                phenotype: Path::new("/tmp/p.tsv"),
                trait_names: &trait_names,
                covariates: &covs,
                maf_cutoff: 0.01,
                run_mode: RunMode::Analyze,
                kinship: &kinship,
                kinship_groups: None,
            }
        };
        let h2 = compute_config_hash(&inputs2);
        assert_ne!(h1, h2);

        // Covariate order doesn't matter (sorted before hashing).
        let covs_reversed = vec!["sex".to_string(), "age".to_string()];
        let inputs3 = ConfigHashInputs {
            covariates: &covs_reversed,
            ..ConfigHashInputs {
                genotypes: Path::new("/tmp/g.vcf.gz"),
                annotations: Path::new("/tmp/a"),
                phenotype: Path::new("/tmp/p.tsv"),
                trait_names: &trait_names,
                covariates: &covs,
                maf_cutoff: 0.01,
                run_mode: RunMode::Analyze,
                kinship: &kinship,
                kinship_groups: None,
            }
        };
        assert_eq!(h1, compute_config_hash(&inputs3));
    }
}
