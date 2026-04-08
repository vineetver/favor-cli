use std::path::PathBuf;

use crate::config::{Config, Environment, Tier};

use super::types::{
    ArtifactKind, FormError, FormField, FormSchema, FormValues, PathKind, RunRequest, SessionCtx,
    SetupConfig, StageId,
};
use super::Stage;

pub struct SetupStage;

const MEMORY_PRESETS: &[&str] = &["auto", "8GB", "16GB", "32GB", "64GB", "128GB", "256GB"];

pub const OPTIONAL_PACK_IDS: &[&str] = &[
    "eqtl",
    "eqtl-catalogue",
    "sc-eqtl",
    "regulatory",
    "enhancer-gene",
    "tissue-scores",
    "pgs",
    "genotypes",
];

const NO_INPUTS: &[ArtifactKind] = &[];
const NO_OUTPUTS: &[ArtifactKind] = &[];

impl Stage for SetupStage {
    fn id(&self) -> StageId {
        StageId("setup")
    }
    fn label(&self) -> &'static str {
        "Setup"
    }
    fn inputs(&self) -> &'static [ArtifactKind] {
        NO_INPUTS
    }
    fn outputs(&self) -> &'static [ArtifactKind] {
        NO_OUTPUTS
    }

    fn form_schema(&self, _ctx: &SessionCtx) -> FormSchema {
        let cfg = Config::load().unwrap_or_default();
        let root = if cfg.data.root_dir.is_empty() {
            Config::default_root_dir()
        } else {
            PathBuf::from(&cfg.data.root_dir)
        };
        let memory_default = cfg
            .resources
            .memory_budget
            .as_deref()
            .and_then(|m| MEMORY_PRESETS.iter().copied().find(|p| *p == m))
            .unwrap_or("auto");
        let env_default = match cfg.resources.environment {
            None => "auto",
            Some(Environment::Hpc) => "hpc",
            Some(Environment::Workstation) => "workstation",
        };
        let tier_default = match cfg.data.tier {
            Tier::Base => "base",
            Tier::Full => "full",
        };
        FormSchema {
            fields: vec![
                FormField::Choice {
                    id: "tier",
                    label: "tier",
                    options: &["base", "full"],
                    default: Some(tier_default),
                },
                FormField::Path {
                    id: "root",
                    label: "data root",
                    kind: PathKind::Dir,
                    default: Some(root),
                },
            ],
            advanced: vec![
                FormField::Choice {
                    id: "env",
                    label: "env",
                    options: &["auto", "hpc", "workstation"],
                    default: Some(env_default),
                },
                FormField::Choice {
                    id: "memory",
                    label: "memory",
                    options: MEMORY_PRESETS,
                    default: Some(memory_default),
                },
                FormField::MultiSelect {
                    id: "packs",
                    label: "packs",
                    options: OPTIONAL_PACK_IDS,
                    default: &[],
                },
            ],
        }
    }

    fn build_command(&self, values: &FormValues) -> Result<RunRequest, FormError> {
        let tier = match values.choice("tier") {
            Some("full") => Tier::Full,
            _ => Tier::Base,
        };
        let root_dir = values
            .path("root")
            .cloned()
            .ok_or(FormError::Missing("root"))?;
        if !root_dir.is_dir() {
            return Err(FormError::Invalid {
                field: "root",
                reason: format!("not a directory: {}", root_dir.display()),
            });
        }
        let environment = match values.choice("env") {
            Some("hpc") => Some(Environment::Hpc),
            Some("workstation") => Some(Environment::Workstation),
            _ => None,
        };
        let memory_budget = match values.choice("memory") {
            Some(v) if v != "auto" && !v.is_empty() => Some(v.to_string()),
            _ => None,
        };
        let packs = values.multi("packs").cloned().unwrap_or_default();
        Ok(RunRequest::Setup(SetupConfig {
            tier,
            root_dir,
            packs,
            environment,
            memory_budget,
        }))
    }
}
