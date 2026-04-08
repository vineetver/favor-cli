use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::commands::{AnnotateConfig, IngestConfig, MetaStaarConfig};
use crate::config::{Environment, Tier};
use crate::staar::pipeline::StaarConfig;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StageId(pub &'static str);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ArtifactKind {
    VcfDir,
    IngestedSet,
    AnnotatedSet,
    GenotypeStore,
    StaarResults,
    SumStats,
    MetaStaarResults,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathKind {
    File,
    Dir,
    Any,
}

#[derive(Debug, Clone)]
pub enum FormField {
    Path {
        id: &'static str,
        label: &'static str,
        kind: PathKind,
        default: Option<PathBuf>,
    },
    Choice {
        id: &'static str,
        label: &'static str,
        options: &'static [&'static str],
        default: Option<&'static str>,
    },
    Toggle {
        id: &'static str,
        label: &'static str,
        default: bool,
    },
    Text {
        id: &'static str,
        label: &'static str,
        default: Option<String>,
    },
    Number {
        id: &'static str,
        label: &'static str,
        default: Option<f64>,
    },
    MultiSelect {
        id: &'static str,
        label: &'static str,
        options: &'static [&'static str],
        default: &'static [&'static str],
    },
}

impl FormField {
    pub fn id(&self) -> &'static str {
        match self {
            Self::Path { id, .. }
            | Self::Choice { id, .. }
            | Self::Toggle { id, .. }
            | Self::Text { id, .. }
            | Self::Number { id, .. }
            | Self::MultiSelect { id, .. } => id,
        }
    }
}

#[derive(Debug, Clone)]
pub struct FormSchema {
    pub fields: Vec<FormField>,
    pub advanced: Vec<FormField>,
}

#[derive(Debug, Clone)]
pub enum FieldValue<T> {
    Inferred(T),
    Edited(T),
}

impl<T> FieldValue<T> {
    pub fn get(&self) -> &T {
        match self {
            Self::Inferred(v) | Self::Edited(v) => v,
        }
    }
}

#[derive(Debug, Clone)]
pub enum FieldData {
    Path(FieldValue<Option<PathBuf>>),
    Choice(FieldValue<Option<String>>),
    Toggle(FieldValue<bool>),
    Text(FieldValue<String>),
    Number(FieldValue<f64>),
    Multi(FieldValue<Vec<String>>),
}

#[derive(Debug, Clone, Default)]
pub struct FormValues {
    pub values: HashMap<&'static str, FieldData>,
}

impl FormValues {
    pub fn path(&self, id: &str) -> Option<&PathBuf> {
        match self.values.get(id)? {
            FieldData::Path(v) => v.get().as_ref(),
            _ => None,
        }
    }

    pub fn toggle(&self, id: &str) -> Option<bool> {
        match self.values.get(id)? {
            FieldData::Toggle(v) => Some(*v.get()),
            _ => None,
        }
    }

    pub fn choice(&self, id: &str) -> Option<&str> {
        match self.values.get(id)? {
            FieldData::Choice(v) => v.get().as_deref(),
            _ => None,
        }
    }

    pub fn text(&self, id: &str) -> Option<&str> {
        match self.values.get(id)? {
            FieldData::Text(v) => Some(v.get().as_str()),
            _ => None,
        }
    }

    pub fn number(&self, id: &str) -> Option<f64> {
        match self.values.get(id)? {
            FieldData::Number(v) => Some(*v.get()),
            _ => None,
        }
    }

    pub fn multi(&self, id: &str) -> Option<&Vec<String>> {
        match self.values.get(id)? {
            FieldData::Multi(v) => Some(v.get()),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub enum FormError {
    Missing(&'static str),
    Invalid { field: &'static str, reason: String },
}

pub struct SessionCtx<'a> {
    pub data_root: &'a Path,
    pub tier: Tier,
    pub focused: Option<&'a Path>,
}

#[derive(Debug, Clone)]
pub struct SetupConfig {
    pub tier: Tier,
    pub root_dir: PathBuf,
    pub packs: Vec<String>,
    pub environment: Option<Environment>,
    pub memory_budget: Option<String>,
}

pub enum RunRequest {
    Ingest(IngestConfig),
    Annotate(AnnotateConfig),
    Staar(Box<StaarConfig>),
    MetaStaar(MetaStaarConfig),
    Setup(SetupConfig),
}

impl RunRequest {
    pub fn description(&self) -> String {
        match self {
            RunRequest::Ingest(cfg) => {
                let first = cfg
                    .inputs
                    .first()
                    .map(|p| p.file_name().unwrap_or_default().to_string_lossy().into_owned())
                    .unwrap_or_default();
                if cfg.inputs.len() > 1 {
                    format!("Ingesting {first} (+{} more)", cfg.inputs.len() - 1)
                } else {
                    format!("Ingesting {first}")
                }
            }
            RunRequest::Annotate(cfg) => {
                let name = cfg.input.file_name().unwrap_or_default().to_string_lossy();
                format!("Annotating {name} ({} tier)", cfg.tier.as_str())
            }
            RunRequest::Staar(cfg) => {
                let trait_name = cfg.trait_names.first().cloned().unwrap_or_default();
                format!("STAAR {trait_name}")
            }
            RunRequest::MetaStaar(cfg) => {
                format!("Meta-STAAR ({} studies)", cfg.study_dirs.len())
            }
            RunRequest::Setup(cfg) => format!("Setup → {}", cfg.root_dir.display()),
        }
    }

    pub fn expected_artifact(&self) -> PathBuf {
        match self {
            RunRequest::Ingest(cfg) => cfg.output.clone(),
            RunRequest::Annotate(cfg) => cfg.output.clone(),
            RunRequest::Staar(cfg) => cfg.output_dir.clone(),
            RunRequest::MetaStaar(cfg) => cfg.output_dir.clone(),
            RunRequest::Setup(cfg) => cfg.root_dir.clone(),
        }
    }
}
