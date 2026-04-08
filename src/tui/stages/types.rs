use std::collections::HashMap;
use std::path::PathBuf;

use crate::commands::{AnnotateConfig, IngestConfig, MetaStaarConfig};
use crate::config::Tier;
use crate::staar::pipeline::StaarConfig;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct StageId(pub &'static str);

impl StageId {
    pub const fn as_str(&self) -> &'static str {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageGroup {
    Ingest,
    Annotate,
    Test,
    Meta,
    #[allow(dead_code)]
    Report,
}

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
    pub fn is_edited(&self) -> bool {
        matches!(self, Self::Edited(_))
    }
}

#[derive(Debug, Clone)]
pub enum FieldData {
    Path(FieldValue<Option<PathBuf>>),
    Paths(FieldValue<Vec<PathBuf>>),
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
    pub fn get(&self, id: &str) -> Option<&FieldData> {
        self.values.get(id)
    }

    pub fn path(&self, id: &str) -> Option<&PathBuf> {
        match self.values.get(id)? {
            FieldData::Path(v) => v.get().as_ref(),
            _ => None,
        }
    }

    pub fn paths(&self, id: &str) -> Option<&Vec<PathBuf>> {
        match self.values.get(id)? {
            FieldData::Paths(v) => Some(v.get()),
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
    pub data_root: &'a std::path::Path,
    pub tier: Tier,
    pub artifacts: &'a [PathBuf],
}

pub enum RunRequest {
    Ingest(IngestConfig),
    Annotate(AnnotateConfig),
    Staar(Box<StaarConfig>),
    MetaStaar(MetaStaarConfig),
}
