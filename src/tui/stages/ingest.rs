use std::path::PathBuf;

use crate::commands::IngestConfig;

use super::types::{
    ArtifactKind, FormError, FormField, FormSchema, FormValues, PathKind, RunRequest, SessionCtx,
    StageGroup, StageId,
};
use super::Stage;

pub struct IngestStage;

const INPUTS: &[ArtifactKind] = &[ArtifactKind::VcfDir];
const OUTPUTS: &[ArtifactKind] = &[ArtifactKind::IngestedSet];

impl Stage for IngestStage {
    fn id(&self) -> StageId {
        StageId("ingest")
    }
    fn label(&self) -> &'static str {
        "Ingest VCF"
    }
    fn group(&self) -> StageGroup {
        StageGroup::Ingest
    }
    fn inputs(&self) -> &'static [ArtifactKind] {
        INPUTS
    }
    fn outputs(&self) -> &'static [ArtifactKind] {
        OUTPUTS
    }

    fn form_schema(&self, _ctx: &SessionCtx) -> FormSchema {
        FormSchema {
            fields: vec![
                FormField::Path {
                    id: "inputs",
                    label: "inputs",
                    kind: PathKind::Any,
                    default: None,
                },
                FormField::Path {
                    id: "output",
                    label: "output",
                    kind: PathKind::File,
                    default: None,
                },
            ],
            advanced: vec![
                FormField::Toggle {
                    id: "emit_sql",
                    label: "emit SQL",
                    default: false,
                },
                FormField::Choice {
                    id: "build",
                    label: "build",
                    options: &["auto", "hg38", "hg19"],
                    default: Some("auto"),
                },
            ],
        }
    }

    fn build_command(&self, values: &FormValues) -> Result<RunRequest, FormError> {
        let inputs: Vec<PathBuf> = values
            .paths("inputs")
            .cloned()
            .or_else(|| values.path("inputs").map(|p| vec![p.clone()]))
            .ok_or(FormError::Missing("inputs"))?;
        if inputs.is_empty() {
            return Err(FormError::Missing("inputs"));
        }
        let output = values
            .path("output")
            .cloned()
            .ok_or(FormError::Missing("output"))?;
        let emit_sql = values.toggle("emit_sql").unwrap_or(false);
        let build_override = match values.choice("build") {
            Some("hg38") => Some(crate::cli::GenomeBuild::Hg38),
            Some("hg19") => Some(crate::cli::GenomeBuild::Hg19),
            _ => None,
        };
        Ok(RunRequest::Ingest(IngestConfig {
            inputs,
            output,
            emit_sql,
            build_override,
        }))
    }
}
