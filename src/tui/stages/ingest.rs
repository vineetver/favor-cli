use std::path::{Path, PathBuf};

use crate::commands::IngestConfig;

use super::types::{
    ArtifactKind, FormError, FormField, FormSchema, FormValues, PathKind, RunRequest, SessionCtx,
    StageId,
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
    fn inputs(&self) -> &'static [ArtifactKind] {
        INPUTS
    }
    fn outputs(&self) -> &'static [ArtifactKind] {
        OUTPUTS
    }

    fn form_schema(&self, ctx: &SessionCtx) -> FormSchema {
        let focused_input = ctx.focused.map(|p| p.to_path_buf());
        let default_output = focused_input.as_deref().map(default_ingest_output);
        FormSchema {
            fields: vec![
                FormField::Path {
                    id: "inputs",
                    label: "inputs",
                    kind: PathKind::Any,
                    default: focused_input,
                },
                FormField::Path {
                    id: "output",
                    label: "output",
                    kind: PathKind::Dir,
                    default: default_output,
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
            .path("inputs")
            .map(|p| vec![p.clone()])
            .ok_or(FormError::Missing("inputs"))?;
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

fn default_ingest_output(first: &Path) -> PathBuf {
    let stem = first
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    let stem = stem
        .strip_suffix(".vcf")
        .or_else(|| stem.strip_suffix(".tsv"))
        .or_else(|| stem.strip_suffix(".csv"))
        .unwrap_or(&stem);
    let stem = stem
        .split("_b0_")
        .next()
        .or_else(|| stem.split("_b0.").next())
        .unwrap_or(stem);
    first
        .parent()
        .unwrap_or(first)
        .join(format!("{stem}.ingested"))
}
