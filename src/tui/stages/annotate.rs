use std::path::{Path, PathBuf};

use crate::commands::AnnotateConfig;
use crate::config::Tier;

use super::types::{
    ArtifactKind, FormError, FormField, FormSchema, FormValues, PathKind, RunRequest, SessionCtx,
    StageId,
};
use super::Stage;

pub struct AnnotateStage;

const INPUTS: &[ArtifactKind] = &[ArtifactKind::IngestedSet];
const OUTPUTS: &[ArtifactKind] = &[ArtifactKind::AnnotatedSet];

impl Stage for AnnotateStage {
    fn id(&self) -> StageId {
        StageId("annotate")
    }
    fn label(&self) -> &'static str {
        "Annotate"
    }
    fn inputs(&self) -> &'static [ArtifactKind] {
        INPUTS
    }
    fn outputs(&self) -> &'static [ArtifactKind] {
        OUTPUTS
    }

    fn form_schema(&self, ctx: &SessionCtx) -> FormSchema {
        let focused_input = ctx.focused.map(|p| p.to_path_buf());
        let default_output = focused_input.as_deref().map(default_annotate_output);
        FormSchema {
            fields: vec![
                FormField::Path {
                    id: "input",
                    label: "input",
                    kind: PathKind::Dir,
                    default: focused_input,
                },
                FormField::Path {
                    id: "output",
                    label: "output",
                    kind: PathKind::Dir,
                    default: default_output,
                },
                FormField::Choice {
                    id: "tier",
                    label: "tier",
                    options: &["base", "full"],
                    default: Some(ctx.tier.as_str()),
                },
            ],
            advanced: vec![],
        }
    }

    fn build_command(&self, values: &FormValues) -> Result<RunRequest, FormError> {
        let input = values
            .path("input")
            .cloned()
            .ok_or(FormError::Missing("input"))?;
        let output = values
            .path("output")
            .cloned()
            .ok_or(FormError::Missing("output"))?;
        let tier = match values.choice("tier") {
            Some("full") => Tier::Full,
            _ => Tier::Base,
        };
        Ok(RunRequest::Annotate(AnnotateConfig {
            input,
            output,
            tier,
        }))
    }
}

fn default_annotate_output(input: &Path) -> PathBuf {
    let name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    let stem = name.strip_suffix(".ingested").unwrap_or(&name);
    input
        .parent()
        .unwrap_or(input)
        .join(format!("{stem}.annotated"))
}
