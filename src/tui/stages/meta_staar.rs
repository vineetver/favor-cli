use std::path::PathBuf;

use crate::commands::{parse_mask_categories, MetaStaarConfig};

use super::types::{
    ArtifactKind, FormError, FormField, FormSchema, FormValues, PathKind, RunRequest, SessionCtx,
    StageId,
};
use super::Stage;

pub struct MetaStaarStage;

const INPUTS: &[ArtifactKind] = &[ArtifactKind::SumStats];
const OUTPUTS: &[ArtifactKind] = &[ArtifactKind::MetaStaarResults];

const MASK_DEFAULT: &[&str] = &["coding"];

impl Stage for MetaStaarStage {
    fn id(&self) -> StageId {
        StageId("meta-staar")
    }
    fn label(&self) -> &'static str {
        "Meta-STAAR"
    }
    fn inputs(&self) -> &'static [ArtifactKind] {
        INPUTS
    }
    fn outputs(&self) -> &'static [ArtifactKind] {
        OUTPUTS
    }

    fn form_schema(&self, ctx: &SessionCtx) -> FormSchema {
        let focused = ctx.focused.map(|p| p.to_path_buf());
        FormSchema {
            fields: vec![
                FormField::Path {
                    id: "studies",
                    label: "studies",
                    kind: PathKind::Any,
                    default: focused,
                },
                FormField::MultiSelect {
                    id: "masks",
                    label: "masks",
                    default: MASK_DEFAULT,
                },
                FormField::Path {
                    id: "output_dir",
                    label: "output dir",
                    kind: PathKind::Dir,
                    default: None,
                },
            ],
            advanced: vec![
                FormField::Number {
                    id: "maf_cutoff",
                    label: "MAF cutoff",
                    default: Some(0.01),
                },
                FormField::Number {
                    id: "window_size",
                    label: "window size",
                    default: Some(2000.0),
                },
            ],
        }
    }

    fn build_command(&self, values: &FormValues) -> Result<RunRequest, FormError> {
        let study_dirs: Vec<PathBuf> = values
            .path("studies")
            .map(|p| vec![p.clone()])
            .ok_or(FormError::Missing("studies"))?;
        let masks_raw: Vec<String> = values
            .multi("masks")
            .cloned()
            .unwrap_or_else(|| MASK_DEFAULT.iter().map(|s| s.to_string()).collect());
        let mask_categories = parse_mask_categories(&masks_raw).map_err(|e| FormError::Invalid {
            field: "masks",
            reason: e.to_string(),
        })?;
        let maf_cutoff = values.number("maf_cutoff").unwrap_or(0.01);
        let window_size = values.number("window_size").unwrap_or(2000.0) as u32;
        let output_dir = values
            .path("output_dir")
            .cloned()
            .unwrap_or_else(|| PathBuf::from("meta_staar_out"));
        Ok(RunRequest::MetaStaar(MetaStaarConfig {
            study_dirs,
            mask_categories,
            maf_cutoff,
            window_size,
            output_dir,
        }))
    }
}
