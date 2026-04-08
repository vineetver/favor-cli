use std::collections::HashMap;
use std::path::PathBuf;

use crate::commands::parse_mask_categories;
use crate::staar::masks::ScangParams;
use crate::staar::pipeline::StaarConfig;

use super::types::{
    ArtifactKind, FormError, FormField, FormSchema, FormValues, PathKind, RunRequest, SessionCtx,
    StageId,
};
use super::Stage;

pub struct StaarStage;

const INPUTS: &[ArtifactKind] = &[ArtifactKind::AnnotatedSet, ArtifactKind::GenotypeStore];
const OUTPUTS: &[ArtifactKind] = &[ArtifactKind::StaarResults, ArtifactKind::SumStats];

const MASK_DEFAULT: &[&str] = &["coding"];

impl Stage for StaarStage {
    fn id(&self) -> StageId {
        StageId("staar")
    }
    fn label(&self) -> &'static str {
        "STAAR"
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
                    id: "genotypes",
                    label: "genotypes",
                    kind: PathKind::Any,
                    default: focused.clone(),
                },
                FormField::Path {
                    id: "phenotype",
                    label: "phenotype",
                    kind: PathKind::File,
                    default: None,
                },
                FormField::Path {
                    id: "annotations",
                    label: "annotations",
                    kind: PathKind::Any,
                    default: None,
                },
                FormField::Text {
                    id: "trait_name",
                    label: "trait",
                    default: None,
                },
                FormField::Text {
                    id: "covariates",
                    label: "covariates",
                    default: None,
                },
                FormField::MultiSelect {
                    id: "masks",
                    label: "masks",
                    default: MASK_DEFAULT,
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
                FormField::Toggle {
                    id: "spa",
                    label: "SPA",
                    default: false,
                },
                FormField::Toggle {
                    id: "emit_sumstats",
                    label: "emit sumstats",
                    default: false,
                },
                FormField::Toggle {
                    id: "rebuild_store",
                    label: "rebuild store",
                    default: false,
                },
                FormField::Path {
                    id: "known_loci",
                    label: "known loci",
                    kind: PathKind::File,
                    default: None,
                },
                FormField::Path {
                    id: "output_dir",
                    label: "output dir",
                    kind: PathKind::Dir,
                    default: None,
                },
                FormField::Path {
                    id: "store_dir",
                    label: "store dir",
                    kind: PathKind::Dir,
                    default: None,
                },
            ],
        }
    }

    fn build_command(&self, values: &FormValues) -> Result<RunRequest, FormError> {
        let genotypes = values
            .path("genotypes")
            .cloned()
            .ok_or(FormError::Missing("genotypes"))?;
        let phenotype = values
            .path("phenotype")
            .cloned()
            .ok_or(FormError::Missing("phenotype"))?;
        let annotations = values
            .path("annotations")
            .cloned()
            .ok_or(FormError::Missing("annotations"))?;
        let trait_name = values
            .text("trait_name")
            .ok_or(FormError::Missing("trait_name"))?
            .to_string();
        let covariates: Vec<String> = values
            .text("covariates")
            .map(|s| {
                s.split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(String::from)
                    .collect()
            })
            .unwrap_or_default();
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
        let spa = values.toggle("spa").unwrap_or(false);
        let emit_sumstats = values.toggle("emit_sumstats").unwrap_or(false);
        let rebuild_store = values.toggle("rebuild_store").unwrap_or(false);
        let known_loci = values.path("known_loci").cloned();
        let output_dir = values
            .path("output_dir")
            .cloned()
            .unwrap_or_else(|| PathBuf::from("staar_out"));
        let store_dir = values
            .path("store_dir")
            .cloned()
            .unwrap_or_else(|| output_dir.join("store"));

        Ok(RunRequest::Staar(Box::new(StaarConfig {
            genotypes,
            phenotype,
            annotations,
            trait_names: vec![trait_name],
            covariates,
            mask_categories,
            maf_cutoff,
            window_size,
            individual: false,
            spa,
            ancestry_col: None,
            ai_base_tests: 0,
            ai_seed: 0,
            scang_params: ScangParams::default(),
            kinship: Vec::new(),
            kinship_groups: None,
            known_loci,
            emit_sumstats,
            rebuild_store,
            column_map: HashMap::new(),
            output_dir,
            store_dir,
        })))
    }
}
