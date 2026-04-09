use std::path::PathBuf;

use serde_json::json;

use crate::cli::AttachKind;
use crate::config::Tier;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::store::annotation::refs::{AnnotationKind, AnnotationRegistry};

pub fn attach(
    engine: &Engine,
    name: &str,
    kind: AttachKind,
    path: PathBuf,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let kind = match kind {
        AttachKind::FavorBase => AnnotationKind::Favor { tier: Tier::Base },
        AttachKind::FavorFull => AnnotationKind::Favor { tier: Tier::Full },
        AttachKind::Tissue => AnnotationKind::Tissue,
    };
    let refs_path = engine.store().layout().annotations_refs();
    AnnotationRegistry::attach(&refs_path, name, kind, path.clone())?;
    out.success(&format!(
        "Attached '{}' -> {} ({:?})",
        name,
        path.display(),
        kind,
    ));
    out.result_json(&json!({
        "name": name,
        "kind": kind,
        "path": path.to_string_lossy(),
        "refs_file": refs_path.to_string_lossy(),
    }));
    Ok(())
}
