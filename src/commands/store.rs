use serde_json::json;

use crate::commands::human_bytes;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;

pub fn gc(engine: &Engine, out: &dyn Output) -> Result<(), CohortError> {
    let live = engine.store().list_cohorts()?;
    let summary = engine.store().cache().prune_orphans(&live)?;

    out.success(&format!(
        "store gc: removed {} score cache(s), {} lookup index(es), freed {}",
        summary.removed_score_caches,
        summary.removed_lookup_indexes,
        human_bytes(summary.bytes_freed),
    ));
    out.result_json(&json!({
        "live_cohorts": live.iter().map(|c| c.as_str()).collect::<Vec<_>>(),
        "removed_score_caches": summary.removed_score_caches,
        "removed_lookup_indexes": summary.removed_lookup_indexes,
        "bytes_freed": summary.bytes_freed,
    }));
    Ok(())
}
