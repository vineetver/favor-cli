use std::path::PathBuf;

use serde_json::json;

use crate::column::Col;
use crate::commands::{self, EnrichConfig};
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::output::Output;
use crate::runtime::Engine;
use crate::store::annotation::TissueDb;
use crate::store::list::{parquet_row_count, AnnotatedSet};

pub fn build_config(
    input: PathBuf,
    tissue: String,
    output_path: Option<PathBuf>,
) -> Result<EnrichConfig, CohortError> {
    if !input.exists() {
        return Err(CohortError::Input(format!(
            "Annotated variant set not found: '{}'. Run `cohort annotate` first.",
            input.display()
        )));
    }

    let output = output_path.unwrap_or_else(|| {
        let name = input.file_name().unwrap_or_default().to_string_lossy();
        let stem = name
            .strip_suffix(".annotated")
            .or_else(|| name.strip_suffix("/"))
            .unwrap_or(&name);
        input
            .parent()
            .unwrap_or(&input)
            .join(format!("{stem}.enriched"))
    });

    Ok(EnrichConfig {
        input,
        output,
        tissue_name: tissue,
    })
}

pub fn run_enrich(
    engine: &Engine,
    config: &EnrichConfig,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    if dry_run {
        return emit_dry_run(engine, config, out);
    }

    let annotated = AnnotatedSet::open(&config.input)?;
    annotated.supports(&[Col::Vid])?;
    let registry = engine.annotation_registry()?;
    let tissue_db = registry.open_tissue("favor-tissue")?;
    let tissue_dir = engine.config().tissue_dir();
    let df = engine.df();
    let resources = engine.resources();

    out.status(&format!("Input: {} variants", annotated.variant_count()));
    out.status(&format!(
        "Resources: {}, {} threads ({})",
        resources.memory_human(),
        resources.threads,
        resources.environment()
    ));

    let resolved = resolve_tissue(df, &tissue_dir, &config.tissue_name)?;
    out.status(&format!(
        "Tissue '{}' -> {} subtissues",
        config.tissue_name,
        resolved.len()
    ));

    std::fs::create_dir_all(&config.output).map_err(|e| {
        CohortError::Resource(format!("Cannot create '{}': {e}", config.output.display()))
    })?;

    let tables_written = run_enrichment(
        &annotated,
        &tissue_db,
        &resolved,
        df,
        &config.output,
        out,
    )?;
    write_meta(
        &tables_written,
        &resolved,
        config,
        annotated.variant_count() as i64,
        out,
    );
    report_result(
        &tables_written,
        config,
        annotated.variant_count() as i64,
        out,
    );
    Ok(())
}

fn emit_dry_run(
    engine: &Engine,
    config: &EnrichConfig,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let annotated = AnnotatedSet::open(&config.input)?;
    let tissue_db = TissueDb::open(&engine.config().tissue_dir())?;
    let available_tables: Vec<&str> = tissue_db
        .available_tables()
        .iter()
        .map(|t| t.display_name())
        .collect();
    let plan = commands::DryRunPlan {
        command: "enrich".into(),
        inputs: json!({
            "file": config.input.to_string_lossy(),
            "variant_count": annotated.variant_count(),
            "tissue": config.tissue_name,
            "available_tables": available_tables,
        }),
        memory: commands::MemoryEstimate::default_estimate(),
        runtime: None,
        output_path: config.output.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
}

fn resolve_tissue(
    engine: &DfEngine,
    tissue_dir: &std::path::Path,
    tissue_query: &str,
) -> Result<Vec<String>, CohortError> {
    let tissue_vocab_path = tissue_dir.join("reference/tissue_vocab.parquet");
    if !tissue_vocab_path.exists() {
        return Ok(vec![tissue_query.to_string()]);
    }

    engine.register_parquet_file("_tissue_vocab", &tissue_vocab_path)?;

    let escaped = tissue_query.replace('\'', "''");
    let resolved = engine.query_strings(&format!(
        "SELECT DISTINCT tissue_norm FROM _tissue_vocab \
         WHERE tissue_group ILIKE '%{escaped}%' \
            OR tissue_norm ILIKE '%{escaped}%' \
            OR tissue_raw ILIKE '%{escaped}%'"
    ))?;

    if resolved.is_empty() {
        let groups = engine
            .query_strings("SELECT DISTINCT tissue_group FROM _tissue_vocab ORDER BY tissue_group")
            .unwrap_or_default();
        return Err(CohortError::Input(format!(
            "Unknown tissue '{}'. Available groups: {}",
            tissue_query,
            groups.join(", ")
        )));
    }

    Ok(resolved)
}

fn run_enrichment(
    annotated: &AnnotatedSet,
    tissue_db: &TissueDb,
    resolved_tissues: &[String],
    engine: &DfEngine,
    output_dir: &std::path::Path,
    out: &dyn Output,
) -> Result<Vec<(String, i64)>, CohortError> {
    let tissue_filter: String = resolved_tissues
        .iter()
        .map(|t| format!("'{}'", t.replace('\'', "''")))
        .collect::<Vec<_>>()
        .join(", ");

    engine.register_parquet_dir("_input", annotated.root())?;
    let annotated_out = output_dir.join("annotated.parquet");
    engine.execute(&format!(
        "COPY (SELECT * FROM _input) TO '{}' STORED AS PARQUET OPTIONS (compression 'zstd(4)')",
        annotated_out.display(),
    ))?;
    out.status(&format!(
        "  annotated.parquet ({} variants)",
        annotated.variant_count()
    ));

    engine.execute("CREATE TABLE _input_vids AS SELECT DISTINCT vid FROM _input")?;

    let mut tables_written: Vec<(String, i64)> = Vec::new();

    for &table in tissue_db.available_tables() {
        let table_path = tissue_db.table_path(table);
        let table_name = format!("_enrich_{}", table.display_name());
        if engine
            .register_parquet_dir(&table_name, &table_path)
            .is_err()
        {
            continue;
        }

        let out_path = output_dir.join(format!("{}.parquet", table.display_name()));

        let where_clause = if table.has_tissue_filter() {
            format!(
                "WHERE t.vid IN (SELECT vid FROM _input_vids) \
                     AND t.tissue_name IN ({tissue_filter})"
            )
        } else {
            "WHERE t.vid IN (SELECT vid FROM _input_vids)".to_string()
        };

        let sql = format!(
            "COPY (\
                SELECT t.* EXCEPT(chrom_id) \
                FROM {table_name} t \
                {where_clause}\
            ) TO '{out_path}' STORED AS PARQUET OPTIONS (compression 'zstd(4)')",
            out_path = out_path.display(),
        );

        out.status(&format!("  {}: joining...", table.display_name()));
        if let Err(e) = engine.execute(&sql) {
            out.warn(&format!("  {}: skipped ({e})", table.display_name()));
            continue;
        }

        let row_count = if out_path.exists() {
            parquet_row_count(&out_path).unwrap_or(0) as i64
        } else {
            0
        };

        if row_count == 0 {
            let _ = std::fs::remove_file(&out_path);
        } else {
            out.status(&format!(
                "  {}.parquet ({} rows)",
                table.display_name(),
                row_count
            ));
            tables_written.push((table.display_name().to_string(), row_count));
        }
    }

    Ok(tables_written)
}

fn write_meta(
    tables_written: &[(String, i64)],
    resolved_tissues: &[String],
    config: &EnrichConfig,
    input_count: i64,
    out: &dyn Output,
) {
    let meta = json!({
        "cohort_enrich_version": 2,
        "source": config.input.to_string_lossy(),
        "tissue": config.tissue_name,
        "resolved_tissues": resolved_tissues,
        "output_dir": config.output.to_string_lossy(),
        "tables": tables_written.iter()
            .map(|(name, rows)| json!({"name": name, "rows": rows}))
            .collect::<Vec<_>>(),
        "input_count": input_count,
        "join_key": "vid",
        "usage": "SELECT a.*, e.* FROM 'annotated.parquet' a INNER JOIN 'eqtl.parquet' e ON a.vid = e.vid",
    });
    let meta_path = config.output.join("enriched.meta.json");
    if let Ok(json_str) = serde_json::to_string_pretty(&meta) {
        let _ = std::fs::write(&meta_path, json_str);
    }
    let _ = out;
}

fn report_result(
    tables_written: &[(String, i64)],
    config: &EnrichConfig,
    input_count: i64,
    out: &dyn Output,
) {
    if tables_written.is_empty() {
        out.warn("No enrichment data found for these variants in this tissue.");
    } else {
        out.success(&format!(
            "Enriched -> {} ({} tables)",
            config.output.display(),
            tables_written.len(),
        ));
    }

    out.result_json(&json!({
        "status": "ok",
        "output_dir": config.output.to_string_lossy(),
        "tables": tables_written.iter()
            .map(|(name, rows)| json!({"name": name, "rows": rows}))
            .collect::<Vec<_>>(),
        "tissue": config.tissue_name,
        "input_count": input_count,
    }));
}
