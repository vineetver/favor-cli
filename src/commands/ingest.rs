use std::path::PathBuf;

use serde_json::json;

use crate::cli::GenomeBuild;
use crate::commands::{self, IngestConfig};
use crate::config::Config;
use crate::store::list::{VariantSetKind, VariantSetWriter};
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::ingest::format::FormatRegistry;
use crate::ingest::{self, BuildGuess, InputFormat};
use crate::output::{bail_if_cancelled, Output};
use crate::resource::Resources;

/// Resolve inputs: expand directories to their VCF/parquet contents, validate all exist.
fn resolve_inputs(raw: Vec<PathBuf>) -> Result<Vec<PathBuf>, CohortError> {
    let mut resolved = Vec::new();
    for p in raw {
        if !p.exists() {
            return Err(CohortError::Input(format!(
                "File not found: {}",
                p.display()
            )));
        }
        if p.is_dir() {
            let mut found = Vec::new();
            for entry in std::fs::read_dir(&p)
                .map_err(|e| CohortError::Input(format!("Cannot read dir '{}': {e}", p.display())))?
            {
                let entry =
                    entry.map_err(|e| CohortError::Input(format!("Dir entry error: {e}")))?;
                let name = entry.file_name().to_string_lossy().to_lowercase();
                if name.ends_with(".vcf.gz")
                    || name.ends_with(".vcf.bgz")
                    || name.ends_with(".vcf")
                    || name.ends_with(".bcf")
                    || name.ends_with(".tsv")
                    || name.ends_with(".tsv.gz")
                    || name.ends_with(".csv")
                    || name.ends_with(".csv.gz")
                    || name.ends_with(".parquet")
                {
                    found.push(entry.path());
                }
            }
            if found.is_empty() {
                return Err(CohortError::Input(format!(
                    "No supported files in directory '{}'. Expected .vcf.gz, .tsv, .csv, or .parquet.",
                    p.display()
                )));
            }
            found.sort();
            resolved.extend(found);
        } else {
            resolved.push(p);
        }
    }
    if resolved.is_empty() {
        return Err(CohortError::Input("No input files specified.".into()));
    }
    Ok(resolved)
}

/// Entry point from CLI dispatch. Validates raw args, builds typed config, delegates to `run_ingest`.
pub fn handle(
    raw_inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    emit_sql: bool,
    build_override: Option<GenomeBuild>,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    let inputs = resolve_inputs(raw_inputs)?;

    let first = &inputs[0];
    let output_path = output.unwrap_or_else(|| {
        let stem = first.file_stem().unwrap_or_default().to_string_lossy();
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
    });

    let config = IngestConfig {
        inputs,
        output: output_path,
        emit_sql,
        build_override,
    };

    if dry_run {
        return run_ingest_dry(&config, out);
    }

    run_ingest(&config, out)
}

/// Backward-compatible entry point -- delegates to `handle`.
pub fn run(
    raw_inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    emit_sql: bool,
    build_override: Option<GenomeBuild>,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    handle(raw_inputs, output, emit_sql, build_override, out, dry_run)
}

fn run_ingest_dry(config: &IngestConfig, out: &dyn Output) -> Result<(), CohortError> {
    let first = &config.inputs[0];
    let mut analysis = ingest::analyze(first)?;
    validate_all_same_format(&config.inputs, &analysis)?;
    apply_build_override(&mut analysis, config, first)?;
    report_analysis(&analysis, &config.inputs, out);

    let file_list: Vec<String> = config
        .inputs
        .iter()
        .map(|p| p.to_string_lossy().into())
        .collect();
    let plan = commands::DryRunPlan {
        command: "ingest".into(),
        inputs: json!({
            "files": file_list,
            "file_count": config.inputs.len(),
            "format": format!("{:?}", analysis.format),
            "join_key": format!("{:?}", analysis.join_key),
            "build": format!("{:?}", analysis.build_guess),
            "columns_mapped": analysis.columns.len(),
            "columns_ambiguous": analysis.ambiguous.len(),
            "needs_intervention": analysis.needs_intervention(),
        }),
        memory: commands::MemoryEstimate::default_estimate(),
        output_path: config.output.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
}

/// Core ingest pipeline: detect → validate → ingest → report.
pub fn run_ingest(config: &IngestConfig, out: &dyn Output) -> Result<(), CohortError> {
    let first = &config.inputs[0];
    let mut analysis = ingest::analyze(first)?;
    validate_all_same_format(&config.inputs, &analysis)?;
    apply_build_override(&mut analysis, config, first)?;
    report_analysis(&analysis, &config.inputs, out);

    if analysis.format == InputFormat::Vcf {
        return ingest_vcf(config, out);
    }
    if config.emit_sql || analysis.needs_intervention() {
        return emit_sql_script(config, &analysis, out);
    }
    ingest_tabular(config, &analysis, out)
}

fn validate_all_same_format(
    inputs: &[PathBuf],
    analysis: &ingest::Analysis,
) -> Result<(), CohortError> {
    if inputs.len() <= 1 {
        return Ok(());
    }
    let registry = FormatRegistry::new();
    let first = &inputs[0];
    for p in &inputs[1..] {
        let detected = registry.detect(p)?;
        if detected.format != analysis.format {
            return Err(CohortError::Input(format!(
                "Mixed formats: {} is {:?} but {} is {} ({:?}). All inputs must share the same format.",
                first.display(), analysis.format, p.display(), detected.handler_name, detected.format,
            )));
        }
    }
    Ok(())
}

fn ingest_vcf(config: &IngestConfig, out: &dyn Output) -> Result<(), CohortError> {
    let first = &config.inputs[0];

    if config.emit_sql {
        let analysis = ingest::analyze(first)?;
        let script = ingest::sql::generate_script(&analysis, first, &config.output);
        let script_path = config.output.with_extension("ingest.sql");
        std::fs::write(&script_path, &script).map_err(|e| {
            CohortError::Resource(format!("Cannot write '{}': {e}", script_path.display()))
        })?;
        out.status(&format!("VCF hint script: {}", script_path.display()));
        return Ok(());
    }

    let resources = Resources::detect_configured();
    out.status(&format!(
        "Ingesting {} VCF file(s) ({}, {} threads)",
        config.inputs.len(),
        resources.memory_human(),
        resources.threads
    ));

    let source_str = if config.inputs.len() == 1 {
        first.display().to_string()
    } else {
        format!("{} files", config.inputs.len())
    };
    let mut vs_writer =
        VariantSetWriter::new(&config.output, ingest::JoinKey::ChromPosRefAlt, &source_str)?;
    vs_writer.set_kind(VariantSetKind::Ingested);

    bail_if_cancelled(out)?;
    let result = ingest::vcf::ingest_vcfs(
        &config.inputs,
        &mut vs_writer,
        resources.memory_bytes,
        resources.threads,
        out,
    )?;
    bail_if_cancelled(out)?;
    let vs = vs_writer.finish()?;

    out.success(&format!(
        "Ingested {} variants from {} file(s) -> {}",
        result.variant_count,
        config.inputs.len(),
        vs.root().display()
    ));
    if result.filtered_contigs > 0 {
        out.status(&format!(
            "  {} records on non-standard contigs filtered",
            result.filtered_contigs
        ));
    }
    if result.multiallelic_split > 0 {
        out.status(&format!(
            "  {} multi-allelic sites split to biallelic",
            result.multiallelic_split
        ));
    }

    out.result_json(&json!({
        "status": "ok",
        "output": vs.root().to_string_lossy(),
        "variant_count": result.variant_count,
        "file_count": config.inputs.len(),
        "join_key": "chrom_pos_ref_alt",
    }));

    Ok(())
}

fn emit_sql_script(
    config: &IngestConfig,
    analysis: &ingest::Analysis,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let first = &config.inputs[0];
    let script = ingest::sql::generate_script(analysis, first, &config.output);
    let script_path = config.output.with_extension("ingest.sql");
    std::fs::write(&script_path, &script).map_err(|e| {
        CohortError::Resource(format!("Cannot write '{}': {e}", script_path.display()))
    })?;

    if analysis.needs_intervention() {
        out.warn("Cannot auto-ingest -- ambiguities detected. Edit the script, then re-run:");
        out.warn(&format!("  cohort ingest {}", first.display()));
    } else {
        out.status(&format!("SQL script written to: {}", script_path.display()));
    }

    out.result_json(&json!({
        "status": analysis.status(),
        "format": analysis.format,
        "join_key": analysis.join_key,
        "build": analysis.build_guess,
        "coord_base": analysis.coord_base,
        "columns_mapped": analysis.columns.len(),
        "columns_ambiguous": analysis.ambiguous.len(),
        "columns_unmapped": analysis.unmapped.len(),
        "script_path": script_path.to_string_lossy(),
    }));

    Ok(())
}

fn ingest_tabular(
    config: &IngestConfig,
    analysis: &ingest::Analysis,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let first = &config.inputs[0];
    out.status("Ingesting...");

    let resources = Resources::detect_configured();
    let engine = DfEngine::new(&resources)?;

    bail_if_cancelled(out)?;
    register_tabular_inputs(&engine, &config.inputs, analysis)?;

    let select_sql = ingest::sql::generate_select(analysis);
    let copy_sql = ingest::sql::copy_statement(&select_sql, &config.output);
    engine.execute(&copy_sql)?;

    let mut vs_writer = VariantSetWriter::new(
        &config.output,
        analysis.join_key,
        &first.display().to_string(),
    )?;
    vs_writer.set_kind(VariantSetKind::Ingested);
    vs_writer.scan_and_register()?;
    let vs = vs_writer.finish()?;

    out.success(&format!("Ingested -> {}", vs.root().display()));
    out.status(&format!("  Join key: {:?}", analysis.join_key));

    out.result_json(&json!({
        "status": "ok",
        "output": vs.root().to_string_lossy(),
        "join_key": analysis.join_key,
        "build": analysis.build_guess,
    }));

    Ok(())
}

fn register_tabular_inputs(
    engine: &DfEngine,
    inputs: &[PathBuf],
    analysis: &ingest::Analysis,
) -> Result<(), CohortError> {
    let delimiter = analysis.delimiter.unwrap_or(ingest::Delimiter::Tab);
    match analysis.format {
        InputFormat::Tabular => {
            if inputs.len() == 1 {
                engine.register_csv("_ingest_input", &inputs[0], delimiter.byte())?;
            } else {
                for (i, p) in inputs.iter().enumerate() {
                    engine.register_csv(&format!("_ingest_part_{i}"), p, delimiter.byte())?;
                }
                let union_parts: Vec<String> = (0..inputs.len())
                    .map(|i| format!("SELECT * FROM _ingest_part_{i}"))
                    .collect();
                engine.execute(&format!(
                    "CREATE VIEW _ingest_input AS {}",
                    union_parts.join(" UNION ALL "),
                ))?;
            }
        }
        InputFormat::Parquet => {
            if inputs.len() == 1 {
                engine.register_parquet_dir("_ingest_input", &inputs[0])?;
            } else {
                for (i, p) in inputs.iter().enumerate() {
                    engine.register_parquet_file(&format!("_ingest_part_{i}"), p)?;
                }
                let union_parts: Vec<String> = (0..inputs.len())
                    .map(|i| format!("SELECT * FROM _ingest_part_{i}"))
                    .collect();
                engine.execute(&format!(
                    "CREATE VIEW _ingest_input AS {}",
                    union_parts.join(" UNION ALL "),
                ))?;
            }
        }
        InputFormat::Vcf => unreachable!(),
    }
    Ok(())
}

fn apply_build_override(
    analysis: &mut ingest::Analysis,
    config: &IngestConfig,
    first: &std::path::Path,
) -> Result<(), CohortError> {
    match config.build_override {
        Some(GenomeBuild::Hg19) => {
            analysis.build_guess = BuildGuess::Hg19 {
                match_rate_hg38: 0.0,
                match_rate_hg19: 1.0,
            };
        }
        Some(GenomeBuild::Hg38) => {
            analysis.build_guess = BuildGuess::Hg38;
        }
        None if analysis.format != InputFormat::Vcf => {
            if let Ok(global) = Config::load() {
                if !global.data.root_dir.is_empty() {
                    let resources = Resources::detect_with_config(&global.resources);
                    if let Ok(engine) = DfEngine::new(&resources) {
                        let _ = ingest::detect::detect_build_and_coords(
                            analysis, first, &engine, &global,
                        );
                    }
                }
            }
        }
        None => {}
    }
    Ok(())
}

fn report_analysis(analysis: &ingest::Analysis, inputs: &[PathBuf], out: &dyn Output) {
    let first = &inputs[0];
    if inputs.len() == 1 {
        out.status(&format!("Analyzing {}...", first.display()));
    } else {
        out.status(&format!(
            "Analyzing {} files (first: {})...",
            inputs.len(),
            first.file_name().unwrap_or_default().to_string_lossy()
        ));
    }

    out.status(&format!("  Format: {:?}", analysis.format));
    out.status(&format!("  Join key: {:?}", analysis.join_key));

    for mapping in &analysis.columns {
        if mapping.canonical != mapping.input_name.to_lowercase() {
            out.status(&format!(
                "  {} -> {}",
                mapping.input_name, mapping.canonical
            ));
        }
    }
    for amb in &analysis.ambiguous {
        out.warn(&format!("  {} -- {}", amb.column, amb.reason));
    }

    match &analysis.build_guess {
        BuildGuess::Hg38 => out.status("  Build: hg38"),
        BuildGuess::Hg19 { .. } => out.warn("  Build: likely hg19 -- needs liftover"),
        BuildGuess::Unknown => out.status("  Build: unknown (no annotations to probe)"),
    }

    match analysis.coord_base {
        ingest::CoordBase::OneBased => out.status("  Coordinates: 1-based"),
        ingest::CoordBase::ZeroBased => {
            out.warn("  Coordinates: 0-based (will convert to 1-based)")
        }
        ingest::CoordBase::Unknown => {}
    }
}
