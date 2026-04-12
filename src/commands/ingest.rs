use std::path::PathBuf;

use serde_json::json;

use crate::cli::GenomeBuild;
use crate::commands::{self, IngestConfig};
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::ingest::format::FormatRegistry;
use crate::ingest::{self, BuildGuess, InputFormat};
use crate::output::{bail_if_cancelled, Output};
use crate::runtime::Engine;
use crate::staar::genotype::read_sample_names;
use crate::store::cohort::{BuildOpts, CohortId, CohortSources, ProbeReason, StoreProbe};
use crate::store::list::{VariantSet, VariantSetKind, VariantSetWriter};

pub(crate) fn resolve_inputs(raw: Vec<PathBuf>) -> Result<Vec<PathBuf>, CohortError> {
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

pub fn build_config(
    raw_inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    emit_sql: bool,
    build_override: Option<GenomeBuild>,
    annotations: Option<PathBuf>,
    cohort_id: Option<String>,
    rebuild: bool,
) -> Result<IngestConfig, CohortError> {
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

    Ok(IngestConfig {
        inputs,
        output: output_path,
        emit_sql,
        build_override,
        annotations,
        cohort_id,
        rebuild,
    })
}

pub fn run_ingest(
    engine: &Engine,
    config: &IngestConfig,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    let first = &config.inputs[0];
    let mut analysis = ingest::analyze(first)?;
    validate_all_same_format(&config.inputs, &analysis)?;
    apply_build_override(engine, &mut analysis, config, first)?;
    report_analysis(&analysis, &config.inputs, out);

    if dry_run {
        return run_ingest_dry(config, &analysis, out);
    }

    if analysis.format == InputFormat::Vcf {
        let n_samples = read_sample_names(first)?.len();
        if config.annotations.is_some() && n_samples > 0 {
            return run_cohort_build(engine, config, n_samples, out);
        }
        return ingest_vcf(engine, config, n_samples, out);
    }
    if config.emit_sql || analysis.needs_intervention() {
        return emit_sql_script(config, &analysis, out);
    }
    ingest_tabular(engine, config, &analysis, out)
}

fn run_cohort_build(
    engine: &Engine,
    config: &IngestConfig,
    n_samples: usize,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let vcf_path = &config.inputs[0];
    let annotations_path = config.annotations.as_ref().ok_or_else(|| {
        CohortError::Input(format!(
            "Multi-sample VCF detected ({n_samples} samples in {}). \
             `favor ingest` needs an annotated variant set for the genotype \
             store build:\n\
             \n\
             1. favor ingest <vcf> --output variants.set    (sites only — drops genotypes)\n\
             2. favor annotate variants.set --full\n\
             3. favor ingest <vcf> --annotations variants.set.annotated --cohort-id <id>",
            vcf_path.display()
        ))
    })?;
    if !annotations_path.exists() {
        return Err(CohortError::Input(format!(
            "Annotations not found at '{}'",
            annotations_path.display()
        )));
    }
    let ann_vs = VariantSet::open(annotations_path)?;
    ann_vs.require_annotated()?;
    drop(ann_vs);

    let cohort_id = CohortId::new(
        config
            .cohort_id
            .clone()
            .unwrap_or_else(|| commands::derive_cohort_id(vcf_path)),
    );

    let cohort = engine.cohort(&cohort_id);
    let probe = if config.rebuild {
        if let Err(e) = engine.store().cache().prune_cohort(&cohort_id) {
            out.warn(&format!("  prune cache before rebuild: {e}"));
        }
        StoreProbe {
            store_dir: cohort.dir().to_path_buf(),
            manifest: None,
            miss_reason: Some(ProbeReason::NoManifest),
        }
    } else {
        cohort.probe(&config.inputs, annotations_path)
    };

    let staging_dir = {
        let parent = cohort.dir().parent().ok_or_else(|| {
            CohortError::Resource(format!(
                "cohort dir '{}' has no parent",
                cohort.dir().display()
            ))
        })?;
        let mut name = cohort
            .dir()
            .file_name()
            .ok_or_else(|| {
                CohortError::Resource(format!(
                    "cohort dir '{}' has no file name",
                    cohort.dir().display()
                ))
            })?
            .to_os_string();
        name.push(".geno_staging");
        parent.join(name)
    };

    out.status(&format!(
        "Building cohort '{}' from multi-sample VCF ({n_samples} samples)",
        cohort_id.as_str()
    ));

    let result = cohort.build_or_load(
        CohortSources {
            genotypes: &config.inputs,
            annotations: annotations_path,
        },
        BuildOpts {
            staging_dir: &staging_dir,
            rebuild: config.rebuild,
            probe,
        },
        engine,
        out,
    )?;

    out.success(&format!(
        "Cohort '{}': {} variants × {} samples at {}",
        cohort_id.as_str(),
        result.manifest.n_variants,
        result.manifest.n_samples,
        result.store_dir.display(),
    ));

    out.result_json(&json!({
        "status": "ok",
        "cohort_id": cohort_id.as_str(),
        "store_dir": result.store_dir.to_string_lossy(),
        "n_variants": result.manifest.n_variants,
        "n_samples": result.manifest.n_samples,
    }));

    Ok(())
}

fn run_ingest_dry(
    config: &IngestConfig,
    analysis: &ingest::Analysis,
    out: &dyn Output,
) -> Result<(), CohortError> {
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
        runtime: None,
        output_path: config.output.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
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

fn ingest_vcf(
    engine: &Engine,
    config: &IngestConfig,
    n_samples: usize,
    out: &dyn Output,
) -> Result<(), CohortError> {
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

    let resources = engine.resources();
    let dual_write = n_samples > 0;

    if dual_write {
        out.status(&format!(
            "Ingesting {} VCF file(s) + extracting genotypes ({} samples, {}, {} threads)",
            config.inputs.len(),
            n_samples,
            resources.memory_human(),
            resources.threads
        ));
    } else {
        out.status(&format!(
            "Ingesting {} VCF file(s) ({}, {} threads)",
            config.inputs.len(),
            resources.memory_human(),
            resources.threads
        ));
    }

    let geno_dir = config.output.with_extension("genotypes");

    bail_if_cancelled(out)?;

    let (result, vs) = if config.inputs.len() > 1 {
        // Parallel path: each file gets its own worker.
        let result = ingest::vcf::ingest_vcfs_parallel(
            &config.inputs,
            &config.output,
            if dual_write { Some(geno_dir.as_path()) } else { None },
            n_samples,
            resources.memory_bytes,
            resources.threads,
            out,
        )?;

        let source_str = format!("{} files", config.inputs.len());
        let mut vs_writer =
            VariantSetWriter::new(&config.output, ingest::JoinKey::ChromPosRefAlt, &source_str)?;
        vs_writer.set_kind(VariantSetKind::Ingested);
        vs_writer.scan_and_register()?;
        let vs = vs_writer.finish()?;

        if dual_write {
            let sample_names = read_sample_names(first)?;
            let source_vcfs: Vec<String> = config.inputs.iter().map(|p| p.display().to_string()).collect();
            let meta = crate::staar::genotype::GenotypeMeta {
                version: 1,
                n_samples,
                chromosomes: vs.chromosomes().iter().map(|c| c.to_string()).collect(),
                source_vcfs,
            };
            let meta_path = geno_dir.join("genotypes.json");
            std::fs::write(&meta_path, serde_json::to_string_pretty(&meta)
                .map_err(|e| CohortError::Resource(format!("JSON: {e}")))?
            ).map_err(|e| CohortError::Resource(format!("Cannot write '{}': {e}", meta_path.display())))?;
            let sidecar = geno_dir.join("samples.txt");
            std::fs::write(&sidecar, sample_names.join("\n"))
                .map_err(|e| CohortError::Resource(format!("Cannot write '{}': {e}", sidecar.display())))?;
            out.success(&format!(
                "  Genotypes: {} variants x {} samples -> {}",
                result.genotype_variants, n_samples, geno_dir.display()
            ));
        }

        (result, vs)
    } else {
        // Sequential path: single file, zero overhead.
        let source_str = first.display().to_string();
        let mut vs_writer =
            VariantSetWriter::new(&config.output, ingest::JoinKey::ChromPosRefAlt, &source_str)?;
        vs_writer.set_kind(VariantSetKind::Ingested);

        let mut geno_writer = if dual_write {
            Some(crate::staar::genotype::GenotypeWriter::new(
                n_samples, &geno_dir, resources.memory_bytes,
            )?)
        } else {
            None
        };

        let result = ingest::vcf::ingest_vcfs(
            &config.inputs,
            &mut vs_writer,
            geno_writer.as_mut(),
            resources.memory_bytes,
            resources.threads,
            out,
        )?;
        let vs = vs_writer.finish()?;

        if let Some(gw) = geno_writer {
            let sample_names = read_sample_names(first)?;
            let source_vcfs = config.inputs.iter().map(|p| p.display().to_string()).collect();
            let geno_result = gw.finish(&sample_names, source_vcfs)?;
            out.success(&format!(
                "  Genotypes: {} variants x {} samples -> {}",
                result.genotype_variants, n_samples, geno_result.output_dir.display()
            ));
        }

        (result, vs)
    };

    bail_if_cancelled(out)?;

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

    let mut result_json = json!({
        "status": "ok",
        "output": vs.root().to_string_lossy(),
        "variant_count": result.variant_count,
        "file_count": config.inputs.len(),
        "join_key": "chrom_pos_ref_alt",
    });
    if dual_write {
        result_json["genotype_dir"] = json!(geno_dir.to_string_lossy());
        result_json["genotype_variants"] = json!(result.genotype_variants);
        result_json["n_samples"] = json!(n_samples);
    }
    out.result_json(&result_json);

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
    engine: &Engine,
    config: &IngestConfig,
    analysis: &ingest::Analysis,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let first = &config.inputs[0];
    out.status("Ingesting...");

    let df = engine.df();

    bail_if_cancelled(out)?;
    register_tabular_inputs(df, &config.inputs, analysis)?;

    let select_sql = ingest::sql::generate_select(analysis);
    let copy_sql = ingest::sql::copy_statement(&select_sql, &config.output);
    df.execute(&copy_sql)?;

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
    engine: &Engine,
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
            // Pre-setup ingest is a first-class path: skip build detection when
            // no configured engine is available.
            if engine.config_opt().is_some() {
                let _ = ingest::detect::detect_build_and_coords(engine, analysis, first);
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
