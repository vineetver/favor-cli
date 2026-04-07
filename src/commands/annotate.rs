use std::path::PathBuf;

use serde_json::json;

use crate::commands::{self, AnnotateConfig};
use crate::config::{Config, Tier};
use crate::data::{parquet_column_names, parquet_row_count, AnnotationDb};
use crate::data::{VariantSet, VariantSetKind, VariantSetWriter};
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::ingest::{ColumnContract, ColumnRequirement, JoinKey};
use crate::output::Output;
use crate::resource::Resources;

const ANNOTATE_JOIN_COLUMNS: &[ColumnRequirement] = &[
    ColumnRequirement {
        name: "chromosome",
        source: "ingested variant set",
        used_by: "annotation join",
    },
    ColumnRequirement {
        name: "position",
        source: "ingested variant set",
        used_by: "annotation join",
    },
    ColumnRequirement {
        name: "ref",
        source: "ingested variant set",
        used_by: "annotation join (allele match)",
    },
    ColumnRequirement {
        name: "alt",
        source: "ingested variant set",
        used_by: "annotation join (allele match)",
    },
];

const JOIN_KEY_COLUMNS: &[&str] = &["chromosome", "position", "ref", "alt"];

/// Entry point from CLI dispatch. Validates raw args, builds typed config, delegates to `run_annotate`.
pub fn handle(
    input: PathBuf,
    output_path: Option<PathBuf>,
    full: bool,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    if !input.exists() {
        return Err(CohortError::Input(format!(
            "Variant set not found: '{}'. Run `cohort ingest <file>` first to produce one.",
            input.display()
        )));
    }

    let global_config = Config::load_configured()?;
    let tier = if full {
        Tier::Full
    } else {
        global_config.data.tier
    };

    let output = output_path.unwrap_or_else(|| {
        let name = input.file_name().unwrap_or_default().to_string_lossy();
        let stem = name
            .strip_suffix(".ingested")
            .or_else(|| name.strip_suffix("/"))
            .unwrap_or(&name);
        input
            .parent()
            .unwrap_or(&input)
            .join(format!("{stem}.annotated"))
    });

    let config = AnnotateConfig {
        input,
        output,
        tier,
        data_root: global_config.root_dir(),
    };

    if dry_run {
        return emit_dry_run(&config, out);
    }

    run_annotate(&config, out)
}

/// Backward-compatible entry point — delegates to `handle`.
pub fn run(
    input: PathBuf,
    output_path: Option<PathBuf>,
    full: bool,
    out: &dyn Output,
    dry_run: bool,
) -> Result<(), CohortError> {
    handle(input, output_path, full, out, dry_run)
}

fn emit_dry_run(config: &AnnotateConfig, out: &dyn Output) -> Result<(), CohortError> {
    let input_vs = VariantSet::open(&config.input)?;
    let plan = commands::DryRunPlan {
        command: "annotate".into(),
        inputs: json!({
            "file": config.input.to_string_lossy(),
            "variant_count": input_vs.count(),
            "join_key": format!("{:?}", input_vs.join_key()),
            "tier": config.tier.as_str(),
        }),
        memory: commands::MemoryEstimate::default_estimate(),
        output_path: config.output.to_string_lossy().into(),
    };
    commands::emit(&plan, out);
    Ok(())
}

/// Core annotate pipeline: open → validate → join → report.
pub fn run_annotate(config: &AnnotateConfig, out: &dyn Output) -> Result<(), CohortError> {
    let input = VariantSet::open(&config.input)?;
    validate_input(&input)?;
    let db = AnnotationDb::open_tier(config.tier, &config.data_root)?;
    let resources = Resources::detect_configured();
    let engine = DfEngine::new(&resources)?;

    out.status(&format!(
        "Input: {} variants, join key: {:?}",
        input.count(),
        input.join_key()
    ));
    out.status(&format!(
        "Resources: {}, {} threads ({})",
        resources.memory_human(),
        resources.threads,
        resources.environment()
    ));
    out.status(&format!(
        "Annotating against cohort-{} ({})...",
        config.tier,
        config.tier.size_human()
    ));

    let result = annotate_join(&input, &db, &engine, &config.output, out)?;
    report_result(&result, config, &engine, &input, &db, out)
}

fn validate_input(input: &VariantSet) -> Result<(), CohortError> {
    let contract = ColumnContract {
        command: "annotate",
        required: ANNOTATE_JOIN_COLUMNS,
    };
    let missing = contract.check(input.columns());
    if !missing.is_empty() {
        return Err(CohortError::Input(format!(
            "Missing required columns for annotation join:\n{}\n\
             Re-ingest your data: cohort ingest <file>",
            ColumnContract::format_missing(&missing)
        )));
    }
    if input.count() == 0 {
        return Err(CohortError::Input(format!(
            "Input '{}' has 0 variants. Check that `cohort ingest` completed successfully.",
            input.root().display()
        )));
    }
    Ok(())
}

struct AnnotateResult {
    output: VariantSet,
    input_count: i64,
    output_count: i64,
}

fn annotate_join(
    input: &VariantSet,
    db: &AnnotationDb,
    engine: &DfEngine,
    output_path: &std::path::Path,
    out: &dyn Output,
) -> Result<AnnotateResult, CohortError> {
    let join_key = input.join_key();
    let input_count = input.count() as i64;

    let input_struct_cols: Vec<&str> = input
        .columns()
        .iter()
        .map(|s| s.as_str())
        .filter(|c| !JOIN_KEY_COLUMNS.contains(c))
        .collect();
    let input_struct_expr = build_struct_expr(&input_struct_cols);

    let mut writer =
        VariantSetWriter::new(output_path, join_key, &input.root().display().to_string())?;
    writer.set_kind(VariantSetKind::Annotated { tier: db.tier() });

    match join_key {
        JoinKey::ChromPosRefAlt | JoinKey::ChromPos => {
            join_per_chromosome(input, db, engine, &input_struct_expr, &mut writer, out)?;
        }
        JoinKey::Rsid => {
            join_via_rsid(
                input,
                db,
                engine,
                &input_struct_expr,
                output_path,
                &mut writer,
                out,
            )?;
        }
    }

    let output = writer.finish()?;
    let output_count = output.count() as i64;
    Ok(AnnotateResult {
        output,
        input_count,
        output_count,
    })
}

fn build_struct_expr(non_join_cols: &[&str]) -> String {
    if non_join_cols.is_empty() {
        return String::new();
    }
    let args: Vec<String> = non_join_cols
        .iter()
        .map(|c| format!("'{c}', u.\"{c}\""))
        .collect();
    format!("named_struct({})", args.join(", "))
}

fn join_per_chromosome(
    input: &VariantSet,
    db: &AnnotationDb,
    engine: &DfEngine,
    input_struct_expr: &str,
    writer: &mut VariantSetWriter,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let exact_allele = input.join_key() == JoinKey::ChromPosRefAlt;
    let chroms = input.chromosomes();
    let total = chroms.len();

    let input_select = if input_struct_expr.is_empty() {
        String::new()
    } else {
        format!("{input_struct_expr} AS input, ")
    };

    let allele_condition = if exact_allele {
        "AND u.\"ref\" = a.ref_vcf AND u.alt = a.alt_vcf"
    } else {
        ""
    };

    let mut columns_set = false;
    for (i, chrom) in chroms.iter().enumerate() {
        let chrom_parquet = match db.chrom_parquet(chrom) {
            Some(p) => p,
            None => {
                out.warn(&format!("  chr{chrom}: skipping (no annotation file)"));
                continue;
            }
        };

        out.status(&format!("  chr{chrom} ({}/{})", i + 1, total));

        let input_parquet = input.chrom_dir(chrom).join("data.parquet");
        engine.register_parquet_file("_input", &input_parquet)?;
        engine.register_parquet_file("_ann", &chrom_parquet)?;

        let out_path = writer.chrom_path(chrom)?;

        let join_sql = format!(
            "SELECT {input_select}a.* \
             FROM _input u \
             INNER JOIN _ann a \
                 ON u.position = a.position {allele_condition}",
        );
        if i == 0 {
            let plan = engine.explain(&join_sql)?;
            out.status(&format!("  Plan: {plan}"));
        }

        engine.execute(&format!(
            "COPY ({join_sql}) TO '{out_path}' \
             STORED AS PARQUET OPTIONS (compression 'zstd(4)')",
            out_path = out_path.display(),
        ))?;

        let count = if out_path.exists() {
            parquet_row_count(&out_path)?
        } else {
            0
        };
        let size = std::fs::metadata(&out_path).map_or(0, |m| m.len());
        if count > 0 {
            writer.register_chrom(chrom, count, size);
            if !columns_set {
                writer.set_columns(parquet_column_names(&out_path)?);
                columns_set = true;
            }
        } else {
            let _ = std::fs::remove_file(&out_path);
        }

        let _ = engine.execute("DROP TABLE IF EXISTS _input");
        let _ = engine.execute("DROP TABLE IF EXISTS _ann");
    }
    Ok(())
}

fn join_via_rsid(
    input: &VariantSet,
    db: &AnnotationDb,
    engine: &DfEngine,
    input_struct_expr: &str,
    output_path: &std::path::Path,
    writer: &mut VariantSetWriter,
    out: &dyn Output,
) -> Result<(), CohortError> {
    out.status("  Resolving rsids across all chromosomes...");

    engine.register_parquet_dir("_user", input.root())?;

    let lookup_dir = db
        .root()
        .parent()
        .unwrap_or(db.root())
        .join("lookup/rsid_lookup");
    engine.register_parquet_dir("_rsid_lookup", &lookup_dir)?;
    engine.register_parquet_dir("_annotations", db.root())?;

    let input_select = if input_struct_expr.is_empty() {
        String::new()
    } else {
        format!("{input_struct_expr} AS input, ")
    };

    engine.execute(
        "CREATE TABLE _rsid_resolved AS \
         SELECT u.*, lk.chromosome AS _chr, lk.position AS _pos, \
                lk.ref_vcf AS _ref, lk.alt_vcf AS _alt \
         FROM _user u \
         INNER JOIN _rsid_lookup lk ON u.rsid = lk.rsid",
    )?;

    engine.execute(&format!(
        "COPY (\
             SELECT {input_select}a.* \
             FROM _rsid_resolved r \
             INNER JOIN _annotations a \
                 ON r._chr = a.chromosome \
                 AND r._pos = a.position \
                 AND r._ref = a.ref_vcf \
                 AND r._alt = a.alt_vcf\
         ) TO '{output}/' \
         STORED AS PARQUET \
         PARTITIONED BY (chromosome) \
         OPTIONS (compression 'zstd(4)')",
        output = output_path.display(),
    ))?;

    writer.scan_and_register()?;
    Ok(())
}

fn report_result(
    result: &AnnotateResult,
    config: &AnnotateConfig,
    engine: &DfEngine,
    input: &VariantSet,
    db: &AnnotationDb,
    out: &dyn Output,
) -> Result<(), CohortError> {
    let match_rate = if result.input_count > 0 {
        result.output_count as f64 / result.input_count as f64
    } else {
        0.0
    };

    if result.output_count == 0 {
        let diagnostic = diagnose_zero_matches(engine, input, db.root())?;
        out.warn("0 variants annotated.");
        out.warn(&diagnostic.message);
        out.result_json(&json!({
            "status": "no_matches",
            "input_count": result.input_count,
            "output_count": 0,
            "diagnostic": diagnostic.message,
            "suggested_fix": diagnostic.suggested_fix,
        }));
        return Err(CohortError::Analysis(format!(
            "0 variants matched annotations. {}",
            diagnostic.message
        )));
    }

    if match_rate < 0.05 {
        let diagnostic = diagnose_zero_matches(engine, input, db.root())?;
        out.warn(&format!(
            "Low match rate: {}/{} ({:.1}%). {}",
            result.output_count,
            result.input_count,
            match_rate * 100.0,
            diagnostic.message,
        ));
        out.result_json(&json!({
            "status": "low_match_rate",
            "input_count": result.input_count,
            "output_count": result.output_count,
            "match_rate": match_rate,
            "diagnostic": diagnostic.message,
            "suggested_fix": diagnostic.suggested_fix,
            "output": result.output.root().to_string_lossy(),
        }));
        return Ok(());
    }

    let unmatched = result.input_count - result.output_count;
    out.success(&format!(
        "Annotated {}/{} variants ({:.1}%) -> {}",
        result.output_count,
        result.input_count,
        match_rate * 100.0,
        result.output.root().display(),
    ));
    if unmatched > 0 {
        out.status(&format!("  {unmatched} variants not found in annotations"));
    }

    out.result_json(&json!({
        "status": "ok",
        "output": result.output.root().to_string_lossy(),
        "input_count": result.input_count,
        "annotated_count": result.output_count,
        "match_rate": match_rate,
        "tier": config.tier.as_str(),
    }));

    Ok(())
}

struct Diagnostic {
    message: String,
    suggested_fix: String,
}

fn diagnose_zero_matches(
    engine: &DfEngine,
    input_vs: &VariantSet,
    annotations_dir: &std::path::Path,
) -> Result<Diagnostic, CohortError> {
    let chr1_path = annotations_dir.join("chromosome=1/sorted.parquet");
    if !chr1_path.exists() {
        return Ok(Diagnostic {
            message: "Annotations incomplete — chromosome=1 not found.".into(),
            suggested_fix: "cohort data pull".into(),
        });
    }

    engine.register_parquet_dir("_diag_input", input_vs.root())?;
    engine.register_parquet_file("_diag_ann", &chr1_path)?;

    let chroms = engine
        .query_strings("SELECT DISTINCT chromosome FROM _diag_input LIMIT 5")
        .unwrap_or_default();

    if chroms.iter().any(|c| c.starts_with("chr")) {
        return Ok(Diagnostic {
            message:
                "Chromosome names have 'chr' prefix — annotations use bare numbers. Re-ingest."
                    .into(),
            suggested_fix: "cohort ingest <original_file>".into(),
        });
    }

    let position_hits = engine
        .query_scalar(
            "SELECT COUNT(*) FROM (\
            SELECT position FROM _diag_input WHERE chromosome = '1' LIMIT 10\
        ) s \
        WHERE EXISTS (SELECT 1 FROM _diag_ann a WHERE a.position = s.position)",
        )
        .unwrap_or(0);

    if position_hits == 0 {
        return Ok(Diagnostic {
            message: "No positions match hg38 annotations. Input may be hg19 or different build."
                .into(),
            suggested_fix: "cohort ingest <original_file> --build hg19 --emit-sql".into(),
        });
    }

    let mismatch_samples = engine
        .collect(
            "SELECT u.position, u.\"ref\" AS in_ref, u.alt AS in_alt, \
                a.ref_vcf AS ann_ref, a.alt_vcf AS ann_alt \
         FROM _diag_input u \
         INNER JOIN _diag_ann a ON u.position = a.position \
         WHERE u.chromosome = '1' \
           AND (u.\"ref\" != a.ref_vcf OR u.alt != a.alt_vcf) \
         LIMIT 5",
        )
        .unwrap_or_default();

    let mut examples = Vec::new();
    for batch in &mismatch_samples {
        for row in 0..batch.num_rows() {
            let pos = arrow::util::display::array_value_to_string(batch.column(0), row)
                .unwrap_or_default();
            let in_ref = arrow::util::display::array_value_to_string(batch.column(1), row)
                .unwrap_or_default();
            let in_alt = arrow::util::display::array_value_to_string(batch.column(2), row)
                .unwrap_or_default();
            let ann_ref = arrow::util::display::array_value_to_string(batch.column(3), row)
                .unwrap_or_default();
            let ann_alt = arrow::util::display::array_value_to_string(batch.column(4), row)
                .unwrap_or_default();
            examples.push(format!(
                "  pos {pos}: input {in_ref}/{in_alt}, annotation {ann_ref}/{ann_alt}"
            ));
        }
    }

    let detail = if examples.is_empty() {
        String::new()
    } else {
        format!("\nMismatch samples:\n{}", examples.join("\n"))
    };

    Ok(Diagnostic {
        message: format!(
            "Positions match but alleles differ — likely a normalization mismatch. \
             Re-ingest to apply parsimony normalization.{detail}"
        ),
        suggested_fix: "cohort ingest <original_file>".into(),
    })
}
