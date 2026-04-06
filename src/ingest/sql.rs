//! SQL generator for ingest scripts.
//!
//! Takes an Analysis and produces a DataFusion-compatible SQL SELECT.
//! The programmatic path wraps this in a COPY TO statement.
//! The --emit-sql path writes it as a standalone script for inspection/editing.

use std::fmt::Write;
use std::path::Path;

use super::canonical_type;
use super::{Analysis, BuildGuess, CoordBase, InputFormat};

/// Generate the SELECT transformation SQL from an Analysis.
/// Does NOT include COPY TO — caller wraps this for execution or emits as-is.
pub fn generate_select(analysis: &Analysis) -> String {
    let mut sql = String::with_capacity(2048);

    match analysis.format {
        InputFormat::Tabular => write_tabular_select(&mut sql, analysis),
        InputFormat::Parquet => write_parquet_select(&mut sql, analysis),
        InputFormat::Vcf => write_vcf_stub(&mut sql),
    }

    sql
}

/// Generate a complete DataFusion SQL script for --emit-sql.
/// Includes comments, the SELECT, and COPY TO instruction.
pub fn generate_script(analysis: &Analysis, input_path: &Path, output_path: &Path) -> String {
    let mut script = String::with_capacity(4096);

    write_header(&mut script, analysis, input_path, output_path);

    let select = generate_select(analysis);
    let _ = writeln!(script, "-- Register the input file first:");
    match analysis.format {
        InputFormat::Tabular => {
            let delim = analysis.delimiter.unwrap_or(super::Delimiter::Tab);
            let _ = writeln!(
                script,
                "-- (programmatic: engine.register_csv(\"_ingest_input\", path, {:?}))",
                delim.char()
            );
        }
        InputFormat::Parquet => {
            let _ = writeln!(
                script,
                "-- (programmatic: engine.register_parquet(\"_ingest_input\", path))"
            );
        }
        InputFormat::Vcf => {}
    }
    let _ = writeln!(script);
    let _ = writeln!(script, "-- Transformation query:");
    let _ = writeln!(script, "{select}");
    let _ = writeln!(script);

    if analysis.format != InputFormat::Vcf {
        let _ = writeln!(script, "-- Execute as:");
        let _ = writeln!(script, "-- COPY (<above query>)");
        let _ = writeln!(script, "--   TO '{}/'", output_path.display());
        let _ = writeln!(script, "--   STORED AS PARQUET");
        let _ = writeln!(script, "--   PARTITIONED BY (chromosome)");
        let _ = writeln!(script, "--   OPTIONS (compression 'zstd(4)');");
    }

    if let BuildGuess::Hg19 {
        match_rate_hg38,
        match_rate_hg19,
    } = &analysis.build_guess
    {
        let _ = writeln!(script);
        let _ = writeln!(
            script,
            "-- WARNING: Input appears to be hg19 (hg38 match: {:.0}%, hg19 match: {:.0}%)",
            match_rate_hg38 * 100.0,
            match_rate_hg19 * 100.0
        );
        let _ = writeln!(
            script,
            "-- This data needs liftover to hg38 before annotation."
        );
    }

    script
}

/// Build the COPY TO statement for programmatic execution.
pub fn copy_statement(select_sql: &str, output_path: &Path) -> String {
    format!(
        "COPY ({select_sql}) TO '{output}/' \
         STORED AS PARQUET \
         PARTITIONED BY (chromosome) \
         OPTIONS (compression 'zstd(4)')",
        output = output_path.display(),
    )
}

fn write_header(sql: &mut String, analysis: &Analysis, input_path: &Path, output_path: &Path) {
    let _ = writeln!(
        sql,
        "-- FAVOR ingest: {} -> {}",
        input_path.display(),
        output_path.display()
    );
    let _ = writeln!(
        sql,
        "-- Format: {:?} | Join key: {:?}",
        analysis.format, analysis.join_key
    );

    if analysis.needs_intervention() {
        let _ = writeln!(
            sql,
            "-- STATUS: needs_edit — resolve FIXME comments before running"
        );
    } else {
        let _ = writeln!(sql, "-- STATUS: ok — ready to run as-is");
    }
    let _ = writeln!(sql);
}

fn write_tabular_select(sql: &mut String, analysis: &Analysis) {
    let _ = writeln!(sql, "SELECT");

    if let Some(chr_col) = &analysis.chr_col {
        let _ = writeln!(
            sql,
            "    CASE upper(regexp_replace(CAST(\"{chr_col}\" AS VARCHAR), '^chr', '', 'i'))"
        );
        let _ = writeln!(sql, "        WHEN 'M' THEN 'MT'");
        let _ = writeln!(sql, "        WHEN '23' THEN 'X'");
        let _ = writeln!(sql, "        WHEN '24' THEN 'Y'");
        let _ = writeln!(
            sql,
            "        ELSE upper(regexp_replace(CAST(\"{chr_col}\" AS VARCHAR), '^chr', '', 'i'))"
        );
        let _ = writeln!(sql, "    END AS chromosome,");
    } else {
        let _ = writeln!(sql, "    -- FIXME: no chromosome column detected");
    }

    if let Some(pos_col) = &analysis.pos_col {
        match analysis.coord_base {
            CoordBase::ZeroBased => {
                let _ = writeln!(sql, "    CAST(\"{pos_col}\" AS INT) + 1 AS position,");
            }
            _ => {
                let _ = writeln!(sql, "    CAST(\"{pos_col}\" AS INT) AS position,");
            }
        }
    } else {
        let _ = writeln!(sql, "    -- FIXME: no position column detected");
    }

    if let Some(ref_col) = &analysis.ref_col {
        let _ = writeln!(
            sql,
            "    upper(trim(CAST(\"{ref_col}\" AS VARCHAR))) AS ref,"
        );
    }
    if let Some(alt_col) = &analysis.alt_col {
        let _ = writeln!(
            sql,
            "    upper(trim(CAST(\"{alt_col}\" AS VARCHAR))) AS alt,"
        );
    }

    for amb in &analysis.ambiguous {
        let _ = writeln!(
            sql,
            "    -- FIXME: '{}' is ambiguous — {}",
            amb.column, amb.reason
        );
        if analysis.ref_col.is_none() && analysis.alt_col.is_none() {
            if amb.column.to_lowercase() == "a1" || amb.column.to_lowercase() == "allele1" {
                let _ = writeln!(
                    sql,
                    "    upper(trim(CAST(\"{}\" AS VARCHAR))) AS alt,  -- assumed effect allele",
                    amb.column
                );
            } else {
                let _ = writeln!(
                    sql,
                    "    upper(trim(CAST(\"{}\" AS VARCHAR))) AS ref,  -- assumed other allele",
                    amb.column
                );
            }
        }
    }

    if analysis.ref_col.is_none() && analysis.alt_col.is_none() && analysis.ambiguous.is_empty() {
        let _ = writeln!(sql, "    -- FIXME: no ref/alt columns detected");
    }

    if let Some(rsid_col) = &analysis.rsid_col {
        let _ = writeln!(sql, "    CAST(\"{rsid_col}\" AS VARCHAR) AS rsid,");
    }

    for mapping in &analysis.columns {
        if ["chromosome", "position", "ref", "alt", "rsid"].contains(&mapping.canonical) {
            continue;
        }
        let typ = canonical_type(mapping.canonical);
        let _ = writeln!(
            sql,
            "    CAST(\"{}\" AS {typ}) AS {},",
            mapping.input_name, mapping.canonical
        );
    }

    for col in &analysis.unmapped {
        let _ = writeln!(sql, "    \"{col}\",");
    }

    let _ = writeln!(sql, "FROM _ingest_input");

    if let Some(chr_col) = &analysis.chr_col {
        let _ = writeln!(
            sql,
            "WHERE upper(regexp_replace(CAST(\"{}\" AS VARCHAR), '^chr', '', 'i')) IN (",
            chr_col
        );
        let _ = writeln!(sql, "    '1','2','3','4','5','6','7','8','9','10',");
        let _ = writeln!(
            sql,
            "    '11','12','13','14','15','16','17','18','19','20',"
        );
        let _ = writeln!(sql, "    '21','22','X','Y','MT','M','23','24'");
        let _ = writeln!(sql, ")");
    }

    let _ = writeln!(sql, "ORDER BY position");
}

fn write_parquet_select(sql: &mut String, _analysis: &Analysis) {
    let _ = writeln!(
        sql,
        "-- Parquet re-normalization: inspect schema and rename columns"
    );
    let _ = writeln!(
        sql,
        "-- FIXME: verify column names match your parquet schema"
    );
    let _ = writeln!(sql, "SELECT * FROM _ingest_input");
}

fn write_vcf_stub(sql: &mut String) {
    let _ = writeln!(sql, "-- VCF ingestion uses the streaming parser, not SQL.");
    let _ = writeln!(sql, "-- Run: favor ingest <file.vcf.gz>");
}
