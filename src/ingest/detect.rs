//! Build and coordinate-base detection from a small sample.

use std::path::Path;

use super::{Analysis, BuildGuess, CoordBase, Delimiter};
use crate::config::Config;
use crate::store::annotation::AnnotationDb;
use crate::engine::DfEngine;
use crate::error::CohortError;

const PROBE_CHROMS: &[&str] = &["1", "2", "22", "10", "X"];

/// Run build and coordinate detection on an analysis.
/// Mutates analysis.build_guess and analysis.coord_base in place.
pub fn detect_build_and_coords(
    analysis: &mut Analysis,
    input_path: &Path,
    engine: &DfEngine,
    config: &Config,
) -> Result<(), CohortError> {
    let chr_col = match &analysis.chr_col {
        Some(c) => c.clone(),
        None => return Ok(()),
    };
    let pos_col = match &analysis.pos_col {
        Some(c) => c.clone(),
        None => return Ok(()),
    };

    let ann_db = match AnnotationDb::open(config) {
        Ok(db) => db,
        Err(_) => return Ok(()),
    };

    let delimiter = match analysis.delimiter {
        Some(Delimiter::Tab) => b'\t',
        Some(Delimiter::Comma) => b',',
        Some(Delimiter::Space) => b' ',
        None => {
            engine.register_parquet_file("_ingest_input", input_path)?;
            return probe_all_chroms(engine, &chr_col, &pos_col, &ann_db, analysis);
        }
    };
    engine.register_csv("_ingest_input", input_path, delimiter)?;
    probe_all_chroms(engine, &chr_col, &pos_col, &ann_db, analysis)
}

fn probe_all_chroms(
    engine: &DfEngine,
    chr_col: &str,
    pos_col: &str,
    ann_db: &AnnotationDb,
    analysis: &mut Analysis,
) -> Result<(), CohortError> {
    for probe_chrom in PROBE_CHROMS {
        let probe_parquet = match ann_db.chrom_parquet(probe_chrom) {
            Some(p) => p,
            None => continue,
        };

        let result = probe_chromosome(engine, chr_col, pos_col, probe_chrom, &probe_parquet);

        match result {
            Ok(Some((rate_1based, rate_0based))) => {
                if rate_1based > 0.5 {
                    analysis.coord_base = CoordBase::OneBased;
                } else if rate_0based > 0.5 {
                    analysis.coord_base = CoordBase::ZeroBased;
                }

                let best_rate = rate_1based.max(rate_0based);
                if best_rate > 0.5 {
                    analysis.build_guess = BuildGuess::Hg38;
                } else if best_rate < 0.2 {
                    analysis.build_guess = BuildGuess::Hg19 {
                        match_rate_hg38: best_rate,
                        match_rate_hg19: 1.0 - best_rate,
                    };
                }
                return Ok(());
            }
            Ok(None) => continue,
            Err(_) => continue,
        }
    }
    Ok(())
}

fn probe_chromosome(
    engine: &DfEngine,
    chr_col: &str,
    pos_col: &str,
    chrom: &str,
    annotation_path: &std::path::Path,
) -> Result<Option<(f64, f64)>, CohortError> {
    let ann_table = format!("_ann_probe_{chrom}");
    engine.register_parquet_file(&ann_table, annotation_path)?;

    engine.execute(&format!(
        "CREATE OR REPLACE VIEW _ingest_probe AS \
         SELECT CAST(\"{pos_col}\" AS INT) AS pos \
         FROM _ingest_input \
         WHERE upper(regexp_replace(CAST(\"{chr_col}\" AS VARCHAR), '^chr', '', 'i')) = '{chrom}' \
         LIMIT 100"
    ))?;

    let count = engine.query_scalar("SELECT COUNT(*) FROM _ingest_probe")?;
    if count < 5 {
        return Ok(None);
    }

    let hits_1 = engine.query_scalar(&format!(
        "SELECT COUNT(*) FROM _ingest_probe s \
         WHERE EXISTS (SELECT 1 FROM {ann_table} a WHERE a.position = s.pos)"
    ))?;

    // 0-based input matches if pos + 1 hits the annotation.
    let hits_0 = engine.query_scalar(&format!(
        "SELECT COUNT(*) FROM _ingest_probe s \
         WHERE EXISTS (SELECT 1 FROM {ann_table} a WHERE a.position = s.pos + 1)"
    ))?;

    Ok(Some((
        hits_1 as f64 / count as f64,
        hits_0 as f64 / count as f64,
    )))
}
