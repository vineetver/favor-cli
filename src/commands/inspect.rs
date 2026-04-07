use serde_json::json;

use parquet::file::reader::FileReader;

use crate::config::{Config, DirProbe, Tier};
use crate::data::AnnotationDb;
use crate::error::CohortError;
use crate::output::Output;

pub fn schema(table: Option<String>, out: &dyn Output) -> Result<(), CohortError> {
    let config = Config::load_configured()?;

    match table.as_deref() {
        None => list_tables(&config, out),
        Some("base") => describe_tier(&config, Tier::Base, out),
        Some("full") => describe_tier(&config, Tier::Full, out),
        Some(name) => describe_tissue_table(&config, name, out),
    }
}

fn list_tables(config: &Config, out: &dyn Output) -> Result<(), CohortError> {
    let probe = DirProbe::scan(&config.root_dir());
    let mut tables = Vec::new();

    for tier in &[Tier::Base, Tier::Full] {
        let chroms = match tier {
            Tier::Base => probe.base_chroms,
            Tier::Full => probe.full_chroms,
        };
        let status = if chroms == 0 {
            "not installed".to_string()
        } else {
            format!("{}/24 chromosomes", chroms)
        };
        let entry = json!({
            "name": tier.as_str(),
            "type": "annotation",
            "status": status,
            "size": tier.size_human(),
        });
        out.status(&format!(
            "{}: {} ({})",
            tier.as_str(),
            status,
            tier.size_human()
        ));
        tables.push(entry);
    }

    let mut tissue_names: Vec<String> = probe.tissue_tables.clone();
    tissue_names.sort();
    for name in &tissue_names {
        let entry = json!({
            "name": name,
            "type": "tissue",
            "status": "installed",
        });
        out.status(&format!("{}: installed", name));
        tables.push(entry);
    }

    if tissue_names.is_empty() {
        out.status(
            "No tissue tables installed. Run `cohort data pull --pack eqtl` to add tissue data.",
        );
    }

    out.result_json(&json!({ "tables": tables }));
    Ok(())
}

fn describe_tier(config: &Config, tier: Tier, out: &dyn Output) -> Result<(), CohortError> {
    let ann_db = AnnotationDb::open_tier(tier, &config.root_dir())?;
    let sample = match ann_db.chrom_parquet("1") {
        Some(p) => p,
        None => {
            return Err(CohortError::DataMissing(format!(
                "{} tier chromosome=1 not found. Run `cohort data pull{}` first.",
                tier,
                if tier == Tier::Full { " --full" } else { "" },
            )))
        }
    };

    let columns = describe_parquet(&sample)?;

    for col in &columns {
        out.status(&format!("{}: {}", col.name, col.col_type));
    }

    out.result_json(&json!({
        "table": tier.as_str(),
        "file": sample.to_string_lossy(),
        "columns": columns.iter()
            .map(|c| json!({"name": c.name, "type": c.col_type}))
            .collect::<Vec<_>>(),
    }));
    Ok(())
}

fn describe_tissue_table(config: &Config, name: &str, out: &dyn Output) -> Result<(), CohortError> {
    let table_dir = config.tissue_dir().join(name);
    if !table_dir.is_dir() {
        return Err(CohortError::DataMissing(format!(
            "Tissue table '{}' not found. Run `cohort schema` to list available tables.",
            name,
        )));
    }

    let sample = find_first_parquet(&table_dir)?;
    let columns = describe_parquet(&sample)?;

    for col in &columns {
        out.status(&format!("{}: {}", col.name, col.col_type));
    }

    out.result_json(&json!({
        "table": name,
        "file": sample.to_string_lossy(),
        "columns": columns.iter()
            .map(|c| json!({"name": c.name, "type": c.col_type}))
            .collect::<Vec<_>>(),
    }));
    Ok(())
}

struct ColumnInfo {
    name: String,
    col_type: String,
}

fn describe_parquet(path: &std::path::Path) -> Result<Vec<ColumnInfo>, CohortError> {
    let file = std::fs::File::open(path)
        .map_err(|e| CohortError::Resource(format!("Cannot open {}: {e}", path.display())))?;
    let reader = parquet::file::reader::SerializedFileReader::new(file)
        .map_err(|e| CohortError::Resource(format!("Bad parquet {}: {e}", path.display())))?;
    let schema = reader.metadata().file_metadata().schema_descr();
    Ok(schema
        .columns()
        .iter()
        .map(|c| ColumnInfo {
            name: c.name().to_string(),
            col_type: format!("{:?}", c.physical_type()),
        })
        .collect())
}

fn find_first_parquet(table_dir: &std::path::Path) -> Result<std::path::PathBuf, CohortError> {
    let entries: Vec<_> = std::fs::read_dir(table_dir)
        .map_err(|e| {
            CohortError::Resource(format!(
                "Cannot read directory '{}': {e}",
                table_dir.display()
            ))
        })?
        .filter_map(|e| e.ok())
        .filter(|e| {
            let name = e.file_name();
            let name = name.to_string_lossy();
            name.starts_with("chrom_id=") && e.path().is_dir()
        })
        .collect();

    for entry in &entries {
        let data_file = entry.path().join("data_0.parquet");
        if data_file.exists() {
            return Ok(data_file);
        }
    }

    Err(CohortError::DataMissing(format!(
        "No parquet files found in {}. The table may be empty or corrupted.",
        table_dir.display(),
    )))
}

pub fn manifest(output: &dyn Output) -> Result<(), CohortError> {
    let config = Config::load_configured()?;

    let has_data = config.has_annotations();
    let has_tissue = config.has_tissue();

    let manifest = json!({
        "commands": [
            {"name": "ingest",    "status": "available",   "description": "Normalize VCF/TSV to canonical parquet with vid"},
            {"name": "annotate",  "status": if has_data { "available" } else { "unavailable" }, "requires": "annotation data", "reason": if !has_data { "run `cohort setup` then `cohort data pull`" } else { "" }},
            {"name": "enrich",    "status": if has_tissue { "available" } else { "unavailable" }, "requires": "tissue data", "reason": if !has_tissue { "no tissue/ directory found in root" } else { "" }},
            {"name": "interpret", "status": if has_data { "available" } else { "unavailable" }, "requires": "annotated variants"},
            {"name": "staar",     "status": if has_data { "available" } else { "unavailable" }, "requires": "genotypes + annotated variants"},
        ],
        "data": {
            "root": config.data.root_dir,
            "tier": config.data.tier.as_str(),
            "annotations_present": has_data,
            "tissue_present": has_tissue,
        }
    });

    output.result_json(&manifest);
    Ok(())
}
