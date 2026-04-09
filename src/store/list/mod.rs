//! Chromosome-partitioned variant containers.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use parquet::file::reader::FileReader;
use serde::{Deserialize, Serialize};

use crate::column::Col;
use crate::config::Tier;
use crate::error::CohortError;
use crate::ingest::JoinKey;

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum VariantSetKind {
    Ingested,
    Annotated { tier: Tier },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VariantMeta {
    pub version: u32,
    pub join_key: JoinKey,
    pub variant_count: u64,
    pub per_chrom: HashMap<String, ChromMeta>,
    pub columns: Vec<String>,
    pub source: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub kind: Option<VariantSetKind>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromMeta {
    pub variant_count: u64,
    pub size_bytes: u64,
}

pub struct VariantSet {
    root: PathBuf,
    meta: VariantMeta,
}

impl VariantSet {
    pub fn open(path: &Path) -> Result<Self, CohortError> {
        let meta_path = path.join("meta.json");
        if !meta_path.exists() {
            return Err(CohortError::Input(format!(
                "Not a variant set: {}. Missing meta.json. Run `cohort ingest` to produce one.",
                path.display()
            )));
        }
        let content = std::fs::read_to_string(&meta_path)
            .map_err(|e| CohortError::Input(format!("Cannot read {}: {e}", meta_path.display())))?;
        let meta: VariantMeta = serde_json::from_str(&content).map_err(|e| {
            CohortError::Input(format!("Invalid meta.json in {}: {e}", path.display()))
        })?;
        Ok(Self {
            root: path.to_path_buf(),
            meta,
        })
    }

    pub fn chrom_dir(&self, chrom: &str) -> PathBuf {
        self.root.join(format!("chromosome={chrom}"))
    }

    pub fn chromosomes(&self) -> Vec<&str> {
        let mut chroms: Vec<&str> = self.meta.per_chrom.keys().map(|s| s.as_str()).collect();
        chroms.sort_by_key(|c| chrom_sort_key(c));
        chroms
    }

    pub fn variant_count(&self) -> u64 {
        self.meta.variant_count
    }

    pub fn columns(&self) -> &[String] {
        &self.meta.columns
    }

    pub fn join_key(&self) -> JoinKey {
        self.meta.join_key
    }
    pub fn root(&self) -> &Path {
        &self.root
    }

    pub fn tier(&self) -> Option<Tier> {
        match &self.meta.kind {
            Some(VariantSetKind::Annotated { tier }) => Some(*tier),
            _ => None,
        }
    }

    pub fn require_annotated(&self) -> Result<Tier, CohortError> {
        match &self.meta.kind {
            Some(VariantSetKind::Annotated { tier }) => Ok(*tier),
            Some(VariantSetKind::Ingested) => Err(CohortError::Input(format!(
                "{} is an ingested variant set, not annotated. Run `cohort annotate` first.",
                self.root.display()
            ))),
            // Legacy meta predates the kind tag; preserve AnnotatedSet's old default to Full.
            None => Ok(Tier::Full),
        }
    }

    pub fn supports(&self, required: &[Col]) -> Result<Tier, CohortError> {
        let tier = self.require_annotated()?;
        for &col in required {
            if !tier.has(col) {
                return Err(CohortError::Input(format!(
                    "Column '{}' requires {} tier but '{}' was annotated \
                     with {} tier.\nRe-annotate with: cohort annotate {} --full",
                    col.as_str(),
                    Tier::required_for(&[col]),
                    self.root.display(),
                    tier,
                    self.root.display(),
                )));
            }
        }
        Ok(tier)
    }
}

pub struct VariantSetWriter {
    root: PathBuf,
    join_key: JoinKey,
    source: String,
    columns: Option<Vec<String>>,
    per_chrom: HashMap<String, ChromMeta>,
    kind: Option<VariantSetKind>,
}

impl VariantSetWriter {
    pub fn new(root: &Path, join_key: JoinKey, source: &str) -> Result<Self, CohortError> {
        if root.join("meta.json").exists() {
            return Err(CohortError::Input(format!(
                "VariantSet already exists at {}. Remove it first.",
                root.display()
            )));
        }
        std::fs::create_dir_all(root)
            .map_err(|e| CohortError::Resource(format!("Cannot create {}: {e}", root.display())))?;
        Ok(Self {
            root: root.to_path_buf(),
            join_key,
            source: source.to_string(),
            columns: None,
            per_chrom: HashMap::new(),
            kind: None,
        })
    }

    pub fn chrom_path(&self, chrom: &str) -> Result<PathBuf, CohortError> {
        let dir = self.root.join(format!("chromosome={chrom}"));
        std::fs::create_dir_all(&dir)
            .map_err(|e| CohortError::Resource(format!("Cannot create {}: {e}", dir.display())))?;
        Ok(dir.join("data.parquet"))
    }

    pub fn register_chrom(&mut self, chrom: &str, variant_count: u64, size_bytes: u64) {
        self.per_chrom.insert(
            chrom.to_string(),
            ChromMeta {
                variant_count,
                size_bytes,
            },
        );
    }

    pub fn set_columns(&mut self, columns: Vec<String>) {
        self.columns = Some(columns);
    }
    pub fn set_kind(&mut self, kind: VariantSetKind) {
        self.kind = Some(kind);
    }

    pub fn scan_and_register(&mut self) -> Result<(), CohortError> {
        let entries = std::fs::read_dir(&self.root).map_err(|e| {
            CohortError::Resource(format!("Cannot read {}: {e}", self.root.display()))
        })?;

        let mut first_parquet: Option<PathBuf> = None;

        for entry in entries {
            let entry = entry.map_err(|e| CohortError::Resource(format!("{e}")))?;
            let name = entry.file_name().to_string_lossy().to_string();
            let chrom = match name.strip_prefix("chromosome=") {
                Some(c) if !c.is_empty() => c.to_string(),
                _ => continue,
            };
            if !entry.file_type().is_ok_and(|t| t.is_dir()) {
                continue;
            }

            let dir = entry.path();
            let mut total_count: u64 = 0;
            let mut total_size: u64 = 0;

            let files =
                std::fs::read_dir(&dir).map_err(|e| CohortError::Resource(format!("{e}")))?;
            for file in files {
                let file = file.map_err(|e| CohortError::Resource(format!("{e}")))?;
                let fname = file.file_name().to_string_lossy().to_string();
                if !fname.ends_with(".parquet") {
                    continue;
                }
                let fpath = file.path();
                total_size += std::fs::metadata(&fpath).map_or(0, |m| m.len());
                total_count += parquet_row_count(&fpath)?;
                if first_parquet.is_none() {
                    first_parquet = Some(fpath);
                }
            }

            if total_count > 0 {
                self.register_chrom(&chrom, total_count, total_size);
            }
        }

        if self.columns.is_none() {
            if let Some(path) = first_parquet {
                let cols = parquet_column_names(&path)?;
                if !cols.is_empty() {
                    self.columns = Some(cols);
                }
            }
        }

        Ok(())
    }

    pub fn finish(self) -> Result<VariantSet, CohortError> {
        if self.per_chrom.is_empty() {
            return Err(CohortError::Analysis(
                "VariantSet has no chromosomes. No variants were written.".into(),
            ));
        }
        let variant_count: u64 = self.per_chrom.values().map(|m| m.variant_count).sum();
        let meta = VariantMeta {
            version: 1,
            join_key: self.join_key,
            variant_count,
            per_chrom: self.per_chrom,
            columns: self.columns.unwrap_or_default(),
            source: self.source,
            kind: self.kind,
        };
        let meta_path = self.root.join("meta.json");
        let json = serde_json::to_string_pretty(&meta)
            .map_err(|e| CohortError::Resource(format!("JSON serialize failed: {e}")))?;
        std::fs::write(&meta_path, json).map_err(|e| {
            CohortError::Resource(format!("Cannot write {}: {e}", meta_path.display()))
        })?;
        Ok(VariantSet {
            root: self.root,
            meta,
        })
    }
}

fn chrom_sort_key(chrom: &str) -> (u8, u8) {
    match chrom {
        "1" => (0, 1),
        "2" => (0, 2),
        "3" => (0, 3),
        "4" => (0, 4),
        "5" => (0, 5),
        "6" => (0, 6),
        "7" => (0, 7),
        "8" => (0, 8),
        "9" => (0, 9),
        "10" => (0, 10),
        "11" => (0, 11),
        "12" => (0, 12),
        "13" => (0, 13),
        "14" => (0, 14),
        "15" => (0, 15),
        "16" => (0, 16),
        "17" => (0, 17),
        "18" => (0, 18),
        "19" => (0, 19),
        "20" => (0, 20),
        "21" => (0, 21),
        "22" => (0, 22),
        "X" => (1, 0),
        "Y" => (1, 1),
        "MT" => (1, 2),
        _ => (2, 0),
    }
}

pub fn parquet_row_count(path: &Path) -> Result<u64, CohortError> {
    let file = std::fs::File::open(path)
        .map_err(|e| CohortError::Resource(format!("Cannot open {}: {e}", path.display())))?;
    let reader = parquet::file::reader::SerializedFileReader::new(file)
        .map_err(|e| CohortError::Resource(format!("Bad parquet {}: {e}", path.display())))?;
    Ok(reader.metadata().file_metadata().num_rows() as u64)
}

pub fn parquet_column_names(path: &Path) -> Result<Vec<String>, CohortError> {
    let file = std::fs::File::open(path)
        .map_err(|e| CohortError::Resource(format!("Cannot open {}: {e}", path.display())))?;
    let reader = parquet::file::reader::SerializedFileReader::new(file)
        .map_err(|e| CohortError::Resource(format!("Bad parquet {}: {e}", path.display())))?;
    let schema = reader.metadata().file_metadata().schema_descr();
    Ok(schema
        .root_schema()
        .get_fields()
        .iter()
        .map(|f| f.name().to_string())
        .collect())
}

