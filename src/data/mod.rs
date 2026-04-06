//! Data layer for packs, variant sets, and annotation parquet.

pub mod publish;
pub mod transfer;

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use parquet::file::reader::FileReader;

use crate::column::Col;
use crate::config::{Config, Tier};
use crate::error::FavorError;
use crate::ingest::JoinKey;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SourceType {
    Tissue,
    Root,
}

pub struct Pack {
    pub id: &'static str,
    pub name: &'static str,
    pub description: &'static str,
    pub source: &'static str,
    pub size_human: &'static str,
    pub size_bytes: u64,
    pub always_installed: bool,
    pub version: &'static str,
    pub source_type: SourceType,
    pub base_dir: &'static str,
    pub tables: &'static [&'static str],
    pub remote_url: &'static str,
}

impl Pack {
    pub fn all() -> &'static [Pack] {
        PACKS
    }

    pub fn find(id: &str) -> Option<&'static Pack> {
        PACKS.iter().find(|p| p.id == id)
    }

    pub fn optional() -> Vec<&'static Pack> {
        PACKS.iter().filter(|p| !p.always_installed).collect()
    }

    pub fn required() -> Vec<&'static Pack> {
        PACKS.iter().filter(|p| p.always_installed).collect()
    }

    pub fn source_dir(&self, favor_root: &Path) -> PathBuf {
        match self.source_type {
            SourceType::Tissue => favor_root.join("tissue"),
            SourceType::Root => favor_root.to_path_buf(),
        }
    }

    pub fn local_dir(&self, data_root: &Path) -> PathBuf {
        if self.base_dir.is_empty() {
            data_root.to_path_buf()
        } else {
            data_root.join(self.base_dir)
        }
    }

    pub fn base_url<'a>(&self, default: &'a str) -> &'a str {
        if self.remote_url.is_empty() {
            default
        } else {
            self.remote_url
        }
    }
}

const GB: u64 = 1024 * 1024 * 1024;
const MB: u64 = 1024 * 1024;

pub static PACKS: &[Pack] = &[
    Pack {
        id: "reference",
        name: "Reference",
        description: "Gene index, cCRE regions, tissue vocabulary, evidence types",
        source: "FAVOR",
        size_human: "40 MB",
        size_bytes: 40 * MB,
        always_installed: true,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &["reference"],
        remote_url: "",
    },
    Pack {
        id: "rollups",
        name: "Rollups",
        description: "Gene-level and 25kb region-level summary statistics",
        source: "FAVOR",
        size_human: "49 MB",
        size_bytes: 49 * MB,
        always_installed: true,
        version: "1",
        source_type: SourceType::Root,
        base_dir: "",
        tables: &["rollups"],
        remote_url: "",
    },
    Pack {
        id: "variant-index",
        name: "Variant Index",
        description: "Pre-computed variant-to-region junction table (required for enrichment)",
        source: "FAVOR",
        size_human: "155 GB",
        size_bytes: 155 * GB,
        always_installed: true,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &["variant_in_region"],
        remote_url: "",
    },
    Pack {
        id: "eqtl",
        name: "Bulk eQTL",
        description: "eQTL, sQTL, apaQTL across 50 tissues + SuSiE fine-mapped",
        source: "GTEx v10",
        size_human: "3 GB",
        size_bytes: 3 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &[
            "variant_eqtl",
            "variant_sqtl",
            "variant_apaqtl",
            "variant_eqtl_susie",
            "variant_eqtl_ccre",
        ],
        remote_url: "",
    },
    Pack {
        id: "eqtl-catalogue",
        name: "eQTL Catalogue",
        description: "Multi-study eQTL from 37 datasets, 127 cell types",
        source: "EBI eQTL Catalogue",
        size_human: "2 GB",
        size_bytes: 2 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &["variant_eqtl_catalogue"],
        remote_url: "",
    },
    Pack {
        id: "sc-eqtl",
        name: "Single-cell eQTL",
        description: "sc-eQTL from immune cells, brain, and 15 DICE cell types",
        source: "OneK1K, DICE, PsychENCODE",
        size_human: "48 GB",
        size_bytes: 48 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &[
            "variant_sc_eqtl",
            "variant_sc_eqtl_dice",
            "variant_sc_eqtl_psychencode",
        ],
        remote_url: "",
    },
    Pack {
        id: "regulatory",
        name: "Regulatory Elements",
        description: "cCRE tissue signals, chromatin states, accessibility peaks",
        source: "ENCODE SCREEN v4, Roadmap Epigenomics, VISTA",
        size_human: "18 GB",
        size_bytes: 18 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &[
            "region_ccre_tissue_signals",
            "region_chromatin_states",
            "region_accessibility_peaks",
            "region_validated_enhancers",
        ],
        remote_url: "",
    },
    Pack {
        id: "enhancer-gene",
        name: "Enhancer-Gene Links",
        description: "Enhancer-gene predictions + experimental cCRE-gene links",
        source: "ABC, EPIraction, ENCODE rE2G, EpiMap, ChIA-PET, SCREEN v4",
        size_human: "12 GB",
        size_bytes: 12 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &[
            "region_enhancer_gene",
            "region_epiraction",
            "region_encode_re2g",
            "region_epimap",
            "linkage_ccre_experimental_chiapet",
            "linkage_ccre_experimental_crispr",
            "linkage_ccre_computational_screen_v4",
        ],
        remote_url: "",
    },
    Pack {
        id: "tissue-scores",
        name: "Tissue-Specific Scores",
        description: "Tissue scores, ChromBPNet, allelic imbalance/methylation",
        source: "TLand, ChromBPNet, ENTEx",
        size_human: "5 GB",
        size_bytes: 5 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &[
            "variant_tissue_scores",
            "variant_tissue_scores_positional",
            "variant_chrombpnet",
            "variant_allelic_imbalance",
            "variant_allelic_methylation",
            "region_ase_ccre",
            "region_chromatin_loops",
        ],
        remote_url: "",
    },
    Pack {
        id: "pgs",
        name: "Polygenic Risk Scores",
        description: "Pre-computed polygenic scores across phenotypes",
        source: "PGS Catalog",
        size_human: "75 GB",
        size_bytes: 75 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &["pgs_scores", "pgs_metadata"],
        remote_url: "",
    },
    Pack {
        id: "genotypes",
        name: "Genotype Matrix",
        description: "Reference genotype matrix for fine-mapping and LD",
        source: "TOPMed",
        size_human: "3 GB",
        size_bytes: 3 * GB,
        always_installed: false,
        version: "1",
        source_type: SourceType::Tissue,
        base_dir: "tissue",
        tables: &["genotype_matrix"],
        remote_url: "",
    },
];

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum VariantSetKind {
    Ingested,
    Annotated { tier: crate::config::Tier },
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
    pub fn open(path: &Path) -> Result<Self, FavorError> {
        let meta_path = path.join("meta.json");
        if !meta_path.exists() {
            return Err(FavorError::Input(format!(
                "Not a variant set: {}. Missing meta.json. Run `favor ingest` to produce one.",
                path.display()
            )));
        }
        let content = std::fs::read_to_string(&meta_path)
            .map_err(|e| FavorError::Input(format!("Cannot read {}: {e}", meta_path.display())))?;
        let meta: VariantMeta = serde_json::from_str(&content).map_err(|e| {
            FavorError::Input(format!("Invalid meta.json in {}: {e}", path.display()))
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

    pub fn count(&self) -> u64 {
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

    pub fn kind(&self) -> Option<&VariantSetKind> {
        self.meta.kind.as_ref()
    }

    pub fn require_annotated(&self) -> Result<(), FavorError> {
        match &self.meta.kind {
            Some(VariantSetKind::Annotated { .. }) => Ok(()),
            Some(VariantSetKind::Ingested) => Err(FavorError::Input(format!(
                "{} is an ingested variant set, not annotated. Run `favor annotate` first.",
                self.root.display()
            ))),
            None => Ok(()),
        }
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
    pub fn new(root: &Path, join_key: JoinKey, source: &str) -> Result<Self, FavorError> {
        if root.join("meta.json").exists() {
            return Err(FavorError::Input(format!(
                "VariantSet already exists at {}. Remove it first.",
                root.display()
            )));
        }
        std::fs::create_dir_all(root)
            .map_err(|e| FavorError::Resource(format!("Cannot create {}: {e}", root.display())))?;
        Ok(Self {
            root: root.to_path_buf(),
            join_key,
            source: source.to_string(),
            columns: None,
            per_chrom: HashMap::new(),
            kind: None,
        })
    }

    pub fn chrom_path(&self, chrom: &str) -> Result<PathBuf, FavorError> {
        let dir = self.root.join(format!("chromosome={chrom}"));
        std::fs::create_dir_all(&dir)
            .map_err(|e| FavorError::Resource(format!("Cannot create {}: {e}", dir.display())))?;
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

    /// Scan written parquet files to populate per-chrom metadata.
    /// Uses parquet file metadata directly — no query engine needed.
    pub fn scan_and_register(&mut self) -> Result<(), FavorError> {
        let entries = std::fs::read_dir(&self.root).map_err(|e| {
            FavorError::Resource(format!("Cannot read {}: {e}", self.root.display()))
        })?;

        let mut first_parquet: Option<PathBuf> = None;

        for entry in entries {
            let entry = entry.map_err(|e| FavorError::Resource(format!("{e}")))?;
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
                std::fs::read_dir(&dir).map_err(|e| FavorError::Resource(format!("{e}")))?;
            for file in files {
                let file = file.map_err(|e| FavorError::Resource(format!("{e}")))?;
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

    pub fn finish(self) -> Result<VariantSet, FavorError> {
        if self.per_chrom.is_empty() {
            return Err(FavorError::Analysis(
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
            .map_err(|e| FavorError::Resource(format!("JSON serialize failed: {e}")))?;
        std::fs::write(&meta_path, json).map_err(|e| {
            FavorError::Resource(format!("Cannot write {}: {e}", meta_path.display()))
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

pub fn parquet_row_count(path: &Path) -> Result<u64, FavorError> {
    let file = std::fs::File::open(path)
        .map_err(|e| FavorError::Resource(format!("Cannot open {}: {e}", path.display())))?;
    let reader = parquet::file::reader::SerializedFileReader::new(file)
        .map_err(|e| FavorError::Resource(format!("Bad parquet {}: {e}", path.display())))?;
    Ok(reader.metadata().file_metadata().num_rows() as u64)
}

pub fn parquet_column_names(path: &Path) -> Result<Vec<String>, FavorError> {
    let file = std::fs::File::open(path)
        .map_err(|e| FavorError::Resource(format!("Cannot open {}: {e}", path.display())))?;
    let reader = parquet::file::reader::SerializedFileReader::new(file)
        .map_err(|e| FavorError::Resource(format!("Bad parquet {}: {e}", path.display())))?;
    let schema = reader.metadata().file_metadata().schema_descr();
    // Top-level fields from root group — returns "gencode" not its leaf children.
    Ok(schema
        .root_schema()
        .get_fields()
        .iter()
        .map(|f| f.name().to_string())
        .collect())
}

pub struct AnnotationDb {
    tier: Tier,
    root: PathBuf,
}

impl AnnotationDb {
    pub fn open(config: &Config) -> Result<Self, FavorError> {
        Self::open_tier(config.data.tier, &config.root_dir())
    }

    pub fn open_tier(tier: Tier, data_root: &Path) -> Result<Self, FavorError> {
        let root = data_root.join(tier.as_str());
        if !root.exists() {
            return Err(FavorError::DataMissing(format!(
                "Annotations not found at {}. Run `favor data pull --tier {}` first.",
                root.display(),
                tier
            )));
        }
        let db = Self { tier, root };
        db.validate_installed()?;
        Ok(db)
    }

    /// Verify that all autosomal chromosomes have parquet files.
    pub fn validate_installed(&self) -> Result<(), FavorError> {
        let mut missing = Vec::new();
        for n in 1..=22 {
            if !self
                .root
                .join(format!("chromosome={n}/sorted.parquet"))
                .exists()
            {
                missing.push(n.to_string());
            }
        }
        for c in ["X", "Y"] {
            if !self
                .root
                .join(format!("chromosome={c}/sorted.parquet"))
                .exists()
            {
                missing.push(c.to_string());
            }
        }
        if !missing.is_empty() {
            return Err(FavorError::DataMissing(format!(
                "Annotation data ({} tier) missing for chromosomes: {}. \
                 Run: favor data pull --tier {}",
                self.tier,
                missing.join(", "),
                self.tier,
            )));
        }
        Ok(())
    }

    pub fn tier(&self) -> Tier {
        self.tier
    }

    pub fn chrom_parquet(&self, chrom: &str) -> Option<PathBuf> {
        let p = self.root.join(format!("chromosome={chrom}/sorted.parquet"));
        p.exists().then_some(p)
    }

    pub fn root(&self) -> &Path {
        &self.root
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotatedSetMeta {
    pub source: String,
    pub tier: Tier,
    pub match_rate: f64,
    pub variant_count: u64,
    pub chromosomes: Vec<String>,
    pub columns: Vec<String>,
}

pub struct AnnotatedSet {
    root: PathBuf,
    meta: AnnotatedSetMeta,
}

impl AnnotatedSet {
    pub fn open(path: &Path) -> Result<Self, FavorError> {
        let vs = VariantSet::open(path)?;
        let tier = match vs.kind() {
            Some(VariantSetKind::Annotated { tier }) => *tier,
            Some(VariantSetKind::Ingested) => {
                return Err(FavorError::Input(format!(
                    "'{}' is an ingested variant set, not annotated. \
                     Run `favor annotate` first.",
                    path.display()
                )));
            }
            None => Tier::Full, // legacy sets without kind tag — assume full
        };
        let chromosomes = vs
            .chromosomes()
            .into_iter()
            .map(|s| s.to_string())
            .collect();
        Ok(Self {
            root: path.to_path_buf(),
            meta: AnnotatedSetMeta {
                source: path.display().to_string(),
                tier,
                match_rate: 1.0,
                variant_count: vs.count(),
                chromosomes,
                columns: vs.columns().to_vec(),
            },
        })
    }

    pub fn tier(&self) -> Tier {
        self.meta.tier
    }
    pub fn root(&self) -> &Path {
        &self.root
    }
    pub fn variant_count(&self) -> u64 {
        self.meta.variant_count
    }
    #[allow(dead_code)]
    pub fn columns(&self) -> &[String] {
        &self.meta.columns
    }
    #[allow(dead_code)]
    pub fn chromosomes(&self) -> &[String] {
        &self.meta.chromosomes
    }

    /// Validate that this annotated set supports the given pipeline columns.
    pub fn supports(&self, required: &[Col]) -> Result<(), FavorError> {
        for &col in required {
            if !self.tier().has(col) {
                return Err(FavorError::Input(format!(
                    "Column '{}' requires {} tier but this data was annotated \
                     with {} tier.\nRe-annotate with: favor annotate {} --full",
                    col.as_str(),
                    Tier::required_for(&[col]),
                    self.tier(),
                    self.meta.source,
                )));
            }
        }
        Ok(())
    }

    /// Check that a specific column name exists in the output schema.
    #[allow(dead_code)]
    pub fn has_column(&self, name: &str) -> bool {
        self.meta.columns.iter().any(|c| c == name)
    }
}

#[allow(dead_code)] // enrichment pipeline (v0.2.0)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnrichedSetMeta {
    pub source: String,
    pub tier: Tier,
    pub tissue: String,
    pub resolved_tissues: Vec<String>,
    pub tables_written: Vec<(String, i64)>,
    pub input_count: i64,
}

#[allow(dead_code)] // enrichment pipeline (v0.2.0)
pub struct EnrichedSet {
    root: PathBuf,
    meta: EnrichedSetMeta,
}

#[allow(dead_code)] // enrichment pipeline (v0.2.0)
impl EnrichedSet {
    pub fn new(root: PathBuf, meta: EnrichedSetMeta) -> Self {
        Self { root, meta }
    }

    pub fn root(&self) -> &Path {
        &self.root
    }
    pub fn tables_found(&self) -> &[(String, i64)] {
        &self.meta.tables_written
    }
    pub fn tissue(&self) -> &str {
        &self.meta.tissue
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum TissueTable {
    Eqtl,
    Sqtl,
    ApaQtl,
    EqtlSusie,
    EqtlCatalogue,
    ScEqtl,
    ScEqtlDice,
    ScEqtlPsychencode,
    TissueScores,
    ChromBpnet,
    AllelicImbalance,
    AllelicMethylation,
    EqtlCcre,
}

impl TissueTable {
    pub fn all() -> &'static [TissueTable] {
        &[
            TissueTable::Eqtl,
            TissueTable::Sqtl,
            TissueTable::ApaQtl,
            TissueTable::EqtlSusie,
            TissueTable::EqtlCatalogue,
            TissueTable::ScEqtl,
            TissueTable::ScEqtlDice,
            TissueTable::ScEqtlPsychencode,
            TissueTable::TissueScores,
            TissueTable::ChromBpnet,
            TissueTable::AllelicImbalance,
            TissueTable::AllelicMethylation,
            TissueTable::EqtlCcre,
        ]
    }

    pub fn dir_name(self) -> &'static str {
        match self {
            Self::Eqtl => "variant_eqtl",
            Self::Sqtl => "variant_sqtl",
            Self::ApaQtl => "variant_apaqtl",
            Self::EqtlSusie => "variant_eqtl_susie",
            Self::EqtlCatalogue => "variant_eqtl_catalogue",
            Self::ScEqtl => "variant_sc_eqtl",
            Self::ScEqtlDice => "variant_sc_eqtl_dice",
            Self::ScEqtlPsychencode => "variant_sc_eqtl_psychencode",
            Self::TissueScores => "variant_tissue_scores",
            Self::ChromBpnet => "variant_chrombpnet",
            Self::AllelicImbalance => "variant_allelic_imbalance",
            Self::AllelicMethylation => "variant_allelic_methylation",
            Self::EqtlCcre => "variant_eqtl_ccre",
        }
    }

    pub fn display_name(self) -> &'static str {
        match self {
            Self::Eqtl => "eqtl",
            Self::Sqtl => "sqtl",
            Self::ApaQtl => "apaqtl",
            Self::EqtlSusie => "eqtl_susie",
            Self::EqtlCatalogue => "eqtl_catalogue",
            Self::ScEqtl => "sc_eqtl",
            Self::ScEqtlDice => "sc_eqtl_dice",
            Self::ScEqtlPsychencode => "sc_eqtl_psychencode",
            Self::TissueScores => "tissue_scores",
            Self::ChromBpnet => "chrombpnet",
            Self::AllelicImbalance => "allelic_imbalance",
            Self::AllelicMethylation => "allelic_methylation",
            Self::EqtlCcre => "eqtl_ccre",
        }
    }

    pub fn has_tissue_filter(self) -> bool {
        true // all current tables have tissue filtering
    }

    #[allow(dead_code)]
    pub fn required_pack(self) -> &'static str {
        match self {
            Self::Eqtl | Self::Sqtl | Self::ApaQtl | Self::EqtlSusie | Self::EqtlCcre => "eqtl",
            Self::EqtlCatalogue => "eqtl-catalogue",
            Self::ScEqtl | Self::ScEqtlDice | Self::ScEqtlPsychencode => "sc-eqtl",
            Self::TissueScores
            | Self::ChromBpnet
            | Self::AllelicImbalance
            | Self::AllelicMethylation => "tissue-scores",
        }
    }
}

pub struct TissueDb {
    root: PathBuf,
    available: Vec<TissueTable>,
}

impl TissueDb {
    pub fn open(tissue_dir: &Path) -> Result<Self, FavorError> {
        if !tissue_dir.exists() {
            return Err(FavorError::DataMissing(format!(
                "Tissue data not found at {}. Run `favor data pull --pack eqtl` first.",
                tissue_dir.display()
            )));
        }
        let available = TissueTable::all()
            .iter()
            .filter(|t| tissue_dir.join(t.dir_name()).is_dir())
            .copied()
            .collect();
        Ok(Self {
            root: tissue_dir.to_path_buf(),
            available,
        })
    }

    #[allow(dead_code)]
    pub fn root(&self) -> &Path {
        &self.root
    }
    pub fn available_tables(&self) -> &[TissueTable] {
        &self.available
    }

    pub fn table_path(&self, table: TissueTable) -> PathBuf {
        self.root.join(table.dir_name())
    }

    #[allow(dead_code)]
    pub fn validate_tables(&self, needed: &[TissueTable]) -> Result<(), FavorError> {
        let missing: Vec<_> = needed
            .iter()
            .filter(|t| !self.available.contains(t))
            .collect();
        if !missing.is_empty() {
            return Err(FavorError::DataMissing(format!(
                "Tissue data missing: {}. Install with:\n{}",
                missing
                    .iter()
                    .map(|t| t.display_name())
                    .collect::<Vec<_>>()
                    .join(", "),
                missing
                    .iter()
                    .map(|t| format!("  favor data pull --pack {}", t.required_pack()))
                    .collect::<Vec<_>>()
                    .join("\n")
            )));
        }
        Ok(())
    }
}
