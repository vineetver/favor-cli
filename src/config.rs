use std::fmt;
use std::path::{Path, PathBuf};
use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::column::Col;
use crate::error::CohortError;

/// Annotation tier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Tier {
    Base,
    Full,
}

/// Deployment environment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Environment {
    Hpc,
    Workstation,
}

/// Pipeline columns available from base-tier annotations.
const BASE_TIER_COLS: &[Col] = &[
    Col::Chromosome,
    Col::Position,
    Col::RefAllele,
    Col::AltAllele,
    Col::Vid,
    Col::GeneName,
    Col::RegionType,
    Col::Consequence,
    Col::CaddPhred,
    Col::IsCcrePromoter,
    Col::IsCcreEnhancer,
];

/// Pipeline columns available from full-tier annotations.
const FULL_TIER_COLS: &[Col] = &[
    Col::Chromosome,
    Col::Position,
    Col::RefAllele,
    Col::AltAllele,
    Col::Vid,
    Col::GeneName,
    Col::RegionType,
    Col::Consequence,
    Col::CaddPhred,
    Col::IsCcrePromoter,
    Col::IsCcreEnhancer,
    Col::Revel,
    Col::IsCagePromoter,
    Col::IsCageEnhancer,
    Col::WCadd,
    Col::WLinsight,
    Col::WFathmmXf,
    Col::WApcEpiActive,
    Col::WApcEpiRepressed,
    Col::WApcEpiTranscription,
    Col::WApcConservation,
    Col::WApcProteinFunction,
    Col::WApcLocalNd,
    Col::WApcMutationDensity,
    Col::WApcTf,
];

impl Tier {
    pub fn as_str(self) -> &'static str {
        match self {
            Tier::Base => "base",
            Tier::Full => "full",
        }
    }

    pub fn size_human(self) -> &'static str {
        match self {
            Tier::Base => "~200 GB",
            Tier::Full => "~508 GB",
        }
    }

    /// Pipeline columns guaranteed present after annotating with this tier.
    pub fn columns(self) -> &'static [Col] {
        match self {
            Tier::Base => BASE_TIER_COLS,
            Tier::Full => FULL_TIER_COLS,
        }
    }

    /// Does this tier provide a specific pipeline column?
    pub fn has(self, col: Col) -> bool {
        self.columns().contains(&col)
    }

    /// What tier is required to produce all of the given columns?
    pub fn required_for(cols: &[Col]) -> Tier {
        if cols.iter().all(|c| Tier::Base.has(*c)) {
            Tier::Base
        } else {
            Tier::Full
        }
    }

    /// Top-level struct columns present in this tier's annotation parquets.
    #[allow(dead_code)]
    pub fn source_columns(self) -> &'static [&'static str] {
        match self {
            Tier::Base => &["gencode", "main", "ccre"],
            Tier::Full => &[
                "gencode",
                "main",
                "ccre",
                "cage",
                "apc",
                "dbnsfp",
                "linsight",
                "fathmm_xf",
            ],
        }
    }

    /// Does this tier's annotation parquet include a given source column?
    #[allow(dead_code)]
    pub fn has_source(self, col: &str) -> bool {
        self.source_columns().contains(&col)
    }
}

impl fmt::Display for Tier {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for Tier {
    type Err = CohortError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "base" => Ok(Tier::Base),
            "full" => Ok(Tier::Full),
            other => Err(CohortError::Input(format!(
                "Invalid tier '{other}'. Must be 'base' or 'full'."
            ))),
        }
    }
}

impl Environment {
    pub fn as_str(self) -> &'static str {
        match self {
            Environment::Hpc => "hpc",
            Environment::Workstation => "workstation",
        }
    }
}

impl fmt::Display for Environment {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.as_str())
    }
}

impl FromStr for Environment {
    type Err = CohortError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "hpc" | "cluster" => Ok(Environment::Hpc),
            "workstation" | "local" | "desktop" => Ok(Environment::Workstation),
            other => Err(CohortError::Input(format!(
                "Invalid environment '{other}'. Must be 'hpc' or 'workstation'."
            ))),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Config {
    pub data: DataConfig,
    #[serde(default)]
    pub resources: ResourceConfig,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DataConfig {
    pub tier: Tier,
    /// Single root directory — everything derives from this
    pub root_dir: String,
    /// Installed add-on pack IDs
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub packs: Vec<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct ResourceConfig {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub environment: Option<Environment>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub memory_budget: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub threads: Option<usize>,
    #[serde(default, skip_serializing_if = "String::is_empty")]
    pub temp_dir: String,
}

impl ResourceConfig {
    /// Parse a memory string like "16GB", "16G", "16gb", "16384MB", "16384M"
    /// into bytes. Returns None if the string is empty or unparseable.
    pub fn parse_memory_bytes(s: &str) -> Option<u64> {
        let s = s.trim();
        if s.is_empty() {
            return None;
        }

        let s_upper = s.to_uppercase();
        let (num_str, multiplier) = if let Some(n) = s_upper.strip_suffix("TB") {
            (n, 1024u64 * 1024 * 1024 * 1024)
        } else if let Some(n) = s_upper.strip_suffix("GB") {
            (n, 1024u64 * 1024 * 1024)
        } else if let Some(n) = s_upper.strip_suffix("MB") {
            (n, 1024u64 * 1024)
        } else if let Some(n) = s_upper.strip_suffix('T') {
            (n, 1024u64 * 1024 * 1024 * 1024)
        } else if let Some(n) = s_upper.strip_suffix('G') {
            (n, 1024u64 * 1024 * 1024)
        } else if let Some(n) = s_upper.strip_suffix('M') {
            (n, 1024u64 * 1024)
        } else {
            // Assume bytes if no suffix
            return s.parse::<u64>().ok();
        };

        num_str
            .trim()
            .parse::<f64>()
            .ok()
            .map(|n| (n * multiplier as f64) as u64)
    }

    /// Resolve memory_budget to bytes. Returns None if not set.
    pub fn memory_budget_bytes(&self) -> Option<u64> {
        self.memory_budget
            .as_deref()
            .and_then(Self::parse_memory_bytes)
    }
}

impl Config {
    pub fn config_dir() -> PathBuf {
        dirs::home_dir()
            .map(|h| h.join(".cohort"))
            .unwrap_or_else(|| PathBuf::from(".cohort"))
    }

    pub fn config_path() -> PathBuf {
        Self::config_dir().join("config.toml")
    }

    pub fn default_root_dir() -> PathBuf {
        Self::config_dir().join("data")
    }

    /// Root directory for all FAVOR data
    pub fn root_dir(&self) -> PathBuf {
        PathBuf::from(&self.data.root_dir)
    }

    /// Annotation parquets: <root>/base/ or <root>/full/
    pub fn annotations_dir(&self) -> PathBuf {
        self.root_dir().join(self.data.tier.as_str())
    }

    /// Tissue enrichment data: <root>/tissue/
    pub fn tissue_dir(&self) -> PathBuf {
        self.root_dir().join("tissue")
    }

    pub fn load() -> Result<Self, CohortError> {
        let path = Self::config_path();
        if !path.exists() {
            return Ok(Self::default());
        }
        let content = std::fs::read_to_string(&path).map_err(|e| {
            CohortError::Resource(format!("Cannot read config '{}': {e}", path.display()))
        })?;
        let config: Config = toml::from_str(&content).map_err(|e| {
            CohortError::Resource(format!("Invalid config '{}': {e}", path.display()))
        })?;
        Ok(config)
    }

    /// Load config and require that setup has been run.
    /// Use this in every command that needs data paths.
    pub fn load_configured() -> Result<Self, CohortError> {
        let config = Self::load()?;
        if !config.is_configured() {
            return Err(CohortError::DataMissing(
                "Not configured. Run `cohort setup` first.".to_string(),
            ));
        }
        Ok(config)
    }

    pub fn save(&self) -> Result<(), CohortError> {
        let path = Self::config_path();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!(
                    "Cannot create config directory '{}': {e}",
                    parent.display()
                ))
            })?;
        }
        let content = toml::to_string_pretty(self)?;
        std::fs::write(&path, &content).map_err(|e| {
            CohortError::Resource(format!("Cannot write config '{}': {e}", path.display()))
        })?;
        Ok(())
    }

    pub fn is_configured(&self) -> bool {
        !self.data.root_dir.is_empty()
    }

    pub fn has_annotations(&self) -> bool {
        Self::count_chromosomes(&self.annotations_dir()) > 0
    }

    pub fn has_tissue(&self) -> bool {
        let td = self.tissue_dir();
        td.is_dir()
            && std::fs::read_dir(&td)
                .into_iter()
                .flatten()
                .flatten()
                .any(|entry| {
                    let name = entry.file_name();
                    let name = name.to_string_lossy();
                    !name.starts_with('.') && (entry.path().is_dir() || name.ends_with(".parquet"))
                })
    }

    /// Count chromosome=*/sorted.parquet files in a directory
    pub fn count_chromosomes(dir: &Path) -> usize {
        if !dir.exists() {
            return 0;
        }
        std::fs::read_dir(dir)
            .into_iter()
            .flatten()
            .flatten()
            .filter(|entry| {
                let name = entry.file_name();
                let name = name.to_string_lossy();
                name.starts_with("chromosome=") && entry.path().join("sorted.parquet").exists()
            })
            .count()
    }
}

/// What FAVOR data exists at a given root directory.
/// Cheap to compute — just stat() calls, no file reads.
pub struct DirProbe {
    pub base_chroms: usize,
    pub full_chroms: usize,
    pub tissue_tables: Vec<String>,
    pub rollups_found: bool,
    pub reference_found: bool,
}

impl DirProbe {
    /// Scan a directory to detect FAVOR data. Fast — only readdir + stat.
    pub fn scan(root: &Path) -> Self {
        let base_chroms = Config::count_chromosomes(&root.join("base"));
        let full_chroms = Config::count_chromosomes(&root.join("full"));

        let tissue_dir = root.join("tissue");
        let tissue_tables = if tissue_dir.is_dir() {
            std::fs::read_dir(&tissue_dir)
                .into_iter()
                .flatten()
                .flatten()
                .filter(|e| e.path().is_dir())
                .filter_map(|e| {
                    let name = e.file_name().into_string().ok()?;
                    if name.starts_with('.') {
                        return None;
                    }
                    Some(name)
                })
                .collect()
        } else {
            Vec::new()
        };

        let rollups_found = root.join("rollups").is_dir();
        let reference_found = tissue_dir.join("reference").is_dir();

        Self {
            base_chroms,
            full_chroms,
            tissue_tables,
            rollups_found,
            reference_found,
        }
    }

    pub fn has_any_data(&self) -> bool {
        self.base_chroms > 0
            || self.full_chroms > 0
            || !self.tissue_tables.is_empty()
            || self.rollups_found
    }

    /// Best tier detected, if any annotations exist
    pub fn detected_tier(&self) -> Option<Tier> {
        if self.full_chroms > 0 {
            Some(Tier::Full)
        } else if self.base_chroms > 0 {
            Some(Tier::Base)
        } else {
            None
        }
    }

    /// Summary lines for display in TUI
    pub fn summary_lines(&self) -> Vec<(String, ProbeStatus)> {
        let mut lines = Vec::new();

        // Annotations
        if self.full_chroms == 24 {
            lines.push(("full: 24/24 chromosomes".into(), ProbeStatus::Good));
        } else if self.full_chroms > 0 {
            lines.push((
                format!("full: {}/24 chromosomes", self.full_chroms),
                ProbeStatus::Partial,
            ));
        }
        if self.base_chroms == 24 {
            lines.push(("base: 24/24 chromosomes".into(), ProbeStatus::Good));
        } else if self.base_chroms > 0 {
            lines.push((
                format!("base: {}/24 chromosomes", self.base_chroms),
                ProbeStatus::Partial,
            ));
        }
        if self.full_chroms == 0 && self.base_chroms == 0 {
            lines.push(("annotations: not found".into(), ProbeStatus::Missing));
        }

        // Tissue packs
        if !self.tissue_tables.is_empty() {
            lines.push((
                format!("tissue: {} tables", self.tissue_tables.len()),
                ProbeStatus::Good,
            ));
        } else {
            lines.push(("tissue: not found".into(), ProbeStatus::Missing));
        }

        // Rollups
        if self.rollups_found {
            lines.push(("rollups: found".into(), ProbeStatus::Good));
        }

        // Reference
        if self.reference_found {
            lines.push(("reference: found".into(), ProbeStatus::Good));
        }

        if !self.has_any_data() {
            lines.push((
                "(data will be downloaded after setup)".into(),
                ProbeStatus::Info,
            ));
        }

        lines
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProbeStatus {
    Good,
    Partial,
    Missing,
    Info,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            data: DataConfig {
                tier: Tier::Base,
                root_dir: String::new(),
                packs: Vec::new(),
            },
            resources: ResourceConfig::default(),
        }
    }
}
