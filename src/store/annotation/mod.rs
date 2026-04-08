//! Read-only views over FAVOR annotation data. `AnnotationDb` is the
//! base/full tier reader (one parquet per chromosome); `TissueDb` is the
//! tissue add-on reader (one directory per table). Both point at data
//! living outside the store root — the store does not own the bytes.

pub mod refs;

pub use refs::{open_favor_tier, AnnotationKind, AnnotationRef, AnnotationRegistry};

use std::path::{Path, PathBuf};

use crate::config::{Config, Tier};
use crate::error::CohortError;

pub struct AnnotationDb {
    tier: Tier,
    root: PathBuf,
}

impl AnnotationDb {
    pub fn open(config: &Config) -> Result<Self, CohortError> {
        Self::open_tier(config.data.tier, &config.root_dir())
    }

    pub fn open_tier(tier: Tier, data_root: &Path) -> Result<Self, CohortError> {
        let root = data_root.join(tier.as_str());
        if !root.exists() {
            return Err(CohortError::DataMissing(format!(
                "Annotations not found at {}. Run `cohort data pull --tier {}` first.",
                root.display(),
                tier
            )));
        }
        let db = Self { tier, root };
        db.validate_installed()?;
        Ok(db)
    }

    /// Verify that all autosomal chromosomes have parquet files.
    pub fn validate_installed(&self) -> Result<(), CohortError> {
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
            return Err(CohortError::DataMissing(format!(
                "Annotation data ({} tier) missing for chromosomes: {}. \
                 Run: cohort data pull --tier {}",
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
        true
    }
}

pub struct TissueDb {
    root: PathBuf,
    available: Vec<TissueTable>,
}

impl TissueDb {
    pub fn open(tissue_dir: &Path) -> Result<Self, CohortError> {
        if !tissue_dir.exists() {
            return Err(CohortError::DataMissing(format!(
                "Tissue data not found at {}. Run `cohort data pull --pack eqtl` first.",
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

    pub fn available_tables(&self) -> &[TissueTable] {
        &self.available
    }

    pub fn table_path(&self, table: TissueTable) -> PathBuf {
        self.root.join(table.dir_name())
    }
}
