//! Static catalogue of installable annotation packs. Consumed by
//! `data::transfer`, `data::publish`, and `setup`.

pub mod publish;
pub mod transfer;

use std::path::{Path, PathBuf};

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
    #[allow(dead_code)]
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
