use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use crate::config::{DirProbe, Tier};
use crate::data::{VariantSet, VariantSetKind};
use crate::ingest::format::FormatRegistry;
use crate::ingest::InputFormat;

use crate::tui::theme;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ArtifactKind {
    RawVcf,
    PhenotypeTsv,
    KinshipTsv,
    ParquetFile,
    IngestedSet,
    AnnotatedSet { tier: Tier },
    GenotypeStore,
    StaarResults,
    AnnotationRoot,
}

impl ArtifactKind {
    pub fn title(&self) -> &'static str {
        match self {
            Self::RawVcf => "VCF",
            Self::PhenotypeTsv => "Phenotype TSV",
            Self::KinshipTsv => "Kinship TSV",
            Self::ParquetFile => "Parquet file",
            Self::IngestedSet => "Ingested set",
            Self::AnnotatedSet { .. } => "Annotated set",
            Self::GenotypeStore => "Genotype store",
            Self::StaarResults => "STAAR results",
            Self::AnnotationRoot => "Annotation root",
        }
    }

    pub fn glyph(&self) -> &'static str {
        match self {
            Self::RawVcf => theme::GLYPH_VCF,
            Self::PhenotypeTsv => theme::GLYPH_PHENO,
            Self::KinshipTsv => theme::GLYPH_KINSHIP,
            Self::ParquetFile => theme::GLYPH_PARQUET,
            Self::IngestedSet => theme::GLYPH_INGESTED,
            Self::AnnotatedSet { .. } => theme::GLYPH_ANNOTATED,
            Self::GenotypeStore => theme::GLYPH_GENO_STORE,
            Self::StaarResults => theme::GLYPH_STAAR,
            Self::AnnotationRoot => theme::GLYPH_ANNOROOT,
        }
    }

}

pub fn group_order(kind: &ArtifactKind) -> u8 {
    match kind {
        ArtifactKind::AnnotationRoot => 0,
        ArtifactKind::AnnotatedSet { .. } => 1,
        ArtifactKind::IngestedSet => 2,
        ArtifactKind::GenotypeStore => 3,
        ArtifactKind::StaarResults => 4,
        ArtifactKind::RawVcf => 5,
        ArtifactKind::PhenotypeTsv => 6,
        ArtifactKind::KinshipTsv => 7,
        ArtifactKind::ParquetFile => 8,
    }
}

#[derive(Debug, Clone)]
pub struct Artifact {
    pub path: PathBuf,
    pub kind: ArtifactKind,
    pub size_bytes: u64,
    pub display_name: String,
}

pub fn classify(path: &Path) -> Option<Artifact> {
    let display_name = path
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| path.to_string_lossy().into_owned());

    if path.is_file() {
        if !has_known_extension(path) {
            return None;
        }
        let size_bytes = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
        let detected = FormatRegistry::new().detect(path).ok()?;
        let kind = match detected.format {
            InputFormat::Vcf => ArtifactKind::RawVcf,
            InputFormat::Parquet => ArtifactKind::ParquetFile,
            InputFormat::Tabular => classify_tabular(path)?,
        };
        return Some(Artifact { path: path.to_path_buf(), kind, size_bytes, display_name });
    }

    if path.is_dir() {
        if path.join(".genotype_store").join("manifest.json").exists() {
            return Some(Artifact {
                path: path.to_path_buf(),
                kind: ArtifactKind::GenotypeStore,
                size_bytes: 0,
                display_name,
            });
        }
        if path.join("staar.meta.json").exists() {
            return Some(Artifact {
                path: path.to_path_buf(),
                kind: ArtifactKind::StaarResults,
                size_bytes: 0,
                display_name,
            });
        }
        if path.join("meta.json").exists() {
            if let Ok(set) = VariantSet::open(path) {
                let kind = match set.kind() {
                    Some(VariantSetKind::Annotated { tier }) => {
                        ArtifactKind::AnnotatedSet { tier: *tier }
                    }
                    Some(VariantSetKind::Ingested) | None => ArtifactKind::IngestedSet,
                };
                return Some(Artifact {
                    path: path.to_path_buf(),
                    kind,
                    size_bytes: 0,
                    display_name,
                });
            }
        }
        if DirProbe::scan(path).has_any_data() {
            return Some(Artifact {
                path: path.to_path_buf(),
                kind: ArtifactKind::AnnotationRoot,
                size_bytes: 0,
                display_name,
            });
        }
    }

    None
}

fn classify_tabular(path: &Path) -> Option<ArtifactKind> {
    let file = File::open(path).ok()?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();
    reader.read_line(&mut header).ok()?;
    let header = header.trim_end_matches(['\r', '\n']);
    let delim = if header.contains('\t') { '\t' } else { ',' };
    let cols: Vec<&str> = header.split(delim).collect();
    if cols.is_empty() {
        return None;
    }

    let first_lower = cols[0].trim().to_lowercase();
    let is_id = matches!(first_lower.as_str(), "sample_id" | "iid" | "sample" | "id");
    if is_id && cols.len() > 1 {
        return Some(ArtifactKind::PhenotypeTsv);
    }

    if cols.len() == 3 {
        let mut data = String::new();
        if reader.read_line(&mut data).ok()? > 0 {
            let data = data.trim_end_matches(['\r', '\n']);
            let fields: Vec<&str> = data.split(delim).collect();
            if fields.len() == 3 && fields[2].trim().parse::<f64>().is_ok() {
                return Some(ArtifactKind::KinshipTsv);
            }
        }
    }

    None
}


fn has_known_extension(path: &Path) -> bool {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    let lower = name.to_ascii_lowercase();
    lower.ends_with(".vcf")
        || lower.ends_with(".vcf.gz")
        || lower.ends_with(".vcf.bgz")
        || lower.ends_with(".bcf")
        || lower.ends_with(".parquet")
        || lower.ends_with(".tsv")
        || lower.ends_with(".tsv.gz")
        || lower.ends_with(".csv")
        || lower.ends_with(".csv.gz")
        || lower.ends_with(".txt")
}

pub fn keep_entry(path: &Path, depth: usize) -> bool {
    if depth == 0 {
        return true;
    }
    let Some(name) = path.file_name().map(|n| n.to_string_lossy().into_owned()) else {
        return true;
    };
    if name == ".git" || name == "target" {
        return false;
    }
    if name == ".genotype_store" {
        return true;
    }
    if name.starts_with('.') && path.is_dir() {
        return false;
    }
    true
}
