//! Normalize variant inputs into canonical parquet.

pub mod detect;
pub mod format;
pub mod gds;
pub mod reader;
pub mod sql;
pub mod vcf;

pub use reader::{RawRecord, VariantReader};

use std::io::{BufRead, BufReader};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::error::CohortError;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum InputFormat {
    Vcf,
    Tabular,
    Parquet,
    Gds,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Delimiter {
    Tab,
    Comma,
    Space,
}

impl Delimiter {
    pub fn char(self) -> char {
        match self {
            Delimiter::Tab => '\t',
            Delimiter::Comma => ',',
            Delimiter::Space => ' ',
        }
    }

    pub fn byte(self) -> u8 {
        self.char() as u8
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum JoinKey {
    ChromPosRefAlt,
    Rsid,
    ChromPos,
}

#[derive(Debug, Clone, Serialize)]
#[serde(tag = "verdict")]
pub enum BuildGuess {
    #[serde(rename = "hg38")]
    Hg38,
    #[serde(rename = "hg19")]
    Hg19 {
        match_rate_hg38: f64,
        match_rate_hg19: f64,
    },
    #[serde(rename = "unknown")]
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CoordBase {
    OneBased,
    ZeroBased,
    Unknown,
}

#[derive(Debug, Clone, Serialize)]
pub struct ColumnMapping {
    pub input_name: String,
    pub canonical: &'static str,
}

#[derive(Debug, Clone, Serialize)]
pub struct Ambiguity {
    pub column: String,
    pub candidates: Vec<&'static str>,
    pub reason: &'static str,
}

#[derive(Debug, Serialize)]
pub struct Analysis {
    pub format: InputFormat,
    pub delimiter: Option<Delimiter>,
    pub raw_columns: Vec<String>,
    pub columns: Vec<ColumnMapping>,
    pub ambiguous: Vec<Ambiguity>,
    pub unmapped: Vec<String>,
    pub join_key: JoinKey,
    pub build_guess: BuildGuess,
    pub coord_base: CoordBase,
    pub chr_col: Option<String>,
    pub pos_col: Option<String>,
    pub ref_col: Option<String>,
    pub alt_col: Option<String>,
    pub rsid_col: Option<String>,
}

impl Analysis {
    pub fn needs_intervention(&self) -> bool {
        !self.ambiguous.is_empty()
            || matches!(self.build_guess, BuildGuess::Hg19 { .. })
            || matches!(self.coord_base, CoordBase::ZeroBased)
    }

    pub fn status(&self) -> &'static str {
        if self.needs_intervention() {
            "needs_edit"
        } else {
            "ok"
        }
    }
}

pub fn detect_format(path: &Path) -> Result<(InputFormat, Option<Delimiter>), CohortError> {
    let name = path
        .file_name()
        .ok_or_else(|| {
            CohortError::Input(format!(
                "Path '{}' has no file name component.",
                path.display()
            ))
        })?
        .to_string_lossy()
        .to_lowercase();

    if name.ends_with(".vcf")
        || name.ends_with(".vcf.gz")
        || name.ends_with(".vcf.bgz")
        || name.ends_with(".bcf")
    {
        return Ok((InputFormat::Vcf, None));
    }

    if name.ends_with(".parquet") {
        return Ok((InputFormat::Parquet, None));
    }

    if name.ends_with(".gds") {
        return Ok((InputFormat::Gds, None));
    }

    if name.ends_with(".tsv")
        || name.ends_with(".tsv.gz")
        || name.ends_with(".csv")
        || name.ends_with(".csv.gz")
        || name.ends_with(".txt")
        || name.ends_with(".txt.gz")
    {
        let delim = sniff_delimiter(path)?;
        return Ok((InputFormat::Tabular, Some(delim)));
    }

    if path.is_file() {
        let delim = sniff_delimiter(path);
        if let Ok(d) = delim {
            return Ok((InputFormat::Tabular, Some(d)));
        }
    }

    Err(CohortError::Input(format!(
        "Cannot detect format for '{}'. Supported: .vcf, .vcf.gz, .tsv, .csv, .parquet, .txt",
        path.display()
    )))
}

pub(crate) fn sniff_delimiter(path: &Path) -> Result<Delimiter, CohortError> {
    let file = std::fs::File::open(path)
        .map_err(|e| CohortError::Input(format!("Cannot open '{}': {e}", path.display())))?;

    if path.to_string_lossy().ends_with(".gz") {
        return Ok(Delimiter::Tab);
    }

    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    reader
        .read_line(&mut first_line)
        .map_err(|e| CohortError::Input(format!("Cannot read '{}': {e}", path.display())))?;

    if first_line.starts_with("##") || first_line.starts_with("#CHROM") {
        return Ok(Delimiter::Tab);
    }

    let tab_count = first_line.matches('\t').count();
    let comma_count = first_line.matches(',').count();

    if tab_count >= 2 {
        Ok(Delimiter::Tab)
    } else if comma_count >= 2 {
        Ok(Delimiter::Comma)
    } else if first_line.split_whitespace().count() >= 3 {
        Ok(Delimiter::Space)
    } else {
        Ok(Delimiter::Tab)
    }
}

/// Read column headers from a tabular file (first non-comment line).
pub fn read_headers(path: &Path, delimiter: Delimiter) -> Result<Vec<String>, CohortError> {
    let file = std::fs::File::open(path)
        .map_err(|e| CohortError::Input(format!("Cannot open '{}': {e}", path.display())))?;

    let reader = BufReader::new(file);

    for line in reader.lines() {
        let line = line.map_err(|e| CohortError::Input(format!("Read error: {e}")))?;
        let trimmed = line.trim();

        // Skip empty lines and VCF meta-headers
        if trimmed.is_empty() || trimmed.starts_with("##") {
            continue;
        }

        // Strip leading # from #CHROM or #chr
        let header_line = trimmed.strip_prefix('#').unwrap_or(trimmed);

        let cols: Vec<String> = match delimiter {
            Delimiter::Tab => header_line
                .split('\t')
                .map(|s| s.trim().to_string())
                .collect(),
            Delimiter::Comma => header_line
                .split(',')
                .map(|s| s.trim().to_string())
                .collect(),
            Delimiter::Space => header_line
                .split_whitespace()
                .map(|s| s.to_string())
                .collect(),
        };

        return Ok(cols);
    }

    Err(CohortError::Input(format!(
        "No header line found in '{}'",
        path.display()
    )))
}

/// Full analysis: detect format, map columns, detect join key.
/// Build detection is separate (requires annotation parquets).
pub fn analyze(path: &Path) -> Result<Analysis, CohortError> {
    let (format, delimiter) = detect_format(path)?;

    match format {
        InputFormat::Vcf => Ok(Analysis {
            format,
            delimiter: None,
            raw_columns: vec![
                "CHROM".into(),
                "POS".into(),
                "ID".into(),
                "REF".into(),
                "ALT".into(),
            ],
            columns: vec![
                ColumnMapping {
                    input_name: "CHROM".into(),
                    canonical: "chromosome",
                },
                ColumnMapping {
                    input_name: "POS".into(),
                    canonical: "position",
                },
                ColumnMapping {
                    input_name: "REF".into(),
                    canonical: "ref",
                },
                ColumnMapping {
                    input_name: "ALT".into(),
                    canonical: "alt",
                },
                ColumnMapping {
                    input_name: "ID".into(),
                    canonical: "rsid",
                },
            ],
            ambiguous: vec![],
            unmapped: vec![],
            join_key: JoinKey::ChromPosRefAlt,
            build_guess: BuildGuess::Unknown,
            coord_base: CoordBase::OneBased, // VCF is always 1-based
            chr_col: Some("CHROM".into()),
            pos_col: Some("POS".into()),
            ref_col: Some("REF".into()),
            alt_col: Some("ALT".into()),
            rsid_col: Some("ID".into()),
        }),

        InputFormat::Parquet => {
            // For parquet, schema is inspected in the command handler
            Ok(Analysis {
                format,
                delimiter: None,
                raw_columns: vec![],
                columns: vec![],
                ambiguous: vec![],
                unmapped: vec![],
                join_key: JoinKey::ChromPosRefAlt,
                build_guess: BuildGuess::Unknown,
                coord_base: CoordBase::Unknown,
                chr_col: None,
                pos_col: None,
                ref_col: None,
                alt_col: None,
                rsid_col: None,
            })
        }

        InputFormat::Gds => Ok(Analysis {
            format,
            delimiter: None,
            raw_columns: vec![
                "chromosome".into(),
                "position".into(),
                "allele".into(),
            ],
            columns: vec![
                ColumnMapping {
                    input_name: "chromosome".into(),
                    canonical: "chromosome",
                },
                ColumnMapping {
                    input_name: "position".into(),
                    canonical: "position",
                },
            ],
            ambiguous: vec![],
            unmapped: vec![],
            join_key: JoinKey::ChromPosRefAlt,
            build_guess: BuildGuess::Unknown,
            coord_base: CoordBase::OneBased,
            chr_col: Some("chromosome".into()),
            pos_col: Some("position".into()),
            ref_col: None,
            alt_col: None,
            rsid_col: None,
        }),

        InputFormat::Tabular => {
            let delim = delimiter.unwrap_or(Delimiter::Tab);
            let raw_columns = read_headers(path, delim)?;

            let (mapped, ambiguous, unmapped) = map_columns(&raw_columns);

            // Extract key column names from the mapping
            let find_col = |canonical: &str| -> Option<String> {
                mapped
                    .iter()
                    .find(|m| m.canonical == canonical)
                    .map(|m| m.input_name.clone())
            };

            let chr_col = find_col("chromosome");
            let pos_col = find_col("position");
            let ref_col = find_col("ref");
            let alt_col = find_col("alt");
            let rsid_col = find_col("rsid");

            let join_key =
                if chr_col.is_some() && pos_col.is_some() && ref_col.is_some() && alt_col.is_some()
                {
                    JoinKey::ChromPosRefAlt
                } else if rsid_col.is_some() {
                    JoinKey::Rsid
                } else if chr_col.is_some() && pos_col.is_some() {
                    JoinKey::ChromPos
                } else {
                    let has_variant_id = raw_columns.iter().any(|c| {
                        let lower = c.to_lowercase();
                        lower == "variant_id" || lower == "varid" || lower == "variant"
                    });
                    if has_variant_id {
                        JoinKey::ChromPosRefAlt
                    } else {
                        JoinKey::Rsid
                    }
                };

            Ok(Analysis {
                format,
                delimiter: Some(delim),
                raw_columns,
                columns: mapped,
                ambiguous,
                unmapped,
                join_key,
                build_guess: BuildGuess::Unknown,
                coord_base: CoordBase::Unknown,
                chr_col,
                pos_col,
                ref_col,
                alt_col,
                rsid_col,
            })
        }
    }
}

/// Single source of truth for GWAS column name normalization.
/// `(lowercased_input_name, canonical_output_name)`. Ambiguities (e.g.
/// `A1` as effect or other allele) are detected separately.
static ALIASES: &[(&str, &str)] = &[
    ("chromosome", "chromosome"),
    ("chr", "chromosome"),
    ("chrom", "chromosome"),
    ("#chrom", "chromosome"),
    ("#chr", "chromosome"),
    ("chr_name", "chromosome"),
    ("chrom_id", "chromosome"),
    ("hg38chrc", "chromosome"),
    ("position", "position"),
    ("pos", "position"),
    ("bp", "position"),
    ("base_pair_location", "position"),
    ("bp_hg38", "position"),
    ("bp_grch38", "position"),
    ("pos_hg38", "position"),
    ("start", "position"),
    ("chromstart", "position"),
    ("beg", "position"),
    ("ref", "ref"),
    ("reference", "ref"),
    ("ref_allele", "ref"),
    ("other_allele", "ref"),
    ("nea", "ref"),
    ("non_effect_allele", "ref"),
    ("alt", "alt"),
    ("alternate", "alt"),
    ("alt_allele", "alt"),
    ("effect_allele", "alt"),
    ("ea", "alt"),
    ("tested_allele", "alt"),
    ("risk_allele", "alt"),
    ("rsid", "rsid"),
    ("rsids", "rsid"),
    ("snp", "rsid"),
    ("snpid", "rsid"),
    ("variant_id", "rsid"),
    ("markername", "rsid"),
    ("rs", "rsid"),
    ("rs_id", "rsid"),
    ("beta", "beta"),
    ("effect", "beta"),
    ("effect_size", "beta"),
    ("b", "beta"),
    ("log_or", "beta"),
    ("se", "se"),
    ("stderr", "se"),
    ("standard_error", "se"),
    ("sebeta", "se"),
    ("pvalue", "pvalue"),
    ("p", "pvalue"),
    ("pval", "pvalue"),
    ("p_value", "pvalue"),
    ("p.value", "pvalue"),
    ("p_val", "pvalue"),
    ("log10p", "neglog10p"),
    ("neglog10p", "neglog10p"),
    ("mlogp", "neglog10p"),
    ("log10_p", "neglog10p"),
    ("nlog10p", "neglog10p"),
    ("z", "zscore"),
    ("zscore", "zscore"),
    ("z_score", "zscore"),
    ("z_stat", "zscore"),
    ("n", "n"),
    ("neff", "n"),
    ("n_eff", "n"),
    ("sample_size", "n"),
    ("n_total", "n"),
    ("pip", "pip"),
    ("posterior_prob", "pip"),
    ("posterior_inclusion_prob", "pip"),
    ("pp", "pip"),
    ("cs", "cs_id"),
    ("cs_id", "cs_id"),
    ("credible_set", "cs_id"),
    ("cs_index", "cs_id"),
    ("or", "odds_ratio"),
    ("odds_ratio", "odds_ratio"),
    ("af", "allele_freq"),
    ("freq", "allele_freq"),
    ("eaf", "allele_freq"),
    ("maf", "allele_freq"),
    ("allele_frequency", "allele_freq"),
    ("effect_allele_frequency", "allele_freq"),
];

/// Columns that are AMBIGUOUS — A1/A2 could be effect/other or ref/alt
/// depending on the study convention. We flag these instead of guessing.
static AMBIGUOUS_ALLELE_COLS: &[&str] = &["a1", "a2", "allele1", "allele2"];

/// Apply the alias map to raw input column names.
/// Returns (mapped, ambiguous, unmapped) — pure function, no side effects.
pub fn map_columns(raw_columns: &[String]) -> (Vec<ColumnMapping>, Vec<Ambiguity>, Vec<String>) {
    map_columns_with(raw_columns, ALIASES, AMBIGUOUS_ALLELE_COLS)
}

// Structured API over the ALIASES map: resolve input columns to canonical names
// and report what's missing, ambiguous, or unmapped.

#[allow(dead_code)] // public API — consumed by future commands
pub struct ResolvedColumns {
    pub mapping: Vec<ColumnMapping>,
    pub unmapped: Vec<String>,
    pub missing_required: Vec<&'static str>,
    pub ambiguous: Vec<Ambiguity>,
}

#[allow(dead_code)] // public API — presets for variant, GWAS, credible set resolution
pub struct ColumnResolver {
    aliases: &'static [(&'static str, &'static str)],
    ambiguous: &'static [&'static str],
    required: &'static [&'static str],
}

#[allow(dead_code)]
impl ColumnResolver {
    /// Standard variant resolver (chrom, pos, ref, alt required).
    pub fn variant() -> Self {
        Self {
            aliases: ALIASES,
            ambiguous: AMBIGUOUS_ALLELE_COLS,
            required: &["chromosome", "position", "ref", "alt"],
        }
    }

    /// GWAS summary statistics (variant + beta + se + pvalue).
    pub fn gwas_sumstats() -> Self {
        Self {
            aliases: ALIASES,
            ambiguous: AMBIGUOUS_ALLELE_COLS,
            required: &[
                "chromosome",
                "position",
                "ref",
                "alt",
                "beta",
                "se",
                "pvalue",
            ],
        }
    }

    /// Credible set (variant + pip + cs_id).
    pub fn credible_set() -> Self {
        Self {
            aliases: ALIASES,
            ambiguous: AMBIGUOUS_ALLELE_COLS,
            required: &["chromosome", "position", "ref", "alt", "pip", "cs_id"],
        }
    }

    pub fn resolve(&self, input_columns: &[String]) -> ResolvedColumns {
        let (mapped, ambiguous, unmapped) =
            map_columns_with(input_columns, self.aliases, self.ambiguous);

        let mapped_canonicals: Vec<&str> = mapped.iter().map(|m| m.canonical).collect();

        let missing_required: Vec<&str> = self
            .required
            .iter()
            .filter(|r| !mapped_canonicals.contains(r))
            .copied()
            .collect();

        ResolvedColumns {
            mapping: mapped,
            unmapped,
            missing_required,
            ambiguous,
        }
    }
}

/// Internal: map columns using a given alias table and ambiguity list.
fn map_columns_with(
    raw_columns: &[String],
    aliases: &'static [(&'static str, &'static str)],
    ambiguous_cols: &[&str],
) -> (Vec<ColumnMapping>, Vec<Ambiguity>, Vec<String>) {
    let mut mapped = Vec::new();
    let mut ambiguous = Vec::new();
    let mut unmapped = Vec::new();

    for col in raw_columns {
        let lower = col.to_lowercase().trim().to_string();

        if ambiguous_cols.contains(&lower.as_str()) {
            ambiguous.push(Ambiguity {
                column: col.clone(),
                candidates: vec!["ref", "alt"],
                reason: "A1/A2 convention varies across studies. \
                         Most common: A1=effect(alt), A2=other(ref). \
                         Check your study's README.",
            });
            continue;
        }

        if let Some((_, canonical)) = aliases.iter().find(|(alias, _)| *alias == lower) {
            mapped.push(ColumnMapping {
                input_name: col.clone(),
                canonical,
            });
        } else {
            unmapped.push(col.clone());
        }
    }

    (mapped, ambiguous, unmapped)
}

// Upfront schema validation: every command checks required columns before compute.

pub struct ColumnRequirement {
    pub name: &'static str,
    pub source: &'static str,
    pub used_by: &'static str,
}

pub struct ColumnContract {
    #[allow(dead_code)] // documentation/debugging field
    pub command: &'static str,
    pub required: &'static [ColumnRequirement],
}

impl ColumnContract {
    /// Check that all required columns exist in the schema.
    /// Returns Vec of missing columns (empty = valid).
    pub fn check(&self, available: &[String]) -> Vec<&ColumnRequirement> {
        self.required
            .iter()
            .filter(|r| !available.iter().any(|c| c == r.name))
            .collect()
    }

    /// Format missing columns into an actionable error message.
    pub fn format_missing(missing: &[&ColumnRequirement]) -> String {
        missing
            .iter()
            .map(|r| {
                format!(
                    "  '{}' (from {}, needed by {})",
                    r.name, r.source, r.used_by
                )
            })
            .collect::<Vec<_>>()
            .join("\n")
    }
}

/// SQL type for CAST in SELECT statements.
pub fn canonical_type(canonical: &str) -> &'static str {
    match canonical {
        "chromosome" => "VARCHAR",
        "position" => "INTEGER",
        "ref" | "alt" | "rsid" => "VARCHAR",
        "beta" | "se" | "pvalue" | "zscore" | "neglog10p" | "pip" | "odds_ratio"
        | "allele_freq" => "DOUBLE",
        "n" | "cs_id" => "INTEGER",
        _ => "VARCHAR",
    }
}
