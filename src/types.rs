//! Canonical object model for FAVOR.
//!
//! One struct per concept. Downstream code composes these types — no parallel
//! structs, no field duplication, no `Vec<f64>` where `[f64; 11]` suffices.

use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

// ---------------------------------------------------------------------------
// Chromosome
// ---------------------------------------------------------------------------

/// Strongly-typed chromosome identifier.
///
/// `Ord` gives natural sort order: 1..22 < X < Y < MT.
/// Replaces all string-based chromosome fields and `chrom_sort_key()`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Chromosome {
    Autosome(u8), // 1–22
    X,
    Y,
    MT,
}

impl Chromosome {
    fn sort_key(self) -> (u8, u8) {
        match self {
            Chromosome::Autosome(n) => (0, n),
            Chromosome::X => (1, 0),
            Chromosome::Y => (1, 1),
            Chromosome::MT => (1, 2),
        }
    }
}

impl PartialOrd for Chromosome {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Chromosome {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.sort_key().cmp(&other.sort_key())
    }
}

impl Chromosome {
    /// Bare label without "chr" prefix — matches parquet conventions.
    pub fn label(&self) -> String {
        match self {
            Chromosome::Autosome(n) => n.to_string(),
            Chromosome::X => "X".into(),
            Chromosome::Y => "Y".into(),
            Chromosome::MT => "MT".into(),
        }
    }
}

impl fmt::Display for Chromosome {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Chromosome::Autosome(n) => write!(f, "chr{n}"),
            Chromosome::X => write!(f, "chrX"),
            Chromosome::Y => write!(f, "chrY"),
            Chromosome::MT => write!(f, "chrMT"),
        }
    }
}

impl FromStr for Chromosome {
    type Err = String;

    /// Parse chromosome from common representations:
    /// "chr1", "1", "chrX", "X", "chrM", "MT", "chrMT", etc.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let raw = s.strip_prefix("chr").unwrap_or(s);
        match raw {
            "X" | "x" => Ok(Chromosome::X),
            "Y" | "y" => Ok(Chromosome::Y),
            "M" | "MT" | "Mt" | "m" | "mt" => Ok(Chromosome::MT),
            _ => {
                let n: u8 = raw
                    .parse()
                    .map_err(|_| format!("invalid chromosome: {s}"))?;
                if (1..=22).contains(&n) {
                    Ok(Chromosome::Autosome(n))
                } else {
                    Err(format!("autosome out of range 1–22: {n}"))
                }
            }
        }
    }
}

impl Serialize for Chromosome {
    fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.to_string())
    }
}

impl<'de> Deserialize<'de> for Chromosome {
    fn deserialize<D: serde::Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let s = String::deserialize(deserializer)?;
        Chromosome::from_str(&s).map_err(serde::de::Error::custom)
    }
}

// ---------------------------------------------------------------------------
// AnnotationWeights
// ---------------------------------------------------------------------------

/// Fixed-size annotation weights for the 11 STAAR channels.
///
/// Replaces `Vec<f64>` annotation_weights, `WEIGHT_COL_START` / `WEIGHT_COL_COUNT`
/// magic constants, and `ANNOTATION_CHANNELS` (now `DISPLAY_NAMES`).
#[derive(Clone, Copy, Debug)]
pub struct AnnotationWeights(pub [f64; 11]);

impl Default for AnnotationWeights {
    fn default() -> Self {
        Self([0.0; 11])
    }
}

impl AnnotationWeights {
    /// Column names as written to parquet.
    /// Canonical source of truth is now `column::STAAR_WEIGHTS` — these
    /// constants are retained for backward-compat tests only.
    #[cfg(test)]
    pub const NAMES: [&str; 11] = [
        "w_cadd",
        "w_linsight",
        "w_fathmm_xf",
        "w_apc_epi_active",
        "w_apc_epi_repressed",
        "w_apc_epi_transcription",
        "w_apc_conservation",
        "w_apc_protein_function",
        "w_apc_local_nd",
        "w_apc_mutation_density",
        "w_apc_tf",
    ];

    /// Display names matching STAARpipeline R output.
    #[cfg(test)]
    pub const DISPLAY_NAMES: [&str; 11] = [
        "cadd_phred",
        "linsight",
        "fathmm_xf",
        "apc_epigenetics_active",
        "apc_epigenetics_repressed",
        "apc_epigenetics_transcription",
        "apc_conservation",
        "apc_protein_function",
        "apc_local_nucleotide_diversity",
        "apc_mutation_density",
        "apc_transcription_factor",
    ];
}

// ---------------------------------------------------------------------------
// RegionType
// ---------------------------------------------------------------------------

/// Strongly-typed GENCODE region type. Eliminates per-variant String allocations
/// and replaces `.contains()` string matching in mask predicates with integer
/// comparisons.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum RegionType {
    Exonic,
    Intronic,
    Intergenic,
    Splicing,
    Upstream,
    Downstream,
    UpstreamDownstream,
    Utr3,
    Utr5,
    NcrnaExonic,
    NcrnaIntronic,
    NcrnaSplicing,
    NcrnaExonicSplicing,
    Utr5Utr3,
    Unknown,
}

impl RegionType {
    pub fn contains_splicing(self) -> bool {
        matches!(
            self,
            Self::Splicing | Self::NcrnaSplicing | Self::NcrnaExonicSplicing
        )
    }
    pub fn contains_upstream(self) -> bool {
        matches!(self, Self::Upstream | Self::UpstreamDownstream)
    }
    pub fn contains_downstream(self) -> bool {
        matches!(self, Self::Downstream | Self::UpstreamDownstream)
    }
    pub fn contains_utr(self) -> bool {
        matches!(self, Self::Utr3 | Self::Utr5 | Self::Utr5Utr3)
    }
    pub fn contains_ncrna(self) -> bool {
        matches!(
            self,
            Self::NcrnaExonic
                | Self::NcrnaIntronic
                | Self::NcrnaSplicing
                | Self::NcrnaExonicSplicing
        )
    }
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Exonic => "exonic",
            Self::Intronic => "intronic",
            Self::Intergenic => "intergenic",
            Self::Splicing => "splicing",
            Self::Upstream => "upstream",
            Self::Downstream => "downstream",
            Self::UpstreamDownstream => "upstream;downstream",
            Self::Utr3 => "UTR3",
            Self::Utr5 => "UTR5",
            Self::NcrnaExonic => "ncRNA_exonic",
            Self::NcrnaIntronic => "ncRNA_intronic",
            Self::NcrnaSplicing => "ncRNA_splicing",
            Self::NcrnaExonicSplicing => "ncRNA_exonic;splicing",
            Self::Utr5Utr3 => "UTR5;UTR3",
            Self::Unknown => "",
        }
    }
    pub fn from_str_lossy(s: &str) -> Self {
        match s {
            "exonic" => Self::Exonic,
            "intronic" => Self::Intronic,
            "intergenic" => Self::Intergenic,
            "splicing" => Self::Splicing,
            "upstream" => Self::Upstream,
            "downstream" => Self::Downstream,
            "upstream;downstream" => Self::UpstreamDownstream,
            "UTR3" => Self::Utr3,
            "UTR5" => Self::Utr5,
            "ncRNA_exonic" => Self::NcrnaExonic,
            "ncRNA_intronic" => Self::NcrnaIntronic,
            "ncRNA_splicing" => Self::NcrnaSplicing,
            "ncRNA_exonic;splicing" => Self::NcrnaExonicSplicing,
            "UTR5;UTR3" => Self::Utr5Utr3,
            "" => Self::Unknown,
            other => {
                if other.contains("ncRNA") && other.contains("splicing") {
                    Self::NcrnaExonicSplicing
                } else if other.contains("UTR") {
                    Self::Utr5Utr3
                } else if other.contains("ncRNA") {
                    Self::NcrnaExonic
                } else if other.contains("upstream") && other.contains("downstream") {
                    Self::UpstreamDownstream
                } else if other.contains("upstream") {
                    Self::Upstream
                } else if other.contains("downstream") {
                    Self::Downstream
                } else if other.contains("splicing") {
                    Self::Splicing
                } else {
                    Self::Unknown
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Consequence
// ---------------------------------------------------------------------------

/// Strongly-typed variant consequence. Eliminates per-variant String allocations
/// and turns mask predicates into integer match arms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum Consequence {
    Stopgain,
    StopGained,
    Stoploss,
    StopLost,
    FrameshiftInsertion,
    FrameshiftDeletion,
    FrameshiftSubstitution,
    FrameshiftVariant,
    Splicing,
    SpliceDonorVariant,
    SpliceAcceptorVariant,
    StartLost,
    NonsynonymousSNV,
    MissenseVariant,
    SynonymousSNV,
    SynonymousVariant,
    UpstreamGeneVariant,
    DownstreamGeneVariant,
    FivePrimeUtrVariant,
    ThreePrimeUtrVariant,
    NonCodingTranscriptExonVariant,
    NonCodingTranscriptVariant,
    Unknown,
}

impl Consequence {
    pub fn is_plof(self) -> bool {
        matches!(
            self,
            Self::Stopgain
                | Self::StopGained
                | Self::Stoploss
                | Self::StopLost
                | Self::FrameshiftInsertion
                | Self::FrameshiftDeletion
                | Self::FrameshiftSubstitution
                | Self::FrameshiftVariant
                | Self::Splicing
                | Self::SpliceDonorVariant
                | Self::SpliceAcceptorVariant
                | Self::StartLost
        )
    }
    pub fn is_ptv(self) -> bool {
        matches!(
            self,
            Self::Stopgain
                | Self::StopGained
                | Self::FrameshiftInsertion
                | Self::FrameshiftDeletion
                | Self::FrameshiftSubstitution
                | Self::FrameshiftVariant
        )
    }
    pub fn is_splice(self) -> bool {
        matches!(
            self,
            Self::Splicing | Self::SpliceDonorVariant | Self::SpliceAcceptorVariant
        )
    }
    pub fn is_missense(self) -> bool {
        matches!(self, Self::MissenseVariant | Self::NonsynonymousSNV)
    }
    pub fn is_synonymous(self) -> bool {
        matches!(self, Self::SynonymousSNV | Self::SynonymousVariant)
    }
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Stopgain => "stopgain",
            Self::StopGained => "stop_gained",
            Self::Stoploss => "stoploss",
            Self::StopLost => "stop_lost",
            Self::FrameshiftInsertion => "frameshift insertion",
            Self::FrameshiftDeletion => "frameshift deletion",
            Self::FrameshiftSubstitution => "frameshift substitution",
            Self::FrameshiftVariant => "frameshift_variant",
            Self::Splicing => "splicing",
            Self::SpliceDonorVariant => "splice_donor_variant",
            Self::SpliceAcceptorVariant => "splice_acceptor_variant",
            Self::StartLost => "start_lost",
            Self::NonsynonymousSNV => "nonsynonymous SNV",
            Self::MissenseVariant => "missense_variant",
            Self::SynonymousSNV => "synonymous SNV",
            Self::SynonymousVariant => "synonymous_variant",
            Self::UpstreamGeneVariant => "upstream_gene_variant",
            Self::DownstreamGeneVariant => "downstream_gene_variant",
            Self::FivePrimeUtrVariant => "5_prime_UTR_variant",
            Self::ThreePrimeUtrVariant => "3_prime_UTR_variant",
            Self::NonCodingTranscriptExonVariant => "non_coding_transcript_exon_variant",
            Self::NonCodingTranscriptVariant => "non_coding_transcript_variant",
            Self::Unknown => "",
        }
    }
    pub fn from_str_lossy(s: &str) -> Self {
        match s {
            "stopgain" => Self::Stopgain,
            "stop_gained" => Self::StopGained,
            "stoploss" => Self::Stoploss,
            "stop_lost" => Self::StopLost,
            "frameshift insertion" => Self::FrameshiftInsertion,
            "frameshift deletion" => Self::FrameshiftDeletion,
            "frameshift substitution" => Self::FrameshiftSubstitution,
            "frameshift_variant" => Self::FrameshiftVariant,
            "splicing" => Self::Splicing,
            "splice_donor_variant" => Self::SpliceDonorVariant,
            "splice_acceptor_variant" => Self::SpliceAcceptorVariant,
            "start_lost" => Self::StartLost,
            "nonsynonymous SNV" => Self::NonsynonymousSNV,
            "missense_variant" => Self::MissenseVariant,
            "synonymous SNV" => Self::SynonymousSNV,
            "synonymous_variant" => Self::SynonymousVariant,
            "upstream_gene_variant" => Self::UpstreamGeneVariant,
            "downstream_gene_variant" => Self::DownstreamGeneVariant,
            "5_prime_UTR_variant" => Self::FivePrimeUtrVariant,
            "3_prime_UTR_variant" => Self::ThreePrimeUtrVariant,
            "non_coding_transcript_exon_variant" => Self::NonCodingTranscriptExonVariant,
            "non_coding_transcript_variant" => Self::NonCodingTranscriptVariant,
            "" => Self::Unknown,
            _ => Self::Unknown,
        }
    }
}

// ---------------------------------------------------------------------------
// RegulatoryFlags
// ---------------------------------------------------------------------------

/// Regulatory region flags used by mask predicates.
///
/// Replaces 4 loose `bool` fields repeated in the old AnnotatedVariant and
/// MetaVariant structs.
#[derive(Clone, Copy, Debug, Default)]
pub struct RegulatoryFlags {
    pub cage_promoter: bool,
    pub cage_enhancer: bool,
    pub ccre_promoter: bool,
    pub ccre_enhancer: bool,
}

// ---------------------------------------------------------------------------
// FunctionalAnnotation
// ---------------------------------------------------------------------------

/// Everything STAAR needs to classify (masks) and weight (score tests) a variant.
///
/// Mask predicates access `consequence` / `region_type`.
/// Score tests access `weights`.
#[derive(Clone, Copy, Debug)]
pub struct FunctionalAnnotation {
    pub region_type: RegionType,
    pub consequence: Consequence,
    pub cadd_phred: f64,
    pub revel: f64,
    pub regulatory: RegulatoryFlags,
    pub weights: AnnotationWeights,
}

// ---------------------------------------------------------------------------
// AnnotatedVariant (canonical)
// ---------------------------------------------------------------------------

/// Single canonical variant type used by masks, score tests, sumstats export,
/// result writing, and MetaSTAAR merge.
///
/// Replaces the parallel structs in `masks.rs` and `meta.rs`.
#[derive(Clone, Debug)]
pub struct AnnotatedVariant {
    pub chromosome: Chromosome,
    pub position: u32,
    pub ref_allele: Box<str>,
    pub alt_allele: Box<str>,
    pub maf: f64,
    pub gene_name: Box<str>,
    pub annotation: FunctionalAnnotation,
}

// ---------------------------------------------------------------------------
// Variant ID
// ---------------------------------------------------------------------------

/// Format a universal variant ID: "chr-pos-ref-alt" (e.g., "1-1001-A-T").
/// Single source of truth — all vid construction calls this.
pub fn format_vid(chrom: &str, position: u32, ref_allele: &str, alt_allele: &str) -> String {
    format!("{chrom}-{position}-{ref_allele}-{alt_allele}")
}

// ---------------------------------------------------------------------------
// MetaVariant
// ---------------------------------------------------------------------------

/// A variant in a meta-analysis context. Composes `AnnotatedVariant` instead
/// of duplicating its fields.
///
/// Where masks need `&AnnotatedVariant`, pass `&mv.variant`. Zero cost,
/// zero duplication. Deletes the need for `meta_to_annotated()`.
pub struct MetaVariant {
    pub variant: AnnotatedVariant,
    pub u_meta: f64,
    pub mac_total: i64,
    pub n_total: i64,
    pub study_segments: Vec<(usize, i32)>,
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // -- Chromosome --------------------------------------------------------

    #[test]
    fn parse_autosomes() {
        assert_eq!(
            "chr1".parse::<Chromosome>().unwrap(),
            Chromosome::Autosome(1)
        );
        assert_eq!(
            "22".parse::<Chromosome>().unwrap(),
            Chromosome::Autosome(22)
        );
        assert_eq!(
            "chr10".parse::<Chromosome>().unwrap(),
            Chromosome::Autosome(10)
        );
    }

    #[test]
    fn parse_sex_and_mt() {
        assert_eq!("chrX".parse::<Chromosome>().unwrap(), Chromosome::X);
        assert_eq!("Y".parse::<Chromosome>().unwrap(), Chromosome::Y);
        assert_eq!("chrM".parse::<Chromosome>().unwrap(), Chromosome::MT);
        assert_eq!("MT".parse::<Chromosome>().unwrap(), Chromosome::MT);
        assert_eq!("chrMT".parse::<Chromosome>().unwrap(), Chromosome::MT);
    }

    #[test]
    fn parse_invalid() {
        assert!("chr0".parse::<Chromosome>().is_err());
        assert!("chr23".parse::<Chromosome>().is_err());
        assert!("banana".parse::<Chromosome>().is_err());
    }

    #[test]
    fn display_roundtrip() {
        for s in ["chr1", "chr22", "chrX", "chrY", "chrMT"] {
            let c: Chromosome = s.parse().unwrap();
            assert_eq!(c.to_string(), s);
        }
    }

    #[test]
    fn ordering() {
        let mut chroms: Vec<Chromosome> = vec![
            Chromosome::MT,
            Chromosome::Autosome(22),
            Chromosome::X,
            Chromosome::Autosome(1),
            Chromosome::Y,
            Chromosome::Autosome(3),
        ];
        chroms.sort();
        assert_eq!(
            chroms,
            vec![
                Chromosome::Autosome(1),
                Chromosome::Autosome(3),
                Chromosome::Autosome(22),
                Chromosome::X,
                Chromosome::Y,
                Chromosome::MT,
            ]
        );
    }

    #[test]
    fn serde_roundtrip() {
        let c = Chromosome::Autosome(7);
        let json = serde_json::to_string(&c).unwrap();
        assert_eq!(json, "\"chr7\"");
        let back: Chromosome = serde_json::from_str(&json).unwrap();
        assert_eq!(back, c);
    }

    // -- AnnotationWeights -------------------------------------------------

    #[test]
    fn channel_counts() {
        assert_eq!(AnnotationWeights::NAMES.len(), 11);
        assert_eq!(AnnotationWeights::DISPLAY_NAMES.len(), 11);
    }

    // -- Display names match existing ANNOTATION_CHANNELS -------------------

    #[test]
    fn display_names_match_annotation_channels() {
        // These must stay in sync with weights.rs::ANNOTATION_CHANNELS
        let expected = [
            "cadd_phred",
            "linsight",
            "fathmm_xf",
            "apc_epigenetics_active",
            "apc_epigenetics_repressed",
            "apc_epigenetics_transcription",
            "apc_conservation",
            "apc_protein_function",
            "apc_local_nucleotide_diversity",
            "apc_mutation_density",
            "apc_transcription_factor",
        ];
        assert_eq!(AnnotationWeights::DISPLAY_NAMES, expected);
    }
}
