//! Canonical object model for FAVOR CLI.

use std::fmt;
use std::str::FromStr;

use serde::{Deserialize, Serialize};

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

/// Filterable subset of the 25 canonical chromosomes (1-22, X, Y, MT).
/// Parsed from CLI specs like `"22"`, `"1,2,3"`, `"1-22"`, `"1-22,X,Y,MT"`.
/// Named chromosomes (X, Y, MT) may only appear as single tokens, never in
/// a range — `"X-Y"` is an error.
#[derive(Default, Clone, Debug, PartialEq, Eq)]
pub struct ChromosomeSet {
    bits: [bool; 25],
}

impl ChromosomeSet {
    fn index(c: Chromosome) -> usize {
        match c {
            Chromosome::Autosome(n) => (n - 1) as usize,
            Chromosome::X => 22,
            Chromosome::Y => 23,
            Chromosome::MT => 24,
        }
    }

    pub fn insert(&mut self, c: Chromosome) {
        self.bits[Self::index(c)] = true;
    }

    /// Fast membership check against the canonical strings that
    /// `ingest::vcf::normalize_chrom` emits (`"1"`..`"22"`, `"X"`, `"Y"`, `"MT"`).
    pub fn contains_canonical(&self, s: &str) -> bool {
        match s {
            "X" => self.bits[22],
            "Y" => self.bits[23],
            "MT" => self.bits[24],
            _ => s.parse::<u8>().ok()
                .filter(|n| (1..=22).contains(n))
                .map(|n| self.bits[(n - 1) as usize])
                .unwrap_or(false),
        }
    }

    pub fn is_empty(&self) -> bool {
        !self.bits.iter().any(|b| *b)
    }
}

impl FromStr for ChromosomeSet {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut set = ChromosomeSet::default();
        for token in s.split(',').map(str::trim).filter(|t| !t.is_empty()) {
            if let Some((a, b)) = token.split_once('-') {
                let start: u8 = a.trim().parse()
                    .map_err(|_| format!("range must use autosomes 1-22: '{token}'"))?;
                let end: u8 = b.trim().parse()
                    .map_err(|_| format!("range must use autosomes 1-22: '{token}'"))?;
                if !(1..=22).contains(&start) || !(1..=22).contains(&end) {
                    return Err(format!("range {token} outside autosomes 1-22"));
                }
                if start > end {
                    return Err(format!("reversed range: '{token}'"));
                }
                for n in start..=end {
                    set.insert(Chromosome::Autosome(n));
                }
            } else {
                let c: Chromosome = token.parse()?;
                set.insert(c);
            }
        }
        if set.is_empty() {
            return Err("no chromosomes specified".into());
        }
        Ok(set)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct AnnotationWeights(pub [f64; 11]);

impl Default for AnnotationWeights {
    fn default() -> Self {
        Self([0.0; 11])
    }
}

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

/// Format a universal variant ID: "chr-pos-ref-alt" (e.g., "1-1001-A-T").
/// Single source of truth — all vid construction calls this.
pub fn format_vid(chrom: &str, position: u32, ref_allele: &str, alt_allele: &str) -> String {
    format!("{chrom}-{position}-{ref_allele}-{alt_allele}")
}

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

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn chromset_single_token() {
        let s: ChromosomeSet = "22".parse().unwrap();
        assert!(s.contains_canonical("22"));
        assert!(!s.contains_canonical("21"));
        assert!(!s.contains_canonical("X"));
    }

    #[test]
    fn chromset_comma_list() {
        let s: ChromosomeSet = "1,3,X,MT".parse().unwrap();
        for c in ["1", "3", "X", "MT"] {
            assert!(s.contains_canonical(c), "missing {c}");
        }
        assert!(!s.contains_canonical("2"));
        assert!(!s.contains_canonical("Y"));
    }

    #[test]
    fn chromset_range_expands() {
        let s: ChromosomeSet = "1-5".parse().unwrap();
        for n in 1..=5u8 {
            assert!(s.contains_canonical(&n.to_string()), "missing {n}");
        }
        assert!(!s.contains_canonical("6"));
    }

    #[test]
    fn chromset_mixed() {
        let s: ChromosomeSet = "1-3,22,X,Y".parse().unwrap();
        for c in ["1", "2", "3", "22", "X", "Y"] {
            assert!(s.contains_canonical(c), "missing {c}");
        }
        assert!(!s.contains_canonical("4"));
        assert!(!s.contains_canonical("MT"));
    }

    #[test]
    fn chromset_contains_canonical() {
        let s: ChromosomeSet = "1-22,X,Y,MT".parse().unwrap();
        // Every canonical string normalize_chrom emits should hit.
        for n in 1..=22u8 {
            assert!(s.contains_canonical(&n.to_string()), "missing {n}");
        }
        for c in ["X", "Y", "MT"] {
            assert!(s.contains_canonical(c), "missing {c}");
        }
        // Uncanonical inputs don't match.
        assert!(!s.contains_canonical("chr1"));
        assert!(!s.contains_canonical("banana"));
        assert!(!s.contains_canonical("23"));
    }

    #[test]
    fn chromset_empty_spec_is_error() {
        assert!("".parse::<ChromosomeSet>().is_err());
        assert!(",".parse::<ChromosomeSet>().is_err());
    }

    #[test]
    fn chromset_reversed_range_is_error() {
        let e = "10-5".parse::<ChromosomeSet>().unwrap_err();
        assert!(e.contains("reversed"), "got: {e}");
    }

    #[test]
    fn chromset_named_in_range_is_error() {
        assert!("X-Y".parse::<ChromosomeSet>().is_err());
        assert!("1-X".parse::<ChromosomeSet>().is_err());
    }

    #[test]
    fn chromset_out_of_range_is_error() {
        assert!("0-5".parse::<ChromosomeSet>().is_err());
        assert!("20-25".parse::<ChromosomeSet>().is_err());
        assert!("99".parse::<ChromosomeSet>().is_err());
    }
}
