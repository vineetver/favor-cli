//! Canonical column identifiers for the FAVOR pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Col {
    Chromosome,
    Position,
    EndPosition,
    RefAllele,
    AltAllele,

    Maf,

    GeneName,
    RegionType,
    Consequence,

    CaddPhred,
    Revel,
    MetaSvmPred,
    GeneHancer,

    Vid,
    VariantVcf,

    IsCagePromoter,
    IsCageEnhancer,
    IsCcrePromoter,
    IsCcreEnhancer,

    Linsight,
    FathmmXf,
    ApcEpigeneticsActive,
    ApcEpigeneticsRepressed,
    ApcEpigeneticsTranscription,
    ApcConservation,
    ApcProteinFunction,
    ApcLocalNucleotideDiversity,
    ApcMutationDensity,
    ApcTranscriptionFactor,
}

impl Col {
    pub fn as_str(self) -> &'static str {
        match self {
            Col::Chromosome => "chromosome",
            Col::Position => "position",
            Col::EndPosition => "end_position",
            Col::RefAllele => "ref_allele",
            Col::AltAllele => "alt_allele",
            Col::Vid => "vid",
            Col::VariantVcf => "variant_vcf",
            Col::Maf => "maf",
            Col::GeneName => "gene_name",
            Col::RegionType => "region_type",
            Col::Consequence => "consequence",
            Col::CaddPhred => "cadd_phred",
            Col::Revel => "revel",
            Col::MetaSvmPred => "metasvm_pred",
            Col::GeneHancer => "genehancer",
            Col::IsCagePromoter => "is_cage_promoter",
            Col::IsCageEnhancer => "is_cage_enhancer",
            Col::IsCcrePromoter => "is_ccre_promoter",
            Col::IsCcreEnhancer => "is_ccre_enhancer",
            Col::Linsight => "linsight",
            Col::FathmmXf => "fathmm_xf",
            Col::ApcEpigeneticsActive => "apc_epigenetics_active",
            Col::ApcEpigeneticsRepressed => "apc_epigenetics_repressed",
            Col::ApcEpigeneticsTranscription => "apc_epigenetics_transcription",
            Col::ApcConservation => "apc_conservation",
            Col::ApcProteinFunction => "apc_protein_function",
            Col::ApcLocalNucleotideDiversity => "apc_local_nucleotide_diversity",
            Col::ApcMutationDensity => "apc_mutation_density",
            Col::ApcTranscriptionFactor => "apc_transcription_factor",
        }
    }

    pub fn weight_index(self) -> Option<usize> {
        STAAR_PHRED_CHANNELS.iter().position(|&c| c == self)
    }

    pub fn annotation_name(self) -> Option<&'static str> {
        match self {
            Col::Position => Some("position"),
            Col::RefAllele => Some("ref_vcf"),
            Col::AltAllele => Some("alt_vcf"),
            Col::Vid => Some("vid"),
            Col::Chromosome => Some("chromosome"),
            _ => None,
        }
    }
}

impl std::fmt::Display for Col {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// FAVOR-native PHRED-scaled annotation channels used by STAAR. The score
/// kernel receives one transformed weight per channel per variant; the
/// transform `w = 1 - 10^(-phred/10)` is applied uniformly at load time
/// (see VariantIndex::load), turning PHRED into the [0,1] probability the
/// STAAR tests are parametrised in.
pub const STAAR_PHRED_CHANNELS: [Col; 11] = [
    Col::CaddPhred,
    Col::Linsight,
    Col::FathmmXf,
    Col::ApcEpigeneticsActive,
    Col::ApcEpigeneticsRepressed,
    Col::ApcEpigeneticsTranscription,
    Col::ApcConservation,
    Col::ApcProteinFunction,
    Col::ApcLocalNucleotideDiversity,
    Col::ApcMutationDensity,
    Col::ApcTranscriptionFactor,
];

/// Output column name paired with the DuckDB SQL that extracts the PHRED
/// value from the FAVOR annotation alias `a`. Emitted by
/// `annotation_join_sql` as `<sql> AS <name>`. The PHRED -> probability
/// transform happens later, inside the STAAR score kernel.
pub const STAAR_PHRED_SOURCES: [(&str, &str); 11] = [
    ("cadd_phred", "CAST(a.main.cadd.phred AS DOUBLE)"),
    ("linsight", "CAST(a.linsight AS DOUBLE)"),
    ("fathmm_xf", "CAST(a.fathmm_xf AS DOUBLE)"),
    (
        "apc_epigenetics_active",
        "CAST(a.apc.epigenetics_active AS DOUBLE)",
    ),
    (
        "apc_epigenetics_repressed",
        "CAST(a.apc.epigenetics_repressed AS DOUBLE)",
    ),
    (
        "apc_epigenetics_transcription",
        "CAST(a.apc.epigenetics_transcription AS DOUBLE)",
    ),
    (
        "apc_conservation",
        "CAST(a.apc.conservation_v2 AS DOUBLE)",
    ),
    (
        "apc_protein_function",
        "CAST(a.apc.protein_function_v3 AS DOUBLE)",
    ),
    (
        "apc_local_nucleotide_diversity",
        "CAST(a.apc.local_nucleotide_diversity_v3 AS DOUBLE)",
    ),
    (
        "apc_mutation_density",
        "CAST(a.apc.mutation_density AS DOUBLE)",
    ),
    (
        "apc_transcription_factor",
        "CAST(a.apc.transcription_factor AS DOUBLE)",
    ),
];

pub const STORE_METADATA_COLS: [Col; 15] = [
    Col::Position,
    Col::EndPosition,
    Col::RefAllele,
    Col::AltAllele,
    Col::Maf,
    Col::GeneName,
    Col::RegionType,
    Col::Consequence,
    Col::Revel,
    Col::MetaSvmPred,
    Col::GeneHancer,
    Col::IsCagePromoter,
    Col::IsCageEnhancer,
    Col::IsCcrePromoter,
    Col::IsCcreEnhancer,
];

/// Structural annotation columns beyond the 11 weight channels. Required
/// for mask predicates that match STAARpipeline (disruptive_missense,
/// plof_ds, ptv_ds via MetaSVM; GeneHancer pass-through for downstream
/// tooling). Locked at cohort preflight.
pub const STAAR_STRUCTURAL_COLS: [Col; 2] = [Col::MetaSvmPred, Col::GeneHancer];

pub fn store_columns() -> Vec<Col> {
    let mut cols = STORE_METADATA_COLS.to_vec();
    cols.extend_from_slice(&STAAR_PHRED_CHANNELS);
    cols
}

pub struct AnnotationExtract {
    pub output: Col,
    pub sql: &'static str,
}

pub static ANNOTATION_EXTRACTS: &[AnnotationExtract] = &[
    // 1-based inclusive end. Derived from the annotation table because
    // FAVOR is the authoritative variant identity in the join. When
    // FAVOR adds a native end_position field, swap this expression for a
    // direct column reference — output column name and downstream code
    // do not change.
    //
    // Explicit `CAST(... AS INTEGER)` because DataFusion 53 widens
    // `INTEGER + LENGTH(VARCHAR)` to BIGINT (LENGTH returns BIGINT), and
    // the cohort store builder reads this column as `Int32Array`. The
    // cast keeps the on-disk schema and the in-memory downcast aligned.
    AnnotationExtract {
        output: Col::EndPosition,
        sql: "CAST((a.position + LENGTH(a.ref_vcf) - 1) AS INTEGER)",
    },
    AnnotationExtract {
        output: Col::GeneName,
        sql: "COALESCE(a.gencode.genes[1], '')",
    },
    AnnotationExtract {
        output: Col::RegionType,
        sql: "COALESCE(a.gencode.region_type, '')",
    },
    AnnotationExtract {
        output: Col::Consequence,
        sql: "COALESCE(a.gencode.consequence, '')",
    },
    // REVEL is FLOAT (Float32) in FAVOR; the cohort store builder reads
    // it back as Float64Array so wrap COALESCE in CAST. CaddPhred is now
    // emitted via STAAR_PHRED_SOURCES alongside the other 10 channels.
    AnnotationExtract {
        output: Col::Revel,
        sql: "CAST(COALESCE(a.dbnsfp.revel, 0.0) AS DOUBLE)",
    },
    AnnotationExtract {
        output: Col::MetaSvmPred,
        sql: "COALESCE(a.dbnsfp.metasvm_pred, '')",
    },
    AnnotationExtract {
        output: Col::GeneHancer,
        sql: "COALESCE(a.genehancer.id, '')",
    },
    AnnotationExtract {
        output: Col::IsCagePromoter,
        sql: "a.cage.cage_promoter IS NOT NULL",
    },
    AnnotationExtract {
        output: Col::IsCageEnhancer,
        sql: "a.cage.cage_enhancer IS NOT NULL",
    },
    AnnotationExtract {
        output: Col::IsCcrePromoter,
        sql: "COALESCE(CAST(a.ccre.annotations AS VARCHAR) LIKE '%PLS%', false)",
    },
    AnnotationExtract {
        output: Col::IsCcreEnhancer,
        sql: "COALESCE(CAST(a.ccre.annotations AS VARCHAR) LIKE '%ELS%', false)",
    },
];

pub fn annotation_join_sql() -> String {
    let mut select_parts: Vec<String> = Vec::new();

    select_parts.push(format!(
        "CAST(g.chromosome AS VARCHAR) AS {}",
        Col::Chromosome
    ));
    select_parts.push(format!("g.position AS {}", Col::Position));
    select_parts.push(format!("g.\"ref\" AS {}", Col::RefAllele));
    select_parts.push(format!("g.alt AS {}", Col::AltAllele));
    // The genotype extractor writes `maf` as Float32; the cohort store
    // builder downcasts the column as `Float64Array`. Explicit cast at
    // the join boundary keeps both halves of the contract aligned.
    select_parts.push(format!("CAST(g.maf AS DOUBLE) AS {}", Col::Maf));

    for extract in ANNOTATION_EXTRACTS {
        select_parts.push(format!("{} AS {}", extract.sql, extract.output));
    }

    for (name, sql) in STAAR_PHRED_SOURCES {
        select_parts.push(format!("{sql} AS {name}"));
    }

    format!(
        "CREATE TABLE _rare_all AS\n\
         SELECT\n    {}\n\
         FROM _genotypes g\n\
         INNER JOIN _annotations a\n\
             ON g.chromosome = a.chromosome\n\
             AND g.position = a.position AND g.\"ref\" = a.ref_vcf AND g.alt = a.alt_vcf\n\
         WHERE g.maf > 0",
        select_parts.join(",\n    ")
    )
}

pub fn metadata_select_sql(chrom: &str) -> String {
    let cols: Vec<&str> = store_columns().iter().map(|c| c.as_str()).collect();
    format!(
        "SELECT {} FROM _rare_all WHERE {} = '{}' ORDER BY {}, {}, {}",
        cols.join(", "),
        Col::Chromosome,
        chrom,
        Col::Position,
        Col::RefAllele,
        Col::AltAllele,
    )
}

pub fn carrier_metadata_sql(chrom: &str) -> String {
    format!(
        "SELECT {pos}, {ref_a}, {alt_a}, {gene} \
         FROM _rare_all WHERE {chr} = '{chrom}' \
         ORDER BY {pos}, {ref_a}, {alt_a}",
        pos = Col::Position,
        ref_a = Col::RefAllele,
        alt_a = Col::AltAllele,
        gene = Col::GeneName,
        chr = Col::Chromosome,
    )
}

/// Deduplication query: count duplicate (chromosome, position, ref, alt) tuples.
pub fn dedup_count_sql() -> String {
    format!(
        "SELECT COUNT(*) FROM (\
            SELECT {chr}, {pos}, {ref_a}, {alt_a} \
            FROM _rare_all \
            GROUP BY {chr}, {pos}, {ref_a}, {alt_a} \
            HAVING COUNT(*) > 1\
        )",
        chr = Col::Chromosome,
        pos = Col::Position,
        ref_a = Col::RefAllele,
        alt_a = Col::AltAllele,
    )
}

/// Deduplication query: keep first occurrence (highest MAF) per variant.
pub fn dedup_sql() -> String {
    format!(
        "CREATE TABLE _rare_dedup AS \
         SELECT * FROM (\
             SELECT *, ROW_NUMBER() OVER (\
                 PARTITION BY {chr}, {pos}, {ref_a}, {alt_a} \
                 ORDER BY {maf} DESC\
             ) AS _rn FROM _rare_all\
         ) WHERE _rn = 1",
        chr = Col::Chromosome,
        pos = Col::Position,
        ref_a = Col::RefAllele,
        alt_a = Col::AltAllele,
        maf = Col::Maf,
    )
}

/// Query for distinct chromosomes in _rare_all.
pub fn distinct_chroms_sql() -> String {
    format!(
        "SELECT DISTINCT {} FROM _rare_all ORDER BY {}",
        Col::Chromosome,
        Col::Chromosome,
    )
}

/// Query for gene count in _rare_all.
pub fn gene_count_sql() -> String {
    format!(
        "SELECT COUNT(DISTINCT {gene}) FROM _rare_all WHERE {gene} != ''",
        gene = Col::GeneName,
    )
}

/// Simplified Arrow data type for compile-time schema contracts.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ColType {
    Int32,
    UInt32,
    Float32,
    Float64,
    Utf8,
    Boolean,
}

impl ColType {
    /// Check if an Arrow DataType is compatible with this contract type.
    /// Float64 accepts Float32 (widening). Utf8 accepts LargeUtf8 and Utf8View.
    pub fn matches(self, dt: &arrow::datatypes::DataType) -> bool {
        use arrow::datatypes::DataType;
        match self {
            ColType::Int32 => matches!(dt, DataType::Int32),
            ColType::UInt32 => matches!(dt, DataType::UInt32),
            ColType::Float32 => matches!(dt, DataType::Float32),
            ColType::Float64 => matches!(dt, DataType::Float64 | DataType::Float32),
            ColType::Utf8 => matches!(
                dt,
                DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View
            ),
            ColType::Boolean => matches!(dt, DataType::Boolean),
        }
    }
}

impl std::fmt::Display for ColType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ColType::Int32 => write!(f, "Int32"),
            ColType::UInt32 => write!(f, "UInt32"),
            ColType::Float32 => write!(f, "Float32"),
            ColType::Float64 => write!(f, "Float64"),
            ColType::Utf8 => write!(f, "Utf8"),
            ColType::Boolean => write!(f, "Boolean"),
        }
    }
}

/// A typed schema contract: a set of required (Col, ColType) pairs.
/// Validates actual Arrow schemas at pipeline boundaries with clear errors.
pub struct SchemaContract {
    pub name: &'static str,
    pub required: &'static [(Col, ColType)],
}

#[derive(Debug)]
pub enum SchemaError {
    Missing(Col),
    WrongType {
        col: Col,
        expected: ColType,
        actual: arrow::datatypes::DataType,
    },
}

impl std::fmt::Display for SchemaError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SchemaError::Missing(col) => {
                write!(f, "missing column '{}'", col.as_str())
            }
            SchemaError::WrongType {
                col,
                expected,
                actual,
            } => {
                write!(
                    f,
                    "column '{}': expected {expected}, got {actual:?}",
                    col.as_str()
                )
            }
        }
    }
}

impl SchemaContract {
    /// Validate an Arrow schema against this contract.
    pub fn validate(&self, schema: &arrow::datatypes::Schema) -> Result<(), Vec<SchemaError>> {
        let mut errors = Vec::new();
        for &(col, expected_type) in self.required {
            match schema.field_with_name(col.as_str()) {
                Ok(field) => {
                    if !expected_type.matches(field.data_type()) {
                        errors.push(SchemaError::WrongType {
                            col,
                            expected: expected_type,
                            actual: field.data_type().clone(),
                        });
                    }
                }
                Err(_) => errors.push(SchemaError::Missing(col)),
            }
        }
        if errors.is_empty() {
            Ok(())
        } else {
            Err(errors)
        }
    }

    /// Format validation errors into an actionable message.
    pub fn format_errors(errors: &[SchemaError]) -> String {
        errors
            .iter()
            .map(|e| format!("  {e}"))
            .collect::<Vec<_>>()
            .join("\n")
    }
}

/// Schema contract for the variant store parquet (variants.parquet).
pub static VARIANT_STORE_CONTRACT: SchemaContract = SchemaContract {
    name: "variant store",
    required: &[
        (Col::Position, ColType::Int32),
        (Col::EndPosition, ColType::Int32),
        (Col::RefAllele, ColType::Utf8),
        (Col::AltAllele, ColType::Utf8),
        (Col::Maf, ColType::Float64),
        (Col::GeneName, ColType::Utf8),
        (Col::RegionType, ColType::Utf8),
        (Col::Consequence, ColType::Utf8),
        (Col::CaddPhred, ColType::Float64),
        (Col::Revel, ColType::Float64),
        (Col::MetaSvmPred, ColType::Utf8),
        (Col::GeneHancer, ColType::Utf8),
        (Col::IsCagePromoter, ColType::Boolean),
        (Col::IsCageEnhancer, ColType::Boolean),
        (Col::IsCcrePromoter, ColType::Boolean),
        (Col::IsCcreEnhancer, ColType::Boolean),
        (Col::Linsight, ColType::Float64),
        (Col::FathmmXf, ColType::Float64),
        (Col::ApcEpigeneticsActive, ColType::Float64),
        (Col::ApcEpigeneticsRepressed, ColType::Float64),
        (Col::ApcEpigeneticsTranscription, ColType::Float64),
        (Col::ApcConservation, ColType::Float64),
        (Col::ApcProteinFunction, ColType::Float64),
        (Col::ApcLocalNucleotideDiversity, ColType::Float64),
        (Col::ApcMutationDensity, ColType::Float64),
        (Col::ApcTranscriptionFactor, ColType::Float64),
    ],
};

/// Schema contract for individual variant result parquet.
pub static INDIVIDUAL_RESULT_CONTRACT: SchemaContract = SchemaContract {
    name: "individual results",
    required: &[
        (Col::Chromosome, ColType::Utf8),
        (Col::Position, ColType::Int32),
        (Col::RefAllele, ColType::Utf8),
        (Col::AltAllele, ColType::Utf8),
        (Col::Maf, ColType::Float64),
        (Col::GeneName, ColType::Utf8),
        (Col::RegionType, ColType::Utf8),
        (Col::Consequence, ColType::Utf8),
        (Col::CaddPhred, ColType::Float64),
    ],
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn col_as_str_non_empty() {
        for &col in STAAR_PHRED_CHANNELS.iter().chain(STORE_METADATA_COLS.iter()) {
            let s = col.as_str();
            assert!(!s.is_empty(), "{col:?} has empty as_str()");
        }
    }

    #[test]
    fn staar_phred_channels_count() {
        assert_eq!(STAAR_PHRED_CHANNELS.len(), 11);
        assert_eq!(STAAR_PHRED_SOURCES.len(), 11);
    }

    #[test]
    fn phred_channel_names_match_sources() {
        for (i, col) in STAAR_PHRED_CHANNELS.iter().enumerate() {
            assert_eq!(
                col.as_str(),
                STAAR_PHRED_SOURCES[i].0,
                "channel {i}: Col {col:?} as_str() != STAAR_PHRED_SOURCES name"
            );
        }
    }

    #[test]
    fn weight_index_roundtrip() {
        for (i, col) in STAAR_PHRED_CHANNELS.iter().enumerate() {
            assert_eq!(
                col.weight_index(),
                Some(i),
                "{col:?} weight_index() mismatch"
            );
        }
        assert_eq!(Col::Chromosome.weight_index(), None);
    }

    #[test]
    fn store_columns_complete() {
        let cols = store_columns();
        assert_eq!(cols.len(), 26); // 15 metadata + 11 phred channels
    }

    #[test]
    fn annotation_join_sql_contains_all_columns() {
        let sql = annotation_join_sql();
        for col in &STAAR_PHRED_CHANNELS {
            assert!(
                sql.contains(col.as_str()),
                "annotation_join_sql() missing phred channel column {}",
                col.as_str()
            );
        }
        for extract in ANNOTATION_EXTRACTS {
            assert!(
                sql.contains(extract.output.as_str()),
                "annotation_join_sql() missing extraction column {}",
                extract.output.as_str()
            );
        }
    }

    #[test]
    fn annotation_join_sql_structure() {
        let sql = annotation_join_sql();
        assert!(sql.contains("CREATE TABLE _rare_all AS"));
        assert!(sql.contains("INNER JOIN"));
        assert!(sql.contains("g.chromosome = a.chromosome"));
        assert!(sql.contains("g.position = a.position"));
        assert!(sql.contains("g.maf > 0"));
        // Raw PHRED columns, no transform at SQL level.
        assert!(sql.contains("a.main.cadd.phred"));
        assert!(sql.contains("a.linsight"));
        assert!(sql.contains("a.apc.conservation_v2"));
        assert!(!sql.contains("POW(10.0"));
    }

    #[test]
    fn variant_store_contract_covers_all_store_columns() {
        let contract_cols: Vec<Col> = VARIANT_STORE_CONTRACT
            .required
            .iter()
            .map(|&(c, _)| c)
            .collect();
        for col in store_columns() {
            assert!(
                contract_cols.contains(&col),
                "VARIANT_STORE_CONTRACT missing {col:?}"
            );
        }
    }

    #[test]
    fn col_type_matches_arrow() {
        use arrow::datatypes::DataType;
        assert!(ColType::Int32.matches(&DataType::Int32));
        assert!(!ColType::Int32.matches(&DataType::Float64));
        assert!(ColType::Float64.matches(&DataType::Float64));
        assert!(ColType::Float64.matches(&DataType::Float32)); // widening
        assert!(ColType::Utf8.matches(&DataType::Utf8));
        assert!(ColType::Utf8.matches(&DataType::LargeUtf8));
        assert!(ColType::Boolean.matches(&DataType::Boolean));
        assert!(!ColType::Boolean.matches(&DataType::Int32));
    }

    #[test]
    fn metadata_select_sql_structure() {
        let sql = metadata_select_sql("22");
        assert!(sql.contains("FROM _rare_all"));
        assert!(sql.contains("WHERE chromosome = '22'"));
        assert!(sql.contains("ORDER BY position"));
        // All store columns present
        for col in store_columns() {
            assert!(
                sql.contains(col.as_str()),
                "metadata_select_sql missing {col:?}"
            );
        }
    }

    #[test]
    fn carrier_metadata_sql_structure() {
        let sql = carrier_metadata_sql("X");
        assert!(sql.contains(Col::Position.as_str()));
        assert!(sql.contains(Col::GeneName.as_str()));
        assert!(sql.contains("WHERE chromosome = 'X'"));
        assert!(sql.contains("ORDER BY position"));
    }

    #[test]
    fn dedup_sql_structure() {
        let sql = dedup_count_sql();
        assert!(sql.contains("GROUP BY"));
        assert!(sql.contains(Col::Chromosome.as_str()));
        assert!(sql.contains(Col::Position.as_str()));
        assert!(sql.contains("HAVING COUNT(*) > 1"));
    }

    #[test]
    fn no_duplicate_col_names() {
        let all_cols = store_columns();
        let mut names: Vec<&str> = all_cols.iter().map(|c| c.as_str()).collect();
        let len_before = names.len();
        names.sort();
        names.dedup();
        assert_eq!(
            names.len(),
            len_before,
            "Duplicate column names in store_columns()"
        );
    }
}
