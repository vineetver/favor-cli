use std::path::Path;

use crate::types::{
    AnnotatedVariant, AnnotationWeights, Chromosome, Consequence, FunctionalAnnotation, RegionType,
    RegulatoryFlags,
};

// Ground truth annotation weights from real OR11H1 stopgain on chr22.
// Order matches AnnotationWeights::NAMES exactly.
const GROUND_TRUTH_WEIGHTS: [f64; 11] = [
    0.9957,  // cadd: 1 - 10^(-23.7/10)
    0.2149,  // linsight
    1.7630,  // fathmm_xf
    0.0126,  // apc_epigenetics_active
    0.3103,  // apc_epigenetics_repressed
    0.3234,  // apc_epigenetics_transcription
    6.3541,  // apc_conservation_v2
    2.9695,  // apc_protein_function_v3
    13.3931, // apc_local_nucleotide_diversity_v3
    13.3340, // apc_mutation_density
    3.1428,  // apc_transcription_factor
];

pub fn base_variant() -> AnnotatedVariant {
    AnnotatedVariant {
        chromosome: Chromosome::Autosome(22),
        position: 15528164,
        ref_allele: "C".into(),
        alt_allele: "T".into(),
        maf: 0.0007,
        gene_name: "OR11H1".into(),
        annotation: FunctionalAnnotation {
            region_type: RegionType::Exonic,
            consequence: Consequence::Stopgain,
            cadd_phred: 23.7,
            revel: 0.0,
            weights: AnnotationWeights(GROUND_TRUTH_WEIGHTS),
            regulatory: RegulatoryFlags::default(),
        },
    }
}

pub fn stopgain() -> AnnotatedVariant {
    base_variant()
}

pub fn splice() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528200,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Splicing,
            consequence: Consequence::Unknown,
            cadd_phred: 25.1,
            revel: 0.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn missense_high() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528300,
        annotation: FunctionalAnnotation {
            consequence: Consequence::NonsynonymousSNV,
            cadd_phred: 28.0,
            revel: 0.8,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn missense_low() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528400,
        annotation: FunctionalAnnotation {
            consequence: Consequence::NonsynonymousSNV,
            cadd_phred: 8.0,
            revel: 0.2,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn synonymous() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528500,
        annotation: FunctionalAnnotation {
            consequence: Consequence::SynonymousSNV,
            cadd_phred: 6.4,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn frameshift() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528600,
        annotation: FunctionalAnnotation {
            consequence: Consequence::FrameshiftDeletion,
            cadd_phred: 20.9,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn stoploss() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528700,
        annotation: FunctionalAnnotation {
            consequence: Consequence::Stoploss,
            cadd_phred: 15.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn upstream() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15527000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Upstream,
            consequence: Consequence::Unknown,
            cadd_phred: 5.6,
            revel: 0.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn downstream() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15530000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Downstream,
            consequence: Consequence::Unknown,
            cadd_phred: 4.2,
            revel: 0.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn upstream_downstream() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15531000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::UpstreamDownstream,
            consequence: Consequence::Unknown,
            cadd_phred: 3.0,
            revel: 0.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn utr3() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15529000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Utr3,
            consequence: Consequence::Unknown,
            cadd_phred: 7.4,
            revel: 0.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn cage_promoter() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15526000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Intergenic,
            consequence: Consequence::Unknown,
            cadd_phred: 12.3,
            revel: 0.0,
            regulatory: RegulatoryFlags {
                cage_promoter: true,
                ..RegulatoryFlags::default()
            },
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn cage_enhancer() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15525000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Intergenic,
            consequence: Consequence::Unknown,
            cadd_phred: 8.1,
            revel: 0.0,
            regulatory: RegulatoryFlags {
                cage_enhancer: true,
                ..RegulatoryFlags::default()
            },
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn ccre_pls() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15524000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Intergenic,
            consequence: Consequence::Unknown,
            cadd_phred: 7.0,
            revel: 0.0,
            regulatory: RegulatoryFlags {
                ccre_promoter: true,
                ..RegulatoryFlags::default()
            },
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn ccre_els() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15523000,
        annotation: FunctionalAnnotation {
            region_type: RegionType::Intergenic,
            consequence: Consequence::Unknown,
            cadd_phred: 8.7,
            revel: 0.0,
            regulatory: RegulatoryFlags {
                ccre_enhancer: true,
                ..RegulatoryFlags::default()
            },
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn ncrna() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15532000,
        gene_name: "RF00004".into(),
        annotation: FunctionalAnnotation {
            region_type: RegionType::NcrnaExonic,
            consequence: Consequence::Unknown,
            cadd_phred: 9.9,
            revel: 0.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn intronic() -> AnnotatedVariant {
    AnnotatedVariant {
        position: 15528800,
        gene_name: "POTEH".into(),
        annotation: FunctionalAnnotation {
            region_type: RegionType::Intronic,
            consequence: Consequence::Unknown,
            cadd_phred: 13.9,
            revel: 0.0,
            ..base_variant().annotation
        },
        ..base_variant()
    }
}

pub fn all_variants() -> Vec<AnnotatedVariant> {
    vec![
        stopgain(),
        splice(),
        missense_high(),
        missense_low(),
        synonymous(),
        frameshift(),
        stoploss(),
        upstream(),
        downstream(),
        upstream_downstream(),
        utr3(),
        cage_promoter(),
        cage_enhancer(),
        ccre_pls(),
        ccre_els(),
        ncrna(),
        intronic(),
    ]
}

/// SQL that creates a FAVOR-schema-correct annotation query using DataFusion syntax.
fn annotation_rows_sql() -> String {
    let rows: Vec<String> = all_variants().iter().map(|v| {
        let cage = if v.annotation.regulatory.cage_promoter && v.annotation.regulatory.cage_enhancer {
            "named_struct('cage_enhancer', 'e1@chr22', 'cage_promoter', 'p1@chr22', 'cage_tc', CAST(NULL AS VARCHAR))".to_string()
        } else if v.annotation.regulatory.cage_promoter {
            "named_struct('cage_enhancer', CAST(NULL AS VARCHAR), 'cage_promoter', 'p1@chr22', 'cage_tc', CAST(NULL AS VARCHAR))".to_string()
        } else if v.annotation.regulatory.cage_enhancer {
            "named_struct('cage_enhancer', 'e1@chr22', 'cage_promoter', CAST(NULL AS VARCHAR), 'cage_tc', CAST(NULL AS VARCHAR))".to_string()
        } else {
            "named_struct('cage_enhancer', CAST(NULL AS VARCHAR), 'cage_promoter', CAST(NULL AS VARCHAR), 'cage_tc', CAST(NULL AS VARCHAR))".to_string()
        };

        let ccre_ann = if v.annotation.regulatory.ccre_promoter { "'PLS'" }
            else if v.annotation.regulatory.ccre_enhancer { "'dELS'" }
            else { "CAST(NULL AS VARCHAR)" };

        let consequence_sql = if v.annotation.consequence == Consequence::Unknown { "CAST(NULL AS VARCHAR)".into() }
            else { format!("'{}'", v.annotation.consequence.as_str()) };

        let revel_sql = if v.annotation.revel == 0.0 { "CAST(NULL AS FLOAT)".into() }
            else { format!("CAST({} AS FLOAT)", v.annotation.revel) };

        format!(
            "SELECT \
                CAST(22 AS BIGINT) AS chromosome, \
                {pos} AS position, \
                'T' AS ref_vcf, \
                'A' AS alt_vcf, \
                named_struct('region_type', '{rt}', 'genes', make_array('{gene}'), \
                  'consequence', {csq}) AS gencode, \
                named_struct('cadd', named_struct('raw', CAST(0.0 AS FLOAT), 'phred', CAST({cadd} AS FLOAT))) AS main, \
                named_struct('revel', {revel}) AS dbnsfp, \
                CAST({linsight} AS FLOAT) AS linsight, \
                CAST({fathmm} AS FLOAT) AS fathmm_xf, \
                {cage} AS cage, \
                named_struct('ids', CAST(NULL AS VARCHAR), 'accessions', CAST(NULL AS VARCHAR), \
                  'annotations', {ccre_ann}, 'count', CAST(0 AS TINYINT)) AS ccre, \
                named_struct('conservation_v2', CAST({w6} AS FLOAT), 'epigenetics', CAST(0.0 AS FLOAT), \
                  'epigenetics_active', CAST({w3} AS FLOAT), 'epigenetics_repressed', CAST({w4} AS FLOAT), \
                  'epigenetics_transcription', CAST({w5} AS FLOAT), \
                  'local_nucleotide_diversity_v3', CAST({w8} AS FLOAT), \
                  'mappability', CAST(0.0 AS FLOAT), 'micro_rna', CAST(0.0 AS FLOAT), \
                  'mutation_density', CAST({w9} AS FLOAT), \
                  'protein_function_v3', CAST({w7} AS FLOAT), \
                  'proximity_to_coding_v2', CAST(0.0 AS FLOAT), \
                  'proximity_to_tsstes', CAST(0.0 AS FLOAT), \
                  'transcription_factor', CAST({w10} AS FLOAT)) AS apc",
            pos = v.position,
            rt = v.annotation.region_type.as_str(),
            gene = v.gene_name,
            csq = consequence_sql,
            cadd = v.annotation.cadd_phred,
            revel = revel_sql,
            linsight = v.annotation.weights.0[1],
            fathmm = v.annotation.weights.0[2],
            cage = cage,
            ccre_ann = ccre_ann,
            w3 = v.annotation.weights.0[3],
            w4 = v.annotation.weights.0[4],
            w5 = v.annotation.weights.0[5],
            w6 = v.annotation.weights.0[6],
            w7 = v.annotation.weights.0[7],
            w8 = v.annotation.weights.0[8],
            w9 = v.annotation.weights.0[9],
            w10 = v.annotation.weights.0[10],
        )
    }).collect();

    rows.join(" UNION ALL ")
}

pub fn write_test_annotation_parquet(dir: &Path) -> std::path::PathBuf {
    let out = dir.join("test_annotations.parquet");
    let resources = crate::resource::Resources::detect();
    let engine = crate::engine::DfEngine::new(&resources).unwrap();
    let sql = format!(
        "COPY ({}) TO '{}' STORED AS PARQUET OPTIONS (compression 'zstd(4)')",
        annotation_rows_sql(),
        out.display(),
    );
    engine.execute(&sql).unwrap();
    out
}

pub fn write_test_genotype_parquet(dir: &Path) -> std::path::PathBuf {
    let out = dir.join("test_genotypes.parquet");
    let resources = crate::resource::Resources::detect();
    let engine = crate::engine::DfEngine::new(&resources).unwrap();

    let rows: Vec<String> = all_variants().iter().map(|v| {
        format!(
            "SELECT '22' AS chromosome, {} AS position, 'T' AS ref, 'A' AS alt, \
             CAST({} AS FLOAT) AS maf, make_array(CAST(0.0 AS FLOAT), CAST(1.0 AS FLOAT), CAST(0.0 AS FLOAT), CAST(2.0 AS FLOAT), CAST(0.0 AS FLOAT), CAST(1.0 AS FLOAT), CAST(0.0 AS FLOAT), CAST(0.0 AS FLOAT), CAST(1.0 AS FLOAT), CAST(0.0 AS FLOAT)) AS dosages",
            v.position, v.maf,
        )
    }).collect();

    let sql = format!(
        "COPY ({}) TO '{}' STORED AS PARQUET OPTIONS (compression 'zstd(4)')",
        rows.join(" UNION ALL "),
        out.display(),
    );
    engine.execute(&sql).unwrap();
    out
}

/// The exact SQL STAAR uses to build the _rare table, parameterized for test table names.
/// Uses the same column naming as the production annotation_join_sql(), but with test
/// table names and an upper MAF bound for rare-variant filtering.
pub fn staar_rare_sql() -> String {
    use crate::column::Col;
    format!("CREATE TABLE _rare AS
     SELECT
         CAST(g.chromosome AS VARCHAR) AS {chr},
         g.position AS {pos},
         g.ref AS {ref_a},
         g.alt AS {alt_a},
         g.maf,
         COALESCE(array_element(a.gencode.genes, 1), '') AS {gene},
         COALESCE(a.gencode.region_type, '') AS {region},
         COALESCE(a.gencode.consequence, '') AS {csq},
         COALESCE(a.main.cadd.phred, 0) AS {cadd},
         COALESCE(a.dbnsfp.revel, 0) AS {revel},
         a.cage.cage_promoter IS NOT NULL AS {cage_p},
         a.cage.cage_enhancer IS NOT NULL AS {cage_e},
         COALESCE(CAST(a.ccre.annotations AS VARCHAR) LIKE '%PLS%', false) AS {ccre_p},
         COALESCE(CAST(a.ccre.annotations AS VARCHAR) LIKE '%ELS%', false) AS {ccre_e},
         CASE WHEN a.main.cadd.phred > 0 THEN 1.0 - POW(10.0, -a.main.cadd.phred / 10.0) ELSE 0.0 END AS {w_cadd},
         CAST(COALESCE(a.linsight, 0) AS DOUBLE) AS {w_lin},
         CAST(COALESCE(a.fathmm_xf, 0) AS DOUBLE) AS {w_fxf},
         CAST(COALESCE(a.apc.epigenetics_active, 0) AS DOUBLE) AS {w_ea},
         CAST(COALESCE(a.apc.epigenetics_repressed, 0) AS DOUBLE) AS {w_er},
         CAST(COALESCE(a.apc.epigenetics_transcription, 0) AS DOUBLE) AS {w_et},
         CAST(COALESCE(a.apc.conservation_v2, 0) AS DOUBLE) AS {w_con},
         CAST(COALESCE(a.apc.protein_function_v3, 0) AS DOUBLE) AS {w_pf},
         CAST(COALESCE(a.apc.local_nucleotide_diversity_v3, 0) AS DOUBLE) AS {w_lnd},
         CAST(COALESCE(a.apc.mutation_density, 0) AS DOUBLE) AS {w_md},
         CAST(COALESCE(a.apc.transcription_factor, 0) AS DOUBLE) AS {w_tf}
     FROM _test_geno g
     INNER JOIN _test_ann a
         ON CAST(g.chromosome AS VARCHAR) = CAST(a.chromosome AS VARCHAR)
         AND g.position = a.position AND g.ref = a.ref_vcf AND g.alt = a.alt_vcf
     WHERE g.maf > 0 AND g.maf < 0.01",
        chr = Col::Chromosome, pos = Col::Position,
        ref_a = Col::RefAllele, alt_a = Col::AltAllele,
        gene = Col::GeneName, region = Col::RegionType, csq = Col::Consequence,
        cadd = Col::CaddPhred, revel = Col::Revel,
        cage_p = Col::IsCagePromoter, cage_e = Col::IsCageEnhancer,
        ccre_p = Col::IsCcrePromoter, ccre_e = Col::IsCcreEnhancer,
        w_cadd = Col::WCadd, w_lin = Col::WLinsight, w_fxf = Col::WFathmmXf,
        w_ea = Col::WApcEpiActive, w_er = Col::WApcEpiRepressed,
        w_et = Col::WApcEpiTranscription, w_con = Col::WApcConservation,
        w_pf = Col::WApcProteinFunction, w_lnd = Col::WApcLocalNd,
        w_md = Col::WApcMutationDensity, w_tf = Col::WApcTf,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::DfEngine;
    use crate::resource::Resources;
    use crate::staar::masks;
    use crate::staar::MaskType;
    use arrow::array::{Array, BooleanArray, Float64Array, Int32Array, Int64Array, StringArray};

    // DataFusion 53 may return Utf8View instead of Utf8 for CAST(... AS VARCHAR).
    fn str_col_val(col: &dyn Array, row: usize) -> String {
        if col.is_null(row) {
            return String::new();
        }
        if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
            return a.value(row).to_string();
        }
        arrow::util::display::array_value_to_string(col, row).unwrap_or_default()
    }

    // DataFusion 53 may return different numeric types than expected.
    // Use array_value_to_string + parse as a universal extractor.
    fn num_col_val(col: &dyn Array, row: usize) -> f64 {
        if col.is_null(row) {
            return 0.0;
        }
        if let Some(a) = col.as_any().downcast_ref::<Float64Array>() {
            return a.value(row);
        }
        if let Some(a) = col.as_any().downcast_ref::<Int32Array>() {
            return a.value(row) as f64;
        }
        if let Some(a) = col.as_any().downcast_ref::<Int64Array>() {
            return a.value(row) as f64;
        }
        arrow::util::display::array_value_to_string(col, row)
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .unwrap_or(0.0)
    }

    fn bool_col_val(col: &dyn Array, row: usize) -> bool {
        if col.is_null(row) {
            return false;
        }
        if let Some(a) = col.as_any().downcast_ref::<BooleanArray>() {
            return a.value(row);
        }
        false
    }

    fn read_extracted_variants(engine: &DfEngine) -> Vec<AnnotatedVariant> {
        use crate::column::Col;
        // _rare table has: chromosome, position, ref_allele, alt_allele, maf,
        // gene_name, region_type, consequence, cadd_phred, revel,
        // is_cage_promoter..is_ccre_enhancer, w_cadd..w_apc_tf
        // (25 columns total = chromosome + 13 metadata + 11 weights)
        let batches = engine
            .collect(&format!(
                "SELECT {chr}, {pos}, {ref_a}, {alt_a}, {maf}, \
                 {gene}, {region}, {csq}, {cadd}, {revel}, \
                 {cage_p}, {cage_e}, {ccre_p}, {ccre_e}, \
                 {weights} \
                 FROM _rare ORDER BY {pos}",
                chr = Col::Chromosome,
                pos = Col::Position,
                ref_a = Col::RefAllele,
                alt_a = Col::AltAllele,
                maf = Col::Maf,
                gene = Col::GeneName,
                region = Col::RegionType,
                csq = Col::Consequence,
                cadd = Col::CaddPhred,
                revel = Col::Revel,
                cage_p = Col::IsCagePromoter,
                cage_e = Col::IsCageEnhancer,
                ccre_p = Col::IsCcrePromoter,
                ccre_e = Col::IsCcreEnhancer,
                weights = crate::column::STAAR_WEIGHTS
                    .iter()
                    .map(|c| c.as_str())
                    .collect::<Vec<_>>()
                    .join(", "),
            ))
            .unwrap();

        let mut variants = Vec::new();
        for batch in &batches {
            for row in 0..batch.num_rows() {
                let mut weights = [0.0f64; 11];
                for (i, w) in weights.iter_mut().enumerate() {
                    *w = num_col_val(batch.column(14 + i).as_ref(), row);
                }

                variants.push(AnnotatedVariant {
                    chromosome: str_col_val(batch.column(0).as_ref(), row)
                        .parse()
                        .unwrap_or(Chromosome::Autosome(1)),
                    position: num_col_val(batch.column(1).as_ref(), row) as u32,
                    ref_allele: str_col_val(batch.column(2).as_ref(), row).into(),
                    alt_allele: str_col_val(batch.column(3).as_ref(), row).into(),
                    maf: num_col_val(batch.column(4).as_ref(), row),
                    gene_name: str_col_val(batch.column(5).as_ref(), row).into(),
                    annotation: FunctionalAnnotation {
                        region_type: RegionType::from_str_lossy(&str_col_val(
                            batch.column(6).as_ref(),
                            row,
                        )),
                        consequence: Consequence::from_str_lossy(&str_col_val(
                            batch.column(7).as_ref(),
                            row,
                        )),
                        cadd_phred: num_col_val(batch.column(8).as_ref(), row),
                        revel: num_col_val(batch.column(9).as_ref(), row),
                        regulatory: RegulatoryFlags {
                            cage_promoter: bool_col_val(batch.column(10).as_ref(), row),
                            cage_enhancer: bool_col_val(batch.column(11).as_ref(), row),
                            ccre_promoter: bool_col_val(batch.column(12).as_ref(), row),
                            ccre_enhancer: bool_col_val(batch.column(13).as_ref(), row),
                        },
                        weights: AnnotationWeights(weights),
                    },
                });
            }
        }
        variants
    }

    fn setup_rare_table() -> (tempfile::TempDir, DfEngine, Vec<AnnotatedVariant>) {
        let dir = tempfile::tempdir().unwrap();
        let ann = write_test_annotation_parquet(dir.path());
        let geno = write_test_genotype_parquet(dir.path());

        let resources = Resources::detect();
        let engine = DfEngine::new(&resources).unwrap();
        engine.register_parquet_file("_test_geno", &geno).unwrap();
        engine.register_parquet_file("_test_ann", &ann).unwrap();
        engine.execute(&staar_rare_sql()).unwrap();

        let variants = read_extracted_variants(&engine);
        (dir, engine, variants)
    }

    #[test]
    fn extraction_preserves_all_variants() {
        let (_dir, _engine, variants) = setup_rare_table();
        assert_eq!(variants.len(), all_variants().len());
    }

    #[test]
    fn extraction_stopgain_fields() {
        let (_dir, _engine, variants) = setup_rare_table();
        let v = variants.iter().find(|v| v.position == 15528164).unwrap();
        assert_eq!(&*v.gene_name, "OR11H1");
        assert_eq!(v.annotation.region_type, RegionType::Exonic);
        assert_eq!(v.annotation.consequence, Consequence::Stopgain);
        assert!((v.annotation.cadd_phred - 23.7).abs() < 0.1);
        assert!(!v.annotation.regulatory.cage_promoter);
        assert!(!v.annotation.regulatory.cage_enhancer);
        assert!(!v.annotation.regulatory.ccre_promoter);
        assert!(!v.annotation.regulatory.ccre_enhancer);
    }

    #[test]
    fn extraction_splice_has_empty_consequence() {
        let (_dir, _engine, variants) = setup_rare_table();
        let v = variants
            .iter()
            .find(|v| v.annotation.region_type == RegionType::Splicing)
            .unwrap();
        assert_eq!(v.annotation.consequence, Consequence::Unknown);
        assert!((v.annotation.cadd_phred - 25.1).abs() < 0.1);
    }

    #[test]
    fn extraction_cage_flags() {
        let (_dir, _engine, variants) = setup_rare_table();
        let prom = variants.iter().find(|v| v.position == 15526000).unwrap();
        assert!(prom.annotation.regulatory.cage_promoter);
        assert!(!prom.annotation.regulatory.cage_enhancer);

        let enh = variants.iter().find(|v| v.position == 15525000).unwrap();
        assert!(!enh.annotation.regulatory.cage_promoter);
        assert!(enh.annotation.regulatory.cage_enhancer);
    }

    #[test]
    fn extraction_ccre_flags() {
        let (_dir, _engine, variants) = setup_rare_table();
        let pls = variants.iter().find(|v| v.position == 15524000).unwrap();
        assert!(pls.annotation.regulatory.ccre_promoter);
        assert!(!pls.annotation.regulatory.ccre_enhancer);

        let els = variants.iter().find(|v| v.position == 15523000).unwrap();
        assert!(!els.annotation.regulatory.ccre_promoter);
        assert!(els.annotation.regulatory.ccre_enhancer);
    }

    #[test]
    fn extraction_null_cage_is_false() {
        let (_dir, _engine, variants) = setup_rare_table();
        let intronic = variants
            .iter()
            .find(|v| v.annotation.region_type == RegionType::Intronic)
            .unwrap();
        assert!(!intronic.annotation.regulatory.cage_promoter);
        assert!(!intronic.annotation.regulatory.cage_enhancer);
    }

    #[test]
    fn extraction_annotation_weights_nonzero() {
        let (_dir, _engine, variants) = setup_rare_table();
        let v = variants
            .iter()
            .find(|v| v.annotation.consequence == Consequence::Stopgain)
            .unwrap();
        assert!(v.annotation.weights.0[0] > 0.99); // cadd weight for phred=23.7
        assert!(v.annotation.weights.0[1] > 0.0); // linsight
    }

    #[test]
    fn coding_masks_on_extracted() {
        let (_dir, _engine, variants) = setup_rare_table();
        let coding = masks::build_coding_masks(&variants, 1);

        let plof_genes: Vec<&str> = coding
            .iter()
            .filter(|(mt, _)| *mt == MaskType::PLof)
            .flat_map(|(_, groups)| groups.iter().map(|g| g.name.as_str()))
            .collect();
        assert!(plof_genes.contains(&"OR11H1"));

        let ptv = coding.iter().find(|(mt, _)| *mt == MaskType::Ptv).unwrap();
        for group in &ptv.1 {
            for &idx in &group.variant_indices {
                let csq = variants[idx].annotation.consequence;
                assert!(
                    csq == Consequence::Stopgain || csq == Consequence::FrameshiftDeletion,
                    "ptv should not include: {}",
                    csq.as_str()
                );
            }
        }
    }

    #[test]
    fn noncoding_masks_on_extracted() {
        let (_dir, _engine, variants) = setup_rare_table();
        let noncoding = masks::build_noncoding_masks(&variants, 1);

        let has_mask = |mt: MaskType| -> bool {
            noncoding
                .iter()
                .any(|(m, groups)| *m == mt && !groups.is_empty())
        };

        assert!(has_mask(MaskType::Upstream));
        assert!(has_mask(MaskType::Downstream));
        assert!(has_mask(MaskType::Utr));
        assert!(has_mask(MaskType::PromoterCage));
        assert!(has_mask(MaskType::PromoterDhs));
        assert!(has_mask(MaskType::EnhancerCage));
        assert!(has_mask(MaskType::EnhancerDhs));
        assert!(has_mask(MaskType::Ncrna));
    }

    #[test]
    fn intronic_in_no_mask() {
        let (_dir, _engine, variants) = setup_rare_table();
        let intronic_idx = variants
            .iter()
            .position(|v| v.annotation.region_type == RegionType::Intronic)
            .unwrap();

        let coding = masks::build_coding_masks(&variants, 1);
        let noncoding = masks::build_noncoding_masks(&variants, 1);

        for (_, groups) in coding.iter().chain(noncoding.iter()) {
            for group in groups {
                assert!(
                    !group.variant_indices.contains(&intronic_idx),
                    "intronic variant should not appear in any mask"
                );
            }
        }
    }

    #[test]
    fn sliding_windows_on_extracted() {
        let (_dir, _engine, variants) = setup_rare_table();
        let indices: Vec<usize> = (0..variants.len()).collect();
        let windows = masks::build_sliding_windows(
            &variants,
            &indices,
            crate::types::Chromosome::Autosome(22),
            2000,
            2000,
        );
        assert!(!windows.is_empty());
        for w in &windows {
            assert!(w.variant_indices.len() >= 2);
        }
    }
}
