//! Build sparse_g.bin + variants.parquet + membership.parquet from _rare_all.
//!
//! Variants are sorted by (position, ref_allele, alt_allele) and assigned
//! monotonically increasing variant_vcf indices. Gene membership is stored
//! separately in membership.parquet — the sparse genotype matrix has no
//! gene concept.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Seek, SeekFrom, Write as IoWrite};
use std::path::Path;
use std::sync::Arc;

use arrow::array::{
    Array, ArrayRef, BooleanBuilder, FixedSizeListArray, Float32Array, Float64Array,
    Float64Builder, Int32Array, Int32Builder, ListArray, StringArray, StringBuilder, UInt32Builder,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use crate::column::{self, Col, STAAR_WEIGHTS};
use crate::engine::DfEngine;
use crate::error::FavorError;
use crate::output::Output;
use crate::staar::carrier::encoding::*;

pub struct BuildStats {
    pub n_variants: usize,
    pub n_genes: usize,
    pub total_carriers: u64,
}

/// A unique variant (deduplicated by position, ref, alt).
struct UniqueVariant {
    pos: i32,
    ref_a: String,
    alt_a: String,
    genes: Vec<String>,
}

pub fn build_chromosome(
    engine: &DfEngine,
    geno_parquet_path: &str,
    chrom: &str,
    chrom_dir: &Path,
    n_samples: usize,
    out: &dyn Output,
) -> Result<BuildStats, FavorError> {
    let wide = n_samples > 65535;

    // Query metadata sorted by (position, ref, alt) — canonical variant_vcf ordering.
    // Gene column is carried for membership.parquet.
    let meta_batches = engine.collect(&column::carrier_metadata_sql(chrom))?;

    let str_val = |col: &dyn Array, row: usize| -> String {
        if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
            return a.value(row).to_string();
        }
        arrow::util::display::array_value_to_string(col, row).unwrap_or_default()
    };

    // Collect raw rows sorted by (position, ref, alt)
    struct RawRow {
        pos: i32,
        ref_a: String,
        alt_a: String,
        gene: String,
    }
    let mut raw_rows: Vec<RawRow> = Vec::new();
    for batch in &meta_batches {
        let pos_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or_else(|| FavorError::Analysis("_rare_all missing position".into()))?;
        for i in 0..batch.num_rows() {
            raw_rows.push(RawRow {
                pos: pos_arr.value(i),
                ref_a: str_val(batch.column(1).as_ref(), i),
                alt_a: str_val(batch.column(2).as_ref(), i),
                gene: str_val(batch.column(3).as_ref(), i),
            });
        }
    }

    if raw_rows.is_empty() {
        return Ok(BuildStats {
            n_variants: 0,
            n_genes: 0,
            total_carriers: 0,
        });
    }

    // Deduplicate by (pos, ref, alt). A variant in multiple genes gets one
    // variant_vcf but multiple membership rows. Input is already sorted by
    // (pos, ref, alt) from the SQL ORDER BY.
    let mut unique_variants: Vec<UniqueVariant> = Vec::new();
    for row in &raw_rows {
        if let Some(last) = unique_variants.last_mut() {
            if last.pos == row.pos && last.ref_a == row.ref_a && last.alt_a == row.alt_a {
                if !row.gene.is_empty() && !last.genes.contains(&row.gene) {
                    last.genes.push(row.gene.clone());
                }
                continue;
            }
        }
        let genes = if row.gene.is_empty() {
            Vec::new()
        } else {
            vec![row.gene.clone()]
        };
        unique_variants.push(UniqueVariant {
            pos: row.pos,
            ref_a: row.ref_a.clone(),
            alt_a: row.alt_a.clone(),
            genes,
        });
    }

    let n_variants = unique_variants.len();

    // Unique positions for genotype IN clause
    let mut unique_pos: Vec<i32> = unique_variants.iter().map(|v| v.pos).collect();
    unique_pos.sort_unstable();
    unique_pos.dedup();
    let pos_str = unique_pos
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");

    // Load genotype dosages
    engine.register_parquet_file("_geno_carrier", Path::new(geno_parquet_path))?;
    let dos_batches = engine.collect(&format!(
        "SELECT position, \"ref\", alt, dosages \
         FROM _geno_carrier \
         WHERE position IN ({pos_str}) \
         ORDER BY position, \"ref\", alt"
    ))?;

    let mut dosage_map: HashMap<(i32, String, String), Vec<(u32, u8)>> = HashMap::new();
    for batch in &dos_batches {
        let pos_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or_else(|| FavorError::Analysis("Genotype parquet missing position".into()))?;
        let dos_col = batch.column(3);
        let dos_fsl = dos_col.as_any().downcast_ref::<FixedSizeListArray>();
        let dos_list = dos_col.as_any().downcast_ref::<ListArray>();
        if dos_fsl.is_none() && dos_list.is_none() {
            return Err(FavorError::Analysis(
                "Genotype parquet missing dosages".into(),
            ));
        }

        for i in 0..batch.num_rows() {
            let pos = pos_arr.value(i);
            let ref_a = str_val(batch.column(1).as_ref(), i);
            let alt_a = str_val(batch.column(2).as_ref(), i);
            let dosage_arr = if let Some(fsl) = dos_fsl {
                fsl.value(i)
            } else {
                dos_list.unwrap().value(i)
            };
            let dosages_f32 = dosage_arr.as_any().downcast_ref::<Float32Array>();

            let mut carriers = Vec::new();
            let n_dos = dosage_arr.len();
            for si in 0..n_samples.min(n_dos) {
                if dosage_arr.is_null(si) {
                    continue;
                }
                let d = if let Some(a) = dosages_f32 {
                    a.value(si) as f64
                } else if let Some(a) = dosage_arr.as_any().downcast_ref::<Float64Array>() {
                    a.value(si)
                } else {
                    continue;
                };
                if d > 0.0 && d.is_finite() {
                    carriers.push((si as u32, (d.round() as u8).min(2)));
                }
            }
            dosage_map.insert((pos, ref_a, alt_a), carriers);
        }
    }
    let _ = engine.execute("DROP TABLE IF EXISTS _geno_carrier");

    // ── Write sparse_g.bin ────────────────────────────────────────────────

    let sparse_g_path = chrom_dir.join("sparse_g.bin");
    let file = File::create(&sparse_g_path)
        .map_err(|e| FavorError::Resource(format!("Create sparse_g.bin: {e}")))?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, file);

    // Placeholder header — rewritten at end with correct totals
    let placeholder = SparseGHeader::new(n_samples as u32, n_variants as u32, 0, 0);
    placeholder.write_to(&mut w)?;

    let entry_size = if wide {
        CARRIER_ENTRY_WIDE
    } else {
        CARRIER_ENTRY_NARROW
    };
    let mut offsets: Vec<u64> = Vec::with_capacity(n_variants);
    let mut data_offset: u64 = 0;
    let mut total_carriers: u64 = 0;

    // Write carrier data in variant_vcf order (0, 1, 2, ...)
    for uv in &unique_variants {
        offsets.push(data_offset);

        let carriers = dosage_map
            .get(&(uv.pos, uv.ref_a.clone(), uv.alt_a.clone()))
            .cloned()
            .unwrap_or_default();

        let max_carriers = u16::MAX as usize;
        let n_carriers = carriers.len().min(max_carriers) as u16;
        w.write_all(&n_carriers.to_le_bytes())?;
        data_offset += CARRIER_COUNT_SIZE as u64;

        for &(sample_idx, dosage) in carriers.iter().take(max_carriers) {
            if wide {
                w.write_all(&sample_idx.to_le_bytes())?;
            } else {
                w.write_all(&(sample_idx as u16).to_le_bytes())?;
            }
            w.write_all(&[dosage])?;
            data_offset += entry_size as u64;
        }
        total_carriers += n_carriers as u64;
    }

    // Write offsets table
    let offsets_start = SPARSE_G_HEADER_SIZE as u64 + data_offset;
    for &off in &offsets {
        w.write_all(&off.to_le_bytes())?;
    }

    // Rewrite header with correct totals
    w.seek(SeekFrom::Start(0))?;
    let final_header = SparseGHeader::new(
        n_samples as u32,
        n_variants as u32,
        total_carriers,
        offsets_start,
    );
    final_header.write_to(&mut w)?;
    w.flush()?;

    // ── Write membership.parquet ──────────────────────────────────────────

    let mut mem_variant_vcf = UInt32Builder::new();
    let mut mem_gene = StringBuilder::new();
    let mut n_genes_set = std::collections::HashSet::new();

    for (variant_vcf, uv) in unique_variants.iter().enumerate() {
        for gene in &uv.genes {
            mem_variant_vcf.append_value(variant_vcf as u32);
            mem_gene.append_value(gene);
            n_genes_set.insert(gene.clone());
        }
    }

    let mem_schema = Arc::new(Schema::new(vec![
        Field::new(Col::VariantVcf.as_str(), DataType::UInt32, false),
        Field::new(Col::GeneName.as_str(), DataType::Utf8, false),
    ]));
    let mem_batch = RecordBatch::try_new(
        mem_schema.clone(),
        vec![
            Arc::new(mem_variant_vcf.finish()) as ArrayRef,
            Arc::new(mem_gene.finish()) as ArrayRef,
        ],
    )
    .map_err(|e| FavorError::Resource(format!("membership batch: {e}")))?;

    write_parquet(chrom_dir, "membership.parquet", mem_schema, &mem_batch)?;

    // ── Write variants.parquet ────────────────────────────────────────────
    // Query full metadata for these variants (already sorted by pos, ref, alt)

    let full_batches = engine.collect(&column::metadata_select_sql(chrom))?;
    write_variants_parquet(chrom_dir, &full_batches, chrom, &unique_variants)?;

    let n_genes = n_genes_set.len();
    out.status(&format!(
        "    chr{chrom}: {n_variants} variants, {n_genes} genes, {total_carriers} carriers"
    ));

    Ok(BuildStats {
        n_variants,
        n_genes,
        total_carriers,
    })
}

fn write_parquet(
    dir: &Path,
    filename: &str,
    schema: Arc<Schema>,
    batch: &RecordBatch,
) -> Result<(), FavorError> {
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let file = File::create(dir.join(filename))
        .map_err(|e| FavorError::Resource(format!("Create {filename}: {e}")))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| FavorError::Resource(format!("Parquet writer: {e}")))?;
    writer
        .write(batch)
        .map_err(|e| FavorError::Resource(format!("Write: {e}")))?;
    writer
        .close()
        .map_err(|e| FavorError::Resource(format!("Close: {e}")))?;
    Ok(())
}

/// Write variants.parquet with variant_vcf as primary key.
/// Metadata comes from the full _rare_all query, deduplicated to unique variants.
fn write_variants_parquet(
    chrom_dir: &Path,
    full_batches: &[RecordBatch],
    chrom: &str,
    unique_variants: &[UniqueVariant],
) -> Result<(), FavorError> {
    use arrow::array::{BooleanArray, Float32Array};

    // Build a lookup: (pos, ref, alt) -> first batch row for metadata extraction
    // The full query is already sorted by (pos, ref, alt).
    struct MetaRow {
        maf: f64,
        #[allow(dead_code)]
        gene: String,
        region_type: String,
        consequence: String,
        cadd_phred: f64,
        revel: f64,
        cage_prom: bool,
        cage_enh: bool,
        ccre_prom: bool,
        ccre_enh: bool,
        weights: [f64; 11],
    }

    let float_val = |col: &dyn Array, row: usize| -> f64 {
        if let Some(a) = col.as_any().downcast_ref::<Float64Array>() {
            a.value(row)
        } else if let Some(a) = col.as_any().downcast_ref::<Float32Array>() {
            a.value(row) as f64
        } else {
            0.0
        }
    };
    let str_val = |col: &dyn Array, row: usize| -> String {
        if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
            return a.value(row).to_string();
        }
        arrow::util::display::array_value_to_string(col, row).unwrap_or_default()
    };

    let mut meta_map: HashMap<(i32, String, String), MetaRow> = HashMap::new();
    for batch in full_batches {
        let pos_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let bool_col = |i: usize| {
            batch
                .column(i)
                .as_any()
                .downcast_ref::<BooleanArray>()
                .unwrap()
        };
        let cps = bool_col(9);
        let ces = bool_col(10);
        let crps = bool_col(11);
        let cres = bool_col(12);

        for i in 0..batch.num_rows() {
            let key = (
                pos_arr.value(i),
                str_val(batch.column(1).as_ref(), i),
                str_val(batch.column(2).as_ref(), i),
            );
            if meta_map.contains_key(&key) {
                continue;
            } // keep first
            let mut w = [0.0f64; 11];
            #[allow(clippy::needless_range_loop)]
            for ch in 0..11 {
                w[ch] = float_val(batch.column(13 + ch).as_ref(), i);
            }
            meta_map.insert(
                key,
                MetaRow {
                    maf: float_val(batch.column(3).as_ref(), i),
                    gene: str_val(batch.column(4).as_ref(), i),
                    region_type: str_val(batch.column(5).as_ref(), i),
                    consequence: str_val(batch.column(6).as_ref(), i),
                    cadd_phred: float_val(batch.column(7).as_ref(), i),
                    revel: float_val(batch.column(8).as_ref(), i),
                    cage_prom: cps.value(i),
                    cage_enh: ces.value(i),
                    ccre_prom: crps.value(i),
                    ccre_enh: cres.value(i),
                    weights: w,
                },
            );
        }
    }

    let n = unique_variants.len();
    let mut vvcf_b = UInt32Builder::with_capacity(n);
    let mut pos_b = Int32Builder::with_capacity(n);
    let mut ref_b = StringBuilder::with_capacity(n, n * 4);
    let mut alt_b = StringBuilder::with_capacity(n, n * 4);
    let mut vid_b = StringBuilder::with_capacity(n, n * 20);
    let mut maf_b = Float64Builder::with_capacity(n);
    let mut rt_b = StringBuilder::with_capacity(n, n * 16);
    let mut csq_b = StringBuilder::with_capacity(n, n * 16);
    let mut cadd_b = Float64Builder::with_capacity(n);
    let mut revel_b = Float64Builder::with_capacity(n);
    let mut cp_b = BooleanBuilder::with_capacity(n);
    let mut ce_b = BooleanBuilder::with_capacity(n);
    let mut crp_b = BooleanBuilder::with_capacity(n);
    let mut cre_b = BooleanBuilder::with_capacity(n);
    let mut w_builders: Vec<Float64Builder> =
        (0..11).map(|_| Float64Builder::with_capacity(n)).collect();

    for (variant_vcf, uv) in unique_variants.iter().enumerate() {
        let key = (uv.pos, uv.ref_a.clone(), uv.alt_a.clone());
        let mr = meta_map.get(&key);

        vvcf_b.append_value(variant_vcf as u32);
        pos_b.append_value(uv.pos);
        ref_b.append_value(&uv.ref_a);
        alt_b.append_value(&uv.alt_a);
        vid_b.append_value(crate::types::format_vid(
            chrom,
            uv.pos as u32,
            &uv.ref_a,
            &uv.alt_a,
        ));

        if let Some(m) = mr {
            maf_b.append_value(m.maf);
            rt_b.append_value(&m.region_type);
            csq_b.append_value(&m.consequence);
            cadd_b.append_value(m.cadd_phred);
            revel_b.append_value(m.revel);
            cp_b.append_value(m.cage_prom);
            ce_b.append_value(m.cage_enh);
            crp_b.append_value(m.ccre_prom);
            cre_b.append_value(m.ccre_enh);
            for (ch, wb) in w_builders.iter_mut().enumerate() {
                wb.append_value(m.weights[ch]);
            }
        } else {
            maf_b.append_value(0.0);
            rt_b.append_value("");
            csq_b.append_value("");
            cadd_b.append_value(0.0);
            revel_b.append_value(0.0);
            cp_b.append_value(false);
            ce_b.append_value(false);
            crp_b.append_value(false);
            cre_b.append_value(false);
            for wb in w_builders.iter_mut() {
                wb.append_value(0.0);
            }
        }
    }

    let mut fields = vec![
        Field::new(Col::VariantVcf.as_str(), DataType::UInt32, false),
        Field::new(Col::Position.as_str(), DataType::Int32, false),
        Field::new(Col::RefAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::AltAllele.as_str(), DataType::Utf8, false),
        Field::new(Col::Vid.as_str(), DataType::Utf8, false),
        Field::new(Col::Maf.as_str(), DataType::Float64, false),
        Field::new(Col::RegionType.as_str(), DataType::Utf8, false),
        Field::new(Col::Consequence.as_str(), DataType::Utf8, false),
        Field::new(Col::CaddPhred.as_str(), DataType::Float64, false),
        Field::new(Col::Revel.as_str(), DataType::Float64, false),
        Field::new(Col::IsCagePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCageEnhancer.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcrePromoter.as_str(), DataType::Boolean, false),
        Field::new(Col::IsCcreEnhancer.as_str(), DataType::Boolean, false),
    ];
    for col in &STAAR_WEIGHTS {
        fields.push(Field::new(col.as_str(), DataType::Float64, false));
    }
    let schema = Arc::new(Schema::new(fields));

    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(vvcf_b.finish()),
        Arc::new(pos_b.finish()),
        Arc::new(ref_b.finish()),
        Arc::new(alt_b.finish()),
        Arc::new(vid_b.finish()),
        Arc::new(maf_b.finish()),
        Arc::new(rt_b.finish()),
        Arc::new(csq_b.finish()),
        Arc::new(cadd_b.finish()),
        Arc::new(revel_b.finish()),
        Arc::new(cp_b.finish()),
        Arc::new(ce_b.finish()),
        Arc::new(crp_b.finish()),
        Arc::new(cre_b.finish()),
    ];
    for wb in &mut w_builders {
        columns.push(Arc::new(wb.finish()));
    }

    let batch = RecordBatch::try_new(schema.clone(), columns)
        .map_err(|e| FavorError::Resource(format!("Arrow batch: {e}")))?;

    write_parquet(chrom_dir, "variants.parquet", schema, &batch)
}
