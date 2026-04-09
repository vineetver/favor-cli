//! Build sparse_g.bin + variants.parquet + membership.parquet from _rare_all.

use std::collections::HashMap;
use std::fmt::Write as _;
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
use crate::error::CohortError;
use crate::output::Output;
use super::encoding::*;
use super::variants::col_by_name;

/// Stringify an Arrow cell used as a join key. Fast path for `StringArray`
/// (ingest + genotype parquet both produce these); falls back to the generic
/// Arrow display formatter so an upstream schema change still works, just
/// slower.
fn str_cell(col: &dyn Array, row: usize) -> String {
    if let Some(a) = col.as_any().downcast_ref::<StringArray>() {
        return a.value(row).to_string();
    }
    arrow::util::display::array_value_to_string(col, row).unwrap_or_default()
}

pub struct BuildStats {
    pub n_variants: usize,
    pub n_genes: usize,
    pub total_carriers: u64,
}

struct UniqueVariant {
    pos: i32,
    ref_a: Box<str>,
    alt_a: Box<str>,
    genes: Vec<Box<str>>,
}

pub fn build_chromosome(
    engine: &DfEngine,
    geno_parquet_path: &str,
    chrom: &str,
    chrom_dir: &Path,
    n_samples: usize,
    out: &dyn Output,
) -> Result<BuildStats, CohortError> {
    let wide = n_samples > 65535;

    let meta_batches = engine.collect(&column::carrier_metadata_sql(chrom))?;

    // Build unique_variants directly from meta_batches. The SQL orders rows
    // by (pos, ref, alt, gene) so same-variant rows arrive consecutively and
    // collapse into a single UniqueVariant whose `genes` list holds every
    // distinct gene for that locus.
    let mut unique_variants: Vec<UniqueVariant> = Vec::new();
    for batch in &meta_batches {
        let pos_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or_else(|| CohortError::Analysis("_rare_all missing position".into()))?;
        for i in 0..batch.num_rows() {
            let pos = pos_arr.value(i);
            let ref_a = str_cell(batch.column(1).as_ref(), i);
            let alt_a = str_cell(batch.column(2).as_ref(), i);
            let gene = str_cell(batch.column(3).as_ref(), i);

            if let Some(last) = unique_variants.last_mut() {
                if last.pos == pos && &*last.ref_a == ref_a.as_str() && &*last.alt_a == alt_a.as_str() {
                    if !gene.is_empty() && !last.genes.iter().any(|g| &**g == gene.as_str()) {
                        last.genes.push(gene.into_boxed_str());
                    }
                    continue;
                }
            }
            let genes = if gene.is_empty() {
                Vec::new()
            } else {
                vec![gene.into_boxed_str()]
            };
            unique_variants.push(UniqueVariant {
                pos,
                ref_a: ref_a.into_boxed_str(),
                alt_a: alt_a.into_boxed_str(),
                genes,
            });
        }
    }

    if unique_variants.is_empty() {
        return Ok(BuildStats {
            n_variants: 0,
            n_genes: 0,
            total_carriers: 0,
        });
    }

    let n_variants = unique_variants.len();

    let mut unique_pos: Vec<i32> = unique_variants.iter().map(|v| v.pos).collect();
    unique_pos.sort_unstable();
    unique_pos.dedup();
    let mut pos_str = String::with_capacity(unique_pos.len() * 8);
    for (i, p) in unique_pos.iter().enumerate() {
        if i > 0 {
            pos_str.push(',');
        }
        write!(&mut pos_str, "{p}").expect("write to String cannot fail");
    }

    engine.register_parquet_file("_geno_carrier", Path::new(geno_parquet_path))?;
    let dos_batches = engine.collect(&format!(
        "SELECT position, \"ref\", alt, dosages \
         FROM _geno_carrier \
         WHERE position IN ({pos_str}) \
         ORDER BY position, \"ref\", alt"
    ))?;

    // Sort-merge: both `unique_variants` (from meta_batches) and `dos_batches`
    // are sorted by (position, ref, alt). Walk them in lockstep so we never
    // build a HashMap keyed by the string columns. One cursor into
    // unique_variants advances past every dos row's smaller neighbors, and
    // matching rows write straight into the pre-sized `dosage_vecs`. Dos
    // rows with no matching uv are skipped; uvs with no matching dos keep
    // an empty carrier list.
    let mut dosage_vecs: Vec<Vec<(u32, u8)>> = vec![Vec::new(); n_variants];
    let mut uv_idx = 0usize;
    for batch in &dos_batches {
        let pos_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<Int32Array>()
            .ok_or_else(|| CohortError::Analysis("Genotype parquet missing position".into()))?;
        let ref_arr = batch
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| CohortError::Analysis("Genotype parquet ref not StringArray".into()))?;
        let alt_arr = batch
            .column(2)
            .as_any()
            .downcast_ref::<StringArray>()
            .ok_or_else(|| CohortError::Analysis("Genotype parquet alt not StringArray".into()))?;
        let dos_col = batch.column(3);
        let dos_fsl = dos_col.as_any().downcast_ref::<FixedSizeListArray>();
        let dos_list = dos_col.as_any().downcast_ref::<ListArray>();
        if dos_fsl.is_none() && dos_list.is_none() {
            return Err(CohortError::Analysis(
                "Genotype parquet missing dosages".into(),
            ));
        }

        for i in 0..batch.num_rows() {
            let dos_key = (pos_arr.value(i), ref_arr.value(i), alt_arr.value(i));

            while uv_idx < n_variants {
                let uv = &unique_variants[uv_idx];
                if (uv.pos, &*uv.ref_a, &*uv.alt_a) < dos_key {
                    uv_idx += 1;
                    continue;
                }
                break;
            }
            if uv_idx >= n_variants {
                break;
            }
            let uv = &unique_variants[uv_idx];
            if (uv.pos, &*uv.ref_a, &*uv.alt_a) != dos_key {
                continue;
            }

            let dosage_arr = if let Some(fsl) = dos_fsl {
                fsl.value(i)
            } else {
                dos_list.unwrap().value(i)
            };
            let dosages_f32 = dosage_arr.as_any().downcast_ref::<Float32Array>();
            let n_dos = dosage_arr.len();
            let carriers = &mut dosage_vecs[uv_idx];
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
        }
    }
    let _ = engine.execute("DROP TABLE IF EXISTS _geno_carrier");

    let sparse_g_path = chrom_dir.join("sparse_g.bin");
    let file = File::create(&sparse_g_path)
        .map_err(|e| CohortError::Resource(format!("Create sparse_g.bin: {e}")))?;
    let mut w = BufWriter::with_capacity(4 * 1024 * 1024, file);

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

    for (i, _uv) in unique_variants.iter().enumerate() {
        offsets.push(data_offset);

        // Sort-merge stashed this variant's carriers at `dosage_vecs[i]`.
        // Moving out avoids cloning the carrier list and lets the Vec's
        // backing allocation drop with it.
        let carriers = std::mem::take(&mut dosage_vecs[i]);

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

    let mut mem_variant_vcf = UInt32Builder::new();
    let mut mem_gene = StringBuilder::new();
    let mut n_genes_set: std::collections::HashSet<&str> = std::collections::HashSet::new();

    for (variant_vcf, uv) in unique_variants.iter().enumerate() {
        for gene in &uv.genes {
            mem_variant_vcf.append_value(variant_vcf as u32);
            mem_gene.append_value(&**gene);
            n_genes_set.insert(&**gene);
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
    .map_err(|e| CohortError::Resource(format!("membership batch: {e}")))?;

    write_parquet(chrom_dir, "membership.parquet", mem_schema, &mem_batch)?;

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
) -> Result<(), CohortError> {
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .build();
    let file = File::create(dir.join(filename))
        .map_err(|e| CohortError::Resource(format!("Create {filename}: {e}")))?;
    let mut writer = ArrowWriter::try_new(file, schema, Some(props))
        .map_err(|e| CohortError::Resource(format!("Parquet writer: {e}")))?;
    writer
        .write(batch)
        .map_err(|e| CohortError::Resource(format!("Write: {e}")))?;
    writer
        .close()
        .map_err(|e| CohortError::Resource(format!("Close: {e}")))?;
    Ok(())
}

/// Write variants.parquet with variant_vcf as primary key.
/// Metadata comes from the full _rare_all query, deduplicated to unique variants.
fn write_variants_parquet(
    chrom_dir: &Path,
    full_batches: &[RecordBatch],
    chrom: &str,
    unique_variants: &[UniqueVariant],
) -> Result<(), CohortError> {
    use arrow::array::BooleanArray;

    // (pos, ref, alt) → first batch row. The query is already sorted by
    // (pos, ref, alt) so the first hit is the canonical row.
    struct MetaRow {
        end_position: i32,
        maf: f64,
        region_type: Box<str>,
        consequence: Box<str>,
        cadd_phred: f64,
        revel: f64,
        cage_prom: bool,
        cage_enh: bool,
        ccre_prom: bool,
        ccre_enh: bool,
        weights: [f64; 11],
    }

    let mut meta_map: HashMap<(i32, Box<str>, Box<str>), MetaRow> = HashMap::new();
    for batch in full_batches {
        let pos_arr = col_by_name::<Int32Array>(batch, Col::Position.as_str())?;
        let end_arr = col_by_name::<Int32Array>(batch, Col::EndPosition.as_str())?;
        let ref_arr = col_by_name::<StringArray>(batch, Col::RefAllele.as_str())?;
        let alt_arr = col_by_name::<StringArray>(batch, Col::AltAllele.as_str())?;
        let maf_arr = col_by_name::<Float64Array>(batch, Col::Maf.as_str())?;
        let rt_arr = col_by_name::<StringArray>(batch, Col::RegionType.as_str())?;
        let csq_arr = col_by_name::<StringArray>(batch, Col::Consequence.as_str())?;
        let cadd_arr = col_by_name::<Float64Array>(batch, Col::CaddPhred.as_str())?;
        let revel_arr = col_by_name::<Float64Array>(batch, Col::Revel.as_str())?;
        let cps = col_by_name::<BooleanArray>(batch, Col::IsCagePromoter.as_str())?;
        let ces = col_by_name::<BooleanArray>(batch, Col::IsCageEnhancer.as_str())?;
        let crps = col_by_name::<BooleanArray>(batch, Col::IsCcrePromoter.as_str())?;
        let cres = col_by_name::<BooleanArray>(batch, Col::IsCcreEnhancer.as_str())?;
        let w_arrs: Vec<&Float64Array> = STAAR_WEIGHTS
            .iter()
            .map(|c| col_by_name::<Float64Array>(batch, c.as_str()))
            .collect::<Result<Vec<_>, _>>()?;

        for i in 0..batch.num_rows() {
            let key = (
                pos_arr.value(i),
                Box::<str>::from(ref_arr.value(i)),
                Box::<str>::from(alt_arr.value(i)),
            );
            if meta_map.contains_key(&key) {
                continue;
            } // keep first
            let mut w = [0.0f64; 11];
            for (ch, wa) in w_arrs.iter().enumerate() {
                w[ch] = wa.value(i);
            }
            meta_map.insert(
                key,
                MetaRow {
                    end_position: end_arr.value(i),
                    maf: maf_arr.value(i),
                    region_type: rt_arr.value(i).into(),
                    consequence: csq_arr.value(i).into(),
                    cadd_phred: cadd_arr.value(i),
                    revel: revel_arr.value(i),
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
    let mut end_b = Int32Builder::with_capacity(n);
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
        ref_b.append_value(&*uv.ref_a);
        alt_b.append_value(&*uv.alt_a);
        vid_b.append_value(crate::types::format_vid(
            chrom,
            uv.pos as u32,
            &uv.ref_a,
            &uv.alt_a,
        ));

        if let Some(m) = mr {
            end_b.append_value(m.end_position);
            maf_b.append_value(m.maf);
            rt_b.append_value(&*m.region_type);
            csq_b.append_value(&*m.consequence);
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
            end_b.append_value(uv.pos + uv.ref_a.len() as i32 - 1);
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
        Field::new(Col::EndPosition.as_str(), DataType::Int32, false),
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
        Arc::new(end_b.finish()),
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
        .map_err(|e| CohortError::Resource(format!("Arrow batch: {e}")))?;

    write_parquet(chrom_dir, "variants.parquet", schema, &batch)
}
