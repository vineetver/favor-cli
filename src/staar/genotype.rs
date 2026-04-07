//! Multi-sample VCF → genotype parquet (dense packed).
//!
//! Dosages are packed into a single FLOAT[] list column per variant — not one
//! column per sample. Single compression envelope per chromosome
//! instead of N_SAMPLES separate column chunks. For 3000 samples, that's 3000x
//! less column metadata overhead and sequential I/O instead of random seeks.
//!
//! Schema: chromosome, position, ref, alt, maf, dosages FLOAT[N_SAMPLES]
//! Sample order matches the VCF header (stored in a sidecar).

use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{ArrayRef, FixedSizeListBuilder, Float32Builder, Int32Builder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use serde::{Deserialize, Serialize};

use crate::error::CohortError;
use crate::ingest::vcf::{normalize_chrom, open_vcf, parsimony_normalize};
use crate::output::Output;

#[derive(Serialize, Deserialize)]
pub struct GenotypeMeta {
    pub version: u32,
    pub n_samples: usize,
    pub chromosomes: Vec<String>,
    pub source_vcf: String,
}

pub struct GenotypeResult {
    pub sample_names: Vec<String>,
    pub output_dir: PathBuf,
}

pub fn extract_genotypes(
    vcf_path: &Path,
    output_dir: &Path,
    available_memory: u64,
    threads: usize,
    output: &dyn Output,
) -> Result<GenotypeResult, CohortError> {
    let reader = open_vcf(vcf_path, threads)?;
    let mut vcf_reader = noodles_vcf::io::Reader::new(reader);
    let header = vcf_reader
        .read_header()
        .map_err(|e| CohortError::Input(format!("VCF header: {e}")))?;

    let sample_names: Vec<String> = header
        .sample_names()
        .iter()
        .map(|s| s.to_string())
        .collect();
    let n_samples = sample_names.len();
    if n_samples == 0 {
        return Err(CohortError::Input(
            "VCF has no samples. STAAR requires a multi-sample VCF.".into(),
        ));
    }

    output.status(&format!(
        "Extracting genotypes: {} samples, {} decompression threads, {:.1}G memory",
        n_samples,
        threads,
        available_memory as f64 / (1024.0 * 1024.0 * 1024.0)
    ));

    let bytes_per_variant = (n_samples as u64) * 4 + 200;
    let batch_size = ((available_memory / 4) / bytes_per_variant).clamp(1000, 100_000) as usize;

    // Schema: 5 metadata columns + 1 packed dosage list
    let schema = Arc::new(packed_schema(n_samples));
    let geno_dir = output_dir.join("genotypes");
    std::fs::create_dir_all(&geno_dir).map_err(|e| {
        CohortError::Resource(format!("Cannot create '{}': {e}", geno_dir.display()))
    })?;

    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_row_count(Some(batch_size))
        .build();

    let mut state = ExtractState {
        current_chrom: None,
        writer: None,
        batch: PackedBatchBuilder::new(n_samples, batch_size),
        total_variants: 0,
        chromosomes: Vec::new(),
    };

    let pb = output.progress(0, "extracting genotypes");
    for result in vcf_reader.records() {
        let record = result.map_err(|e| {
            CohortError::Analysis(format!("VCF parse error in '{}': {e}", vcf_path.display()))
        })?;
        process_record(
            &record, n_samples, &mut state, &schema, &props, &geno_dir, output,
        )?;
        pb.inc(1);
    }
    pb.finish(&format!("{} variants extracted", state.total_variants));

    flush(&mut state, &schema)?;
    if let Some(w) = state.writer.take() {
        w.close()
            .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
    }

    // Write sample names sidecar (order matters for unpacking dosages)
    let sidecar = geno_dir.join("samples.txt");
    std::fs::write(&sidecar, sample_names.join("\n"))
        .map_err(|e| CohortError::Resource(format!("Cannot write '{}': {e}", sidecar.display())))?;

    // Write metadata — this is the last step so partial extractions are detectable
    let meta = GenotypeMeta {
        version: 1,
        n_samples,
        chromosomes: state.chromosomes.clone(),
        source_vcf: vcf_path.display().to_string(),
    };
    let meta_path = geno_dir.join("genotypes.json");
    std::fs::write(
        &meta_path,
        serde_json::to_string_pretty(&meta)
            .map_err(|e| CohortError::Resource(format!("JSON serialize: {e}")))?,
    )
    .map_err(|e| CohortError::Resource(format!("Cannot write '{}': {e}", meta_path.display())))?;

    output.success(&format!(
        "Extracted {} variants × {} samples → packed dosage lists",
        state.total_variants, n_samples,
    ));

    Ok(GenotypeResult {
        sample_names,
        output_dir: geno_dir,
    })
}

struct ExtractState {
    current_chrom: Option<String>,
    writer: Option<ArrowWriter<File>>,
    batch: PackedBatchBuilder,
    total_variants: u64,
    chromosomes: Vec<String>,
}

fn flush(state: &mut ExtractState, schema: &Arc<Schema>) -> Result<(), CohortError> {
    if state.batch.count > 0 {
        if let Some(w) = state.writer.as_mut() {
            let rb = state.batch.finish(schema)?;
            w.write(&rb)
                .map_err(|e| CohortError::Resource(format!("Parquet write: {e}")))?;
        }
    }
    Ok(())
}

fn switch_chrom(
    chrom: &str,
    state: &mut ExtractState,
    schema: &Arc<Schema>,
    props: &WriterProperties,
    geno_dir: &Path,
    output: &dyn Output,
) -> Result<(), CohortError> {
    flush(state, schema)?;
    if let Some(w) = state.writer.take() {
        w.close()
            .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
    }
    let chr_dir = geno_dir.join(format!("chromosome={chrom}"));
    std::fs::create_dir_all(&chr_dir)
        .map_err(|e| CohortError::Resource(format!("Cannot create '{}': {e}", chr_dir.display())))?;
    let parquet_path = chr_dir.join("data.parquet");
    let f = File::create(&parquet_path).map_err(|e| {
        CohortError::Resource(format!("Cannot create '{}': {e}", parquet_path.display()))
    })?;
    state.writer = Some(
        ArrowWriter::try_new(f, schema.clone(), Some(props.clone()))
            .map_err(|e| CohortError::Resource(format!("Writer init: {e}")))?,
    );
    state.current_chrom = Some(chrom.to_string());
    state.chromosomes.push(chrom.to_string());
    output.status(&format!("  chr{chrom}..."));
    Ok(())
}

fn process_record(
    record: &noodles_vcf::Record,
    n_samples: usize,
    state: &mut ExtractState,
    schema: &Arc<Schema>,
    props: &WriterProperties,
    geno_dir: &Path,
    output: &dyn Output,
) -> Result<(), CohortError> {
    let raw_chrom = record.reference_sequence_name();
    let chrom = match normalize_chrom(raw_chrom) {
        Some(c) => c,
        None => return Ok(()),
    };

    if state.current_chrom.as_deref() != Some(chrom) {
        switch_chrom(chrom, state, schema, props, geno_dir, output)?;
    }

    let pos = match record.variant_start() {
        Some(Ok(p)) => p.get() as i32,
        _ => return Ok(()),
    };

    let ref_allele = record.reference_bases().to_uppercase();
    let alt_bases = record.alternate_bases();
    let alt_str: &str = alt_bases.as_ref();
    let alts: Vec<&str> = alt_str.split(',').collect();

    let samples_raw = record.samples();
    let samples_str: &str = samples_raw.as_ref();
    let gt_index = 0; // GT is first FORMAT field by VCF spec

    let sample_fields: Vec<&str> = if samples_str.is_empty() || samples_str == "." {
        Vec::new()
    } else {
        samples_str.split('\t').collect()
    };

    for (alt_idx, alt) in alts.iter().enumerate() {
        let alt_upper = alt.trim().to_uppercase();
        if alt_upper == "*" || alt_upper == "." || alt_upper.is_empty() {
            continue;
        }

        let (norm_ref, norm_alt, norm_pos) = parsimony_normalize(&ref_allele, &alt_upper, pos);

        // Parse genotypes — stack buffer, no heap allocation for typical diploid GTs
        let mut dosages = vec![f32::NAN; n_samples];
        let mut ac: f64 = 0.0;
        let mut an: f64 = 0.0;

        for (i, sf) in sample_fields.iter().enumerate().take(n_samples) {
            let gt = extract_gt_field(sf, gt_index);
            let dose = gt_to_dosage(gt.as_bytes(), (alt_idx + 1) as u8);
            dosages[i] = dose;
            if dose.is_finite() {
                ac += dose as f64;
                an += 2.0;
            }
        }

        let af = if an > 0.0 { ac / an } else { 0.0 };
        let maf = af.min(1.0 - af) as f32;

        state
            .batch
            .push(chrom, norm_pos, &norm_ref, &norm_alt, maf, &dosages);
        state.total_variants += 1;

        if state.batch.is_full() {
            flush(state, schema)?;
        }
    }

    Ok(())
}

fn packed_schema(n_samples: usize) -> Schema {
    Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("position", DataType::Int32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
        Field::new("maf", DataType::Float32, true),
        Field::new(
            "dosages",
            DataType::FixedSizeList(
                Arc::new(Field::new("item", DataType::Float32, true)),
                n_samples as i32,
            ),
            false,
        ),
    ])
}

struct PackedBatchBuilder {
    chromosome: StringBuilder,
    position: Int32Builder,
    ref_allele: StringBuilder,
    alt_allele: StringBuilder,
    maf: Float32Builder,
    dosages: FixedSizeListBuilder<Float32Builder>,
    count: usize,
    capacity: usize,
}

impl PackedBatchBuilder {
    fn new(n_samples: usize, capacity: usize) -> Self {
        Self {
            chromosome: StringBuilder::with_capacity(capacity, capacity * 3),
            position: Int32Builder::with_capacity(capacity),
            ref_allele: StringBuilder::with_capacity(capacity, capacity * 4),
            alt_allele: StringBuilder::with_capacity(capacity, capacity * 4),
            maf: Float32Builder::with_capacity(capacity),
            dosages: FixedSizeListBuilder::with_capacity(
                Float32Builder::with_capacity(capacity * n_samples),
                n_samples as i32,
                capacity,
            ),
            count: 0,
            capacity,
        }
    }

    fn push(&mut self, chrom: &str, pos: i32, ref_a: &str, alt_a: &str, maf: f32, doses: &[f32]) {
        self.chromosome.append_value(chrom);
        self.position.append_value(pos);
        self.ref_allele.append_value(ref_a);
        self.alt_allele.append_value(alt_a);
        self.maf.append_value(maf);

        let values = self.dosages.values();
        for &d in doses {
            if d.is_nan() {
                values.append_null();
            } else {
                values.append_value(d);
            }
        }
        self.dosages.append(true);

        self.count += 1;
    }

    fn is_full(&self) -> bool {
        self.count >= self.capacity
    }

    fn finish(&mut self, schema: &Arc<Schema>) -> Result<RecordBatch, CohortError> {
        let columns: Vec<ArrayRef> = vec![
            Arc::new(self.chromosome.finish()),
            Arc::new(self.position.finish()),
            Arc::new(self.ref_allele.finish()),
            Arc::new(self.alt_allele.finish()),
            Arc::new(self.maf.finish()),
            Arc::new(self.dosages.finish()),
        ];
        self.count = 0;
        RecordBatch::try_new(schema.clone(), columns)
            .map_err(|e| CohortError::Resource(format!("Arrow batch: {e}")))
    }
}

/// Extract GT field from a colon-delimited sample field. Zero allocation.
#[inline]
pub fn extract_gt_field(sample_field: &str, gt_index: usize) -> &str {
    let mut start = 0;
    let mut field_idx = 0;
    let bytes = sample_field.as_bytes();
    for i in 0..bytes.len() {
        if bytes[i] == b':' {
            if field_idx == gt_index {
                return &sample_field[start..i];
            }
            field_idx += 1;
            start = i + 1;
        }
    }
    if field_idx == gt_index {
        &sample_field[start..]
    } else {
        "."
    }
}

/// Parse GT bytes to dosage. Branchless for common cases, no allocation.
/// "0/0" → 0.0, "0/1" → 1.0, "1/1" → 2.0, "./." → NaN
#[inline]
pub fn gt_to_dosage(gt: &[u8], alt_index: u8) -> f32 {
    // Fast path: most common genotypes are 3 bytes (0/0, 0/1, 1/1, ./.)
    if gt.len() == 3 {
        let a0 = gt[0];
        let a1 = gt[2];
        if a0 == b'.' || a1 == b'.' {
            return f32::NAN;
        }
        let d0 = if a0.wrapping_sub(b'0') == alt_index {
            1.0f32
        } else {
            0.0
        };
        let d1 = if a1.wrapping_sub(b'0') == alt_index {
            1.0f32
        } else {
            0.0
        };
        return d0 + d1;
    }
    // Slow path: multi-digit alleles or missing
    gt_to_dosage_slow(gt, alt_index)
}

fn gt_to_dosage_slow(gt: &[u8], alt_index: u8) -> f32 {
    let mut dose = 0.0f32;
    for allele in std::str::from_utf8(gt).unwrap_or(".").split(['/', '|']) {
        if allele == "." {
            return f32::NAN;
        }
        if let Ok(idx) = allele.parse::<u8>() {
            if idx == alt_index {
                dose += 1.0;
            }
        }
    }
    dose
}

pub fn read_sample_names(vcf_path: &Path) -> Result<Vec<String>, CohortError> {
    // Header-only read — single decompression thread is sufficient.
    let reader = open_vcf(vcf_path, 1)?;
    let mut vcf_reader = noodles_vcf::io::Reader::new(reader);
    let header = vcf_reader
        .read_header()
        .map_err(|e| CohortError::Input(format!("VCF header: {e}")))?;
    Ok(header
        .sample_names()
        .iter()
        .map(|s| s.to_string())
        .collect())
}

use crate::engine::DfEngine;

/// Load genotypes for specific positions into flat buffer: `buf[variant * n_samples + sample]`.
/// Used by known-loci conditional analysis path.
pub fn load(
    engine: &DfEngine,
    geno_path: &str,
    positions: &[u32],
    n_samples: usize,
) -> Result<Vec<f64>, CohortError> {
    use arrow::array::{Array, Float64Array, ListArray};
    engine.register_parquet_file("_geno_load", std::path::Path::new(geno_path))?;
    let pos_str = positions
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let batches = engine.collect(&format!(
        "SELECT dosages FROM _geno_load WHERE position IN ({pos_str}) ORDER BY position"
    ))?;

    let n_variants = positions.len();
    let mut flat = vec![0.0; n_variants * n_samples];
    let mut vi = 0;
    for batch in &batches {
        let dos_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<ListArray>()
            .ok_or_else(|| {
                CohortError::Analysis("Genotype parquet missing list-typed dosages column".into())
            })?;
        for i in 0..batch.num_rows() {
            if vi >= n_variants {
                break;
            }
            let dosage_list = dos_arr.value(i);
            let dosages = dosage_list
                .as_any()
                .downcast_ref::<Float64Array>()
                .ok_or_else(|| {
                    CohortError::Analysis("Dosage list elements are not Float64".into())
                })?;
            let base = vi * n_samples;
            for si in 0..n_samples.min(dosages.len()) {
                flat[base + si] = if dosages.is_null(si) {
                    0.0
                } else {
                    dosages.value(si)
                };
            }
            vi += 1;
        }
    }
    let _ = engine.execute("DROP TABLE IF EXISTS _geno_load");
    Ok(flat)
}
