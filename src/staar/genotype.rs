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
    pub source_vcfs: Vec<String>,
}

pub struct GenotypeResult {
    pub sample_names: Vec<String>,
    pub output_dir: PathBuf,
}

pub fn extract_genotypes(
    vcf_paths: &[PathBuf],
    output_dir: &Path,
    available_memory: u64,
    threads: usize,
    output: &dyn Output,
) -> Result<GenotypeResult, CohortError> {
    if vcf_paths.is_empty() {
        return Err(CohortError::Input("No VCF files provided.".into()));
    }

    let reader = open_vcf(&vcf_paths[0], threads)?;
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
    drop(vcf_reader);

    output.status(&format!(
        "Extracting genotypes: {} samples, {} file(s), {} threads, {:.1}G memory",
        n_samples,
        vcf_paths.len(),
        threads,
        available_memory as f64 / (1024.0 * 1024.0 * 1024.0)
    ));

    let geno_dir = output_dir.join("genotypes");
    let mut gw = GenotypeWriter::new(n_samples, &geno_dir, available_memory)?;

    let mut finished_chroms: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut ref_buf = String::with_capacity(256);
    let mut alt_buf = String::with_capacity(256);
    let pb = output.progress(0, "extracting genotypes");

    for (file_idx, vcf_path) in vcf_paths.iter().enumerate() {
        let reader = open_vcf(vcf_path, threads)?;
        let mut vcf_reader = noodles_vcf::io::Reader::new(reader);
        let file_header = vcf_reader
            .read_header()
            .map_err(|e| CohortError::Input(format!("VCF header in '{}': {e}", vcf_path.display())))?;

        if file_idx > 0 {
            let file_samples: Vec<String> = file_header
                .sample_names()
                .iter()
                .map(|s| s.to_string())
                .collect();
            if file_samples != sample_names {
                return Err(CohortError::Input(format!(
                    "Sample mismatch: '{}' has {} samples but '{}' has {}. \
                     All VCF files must have identical samples in the same order.",
                    vcf_paths[0].display(),
                    sample_names.len(),
                    vcf_path.display(),
                    file_samples.len()
                )));
            }
        }

        for result in vcf_reader.records() {
            let record = result.map_err(|e| {
                CohortError::Analysis(format!("VCF parse error in '{}': {e}", vcf_path.display()))
            })?;

            let raw_chrom = record.reference_sequence_name();
            if let Some(chrom) = normalize_chrom(raw_chrom) {
                if gw.current_chrom.as_deref() != Some(chrom) && finished_chroms.contains(chrom) {
                    return Err(CohortError::Input(format!(
                        "Chromosome {chrom} appears in '{}' but its writer was already \
                         closed after processing an earlier file. Sort VCF files by \
                         chromosome or use one file per chromosome.",
                        vcf_path.display()
                    )));
                }
            }

            process_record_geno(&record, &mut gw, &mut ref_buf, &mut alt_buf, output)?;
            pb.inc(1);
        }

        if let Some(ref current) = gw.current_chrom {
            for c in &gw.chromosomes {
                if c != current {
                    finished_chroms.insert(c.clone());
                }
            }
        }
    }
    pb.finish(&format!("{} variants extracted", gw.variant_count()));

    let total_variants = gw.variant_count();
    let source_vcfs = vcf_paths.iter().map(|p| p.display().to_string()).collect();
    let result = gw.finish(&sample_names, source_vcfs)?;

    output.success(&format!(
        "Extracted {total_variants} variants × {n_samples} samples from {} file(s)",
        vcf_paths.len(),
    ));

    Ok(result)
}

/// Reusable genotype writer. Accepts already-normalized variants and raw
/// VCF sample text, writes per-chromosome dosage parquets. Used both by
/// the standalone `extract_genotypes` path and by the single-pass ingest
/// path in `ingest/vcf.rs`.
pub struct GenotypeWriter {
    current_chrom: Option<String>,
    writer: Option<ArrowWriter<File>>,
    batch: PackedBatchBuilder,
    total_variants: u64,
    chromosomes: Vec<String>,
    dosages: Vec<f32>,
    n_samples: usize,
    schema: Arc<Schema>,
    props: WriterProperties,
    geno_dir: PathBuf,
}

impl GenotypeWriter {
    pub fn new(
        n_samples: usize,
        output_dir: &Path,
        available_memory: u64,
    ) -> Result<Self, CohortError> {
        let bytes_per_variant = (n_samples as u64) * 4 + 200;
        let batch_size =
            ((available_memory / 4) / bytes_per_variant).clamp(1000, 100_000) as usize;

        let geno_dir = output_dir.to_path_buf();
        std::fs::create_dir_all(&geno_dir).map_err(|e| {
            CohortError::Resource(format!("Cannot create '{}': {e}", geno_dir.display()))
        })?;

        let schema = Arc::new(packed_schema(n_samples));
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .set_max_row_group_row_count(Some(batch_size))
            .build();

        Ok(Self {
            current_chrom: None,
            writer: None,
            batch: PackedBatchBuilder::new(n_samples, batch_size),
            total_variants: 0,
            chromosomes: Vec::new(),
            dosages: vec![f32::NAN; n_samples],
            n_samples,
            schema,
            props,
            geno_dir,
        })
    }

    /// Push one biallelic variant with dosages extracted from raw VCF sample text.
    /// `chrom` and `pos/ref/alt` must already be normalized.
    /// `alt_idx` is the 1-based index of this alt allele in the original record.
    #[allow(clippy::too_many_arguments)]
    pub fn push(
        &mut self,
        chrom: &str,
        pos: i32,
        ref_allele: &str,
        alt_allele: &str,
        samples_str: &str,
        alt_idx: u8,
        output: &dyn Output,
    ) -> Result<(), CohortError> {
        if self.current_chrom.as_deref() != Some(chrom) {
            self.switch_chrom(chrom, output)?;
        }

        self.dosages.fill(f32::NAN);
        let mut ac: f64 = 0.0;
        let mut an: f64 = 0.0;

        if !samples_str.is_empty() && samples_str != "." {
            for (i, sf) in samples_str.split('\t').enumerate().take(self.n_samples) {
                let gt = extract_gt_field(sf, 0);
                let dose = gt_to_dosage(gt.as_bytes(), alt_idx);
                self.dosages[i] = dose;
                if dose.is_finite() {
                    ac += dose as f64;
                    an += 2.0;
                }
            }
        }

        let af = if an > 0.0 { ac / an } else { 0.0 };
        let maf = af.min(1.0 - af) as f32;

        let dosages: &[f32] = &self.dosages;
        self.batch
            .push(chrom, pos, ref_allele, alt_allele, maf, dosages);
        self.total_variants += 1;

        if self.batch.is_full() {
            self.flush()?;
        }
        Ok(())
    }

    pub fn variant_count(&self) -> u64 {
        self.total_variants
    }

    /// Finalize: flush remaining batches, close writers, write samples sidecar
    /// and metadata. Returns the output directory path and chromosome list.
    pub fn finish(
        mut self,
        sample_names: &[String],
        source_vcfs: Vec<String>,
    ) -> Result<GenotypeResult, CohortError> {
        self.flush()?;
        if let Some(w) = self.writer.take() {
            w.close()
                .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
        }

        let sidecar = self.geno_dir.join("samples.txt");
        std::fs::write(&sidecar, sample_names.join("\n"))
            .map_err(|e| {
                CohortError::Resource(format!("Cannot write '{}': {e}", sidecar.display()))
            })?;

        let meta = GenotypeMeta {
            version: 1,
            n_samples: self.n_samples,
            chromosomes: self.chromosomes.clone(),
            source_vcfs,
        };
        let meta_path = self.geno_dir.join("genotypes.json");
        std::fs::write(
            &meta_path,
            serde_json::to_string_pretty(&meta)
                .map_err(|e| CohortError::Resource(format!("JSON serialize: {e}")))?,
        )
        .map_err(|e| {
            CohortError::Resource(format!("Cannot write '{}': {e}", meta_path.display()))
        })?;

        Ok(GenotypeResult {
            sample_names: sample_names.to_vec(),
            output_dir: self.geno_dir,
        })
    }

    fn flush(&mut self) -> Result<(), CohortError> {
        if self.batch.count > 0 {
            if let Some(w) = self.writer.as_mut() {
                let rb = self.batch.finish(&self.schema)?;
                w.write(&rb)
                    .map_err(|e| CohortError::Resource(format!("Parquet write: {e}")))?;
            }
        }
        Ok(())
    }

    fn switch_chrom(&mut self, chrom: &str, output: &dyn Output) -> Result<(), CohortError> {
        self.flush()?;
        if let Some(w) = self.writer.take() {
            w.close()
                .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
        }
        let chr_dir = self.geno_dir.join(format!("chromosome={chrom}"));
        std::fs::create_dir_all(&chr_dir).map_err(|e| {
            CohortError::Resource(format!("Cannot create '{}': {e}", chr_dir.display()))
        })?;
        let parquet_path = chr_dir.join("data.parquet");
        let f = File::create(&parquet_path).map_err(|e| {
            CohortError::Resource(format!("Cannot create '{}': {e}", parquet_path.display()))
        })?;
        self.writer = Some(
            ArrowWriter::try_new(f, self.schema.clone(), Some(self.props.clone()))
                .map_err(|e| CohortError::Resource(format!("Writer init: {e}")))?,
        );
        self.current_chrom = Some(chrom.to_string());
        self.chromosomes.push(chrom.to_string());
        output.status(&format!("  chr{chrom} (genotypes)..."));
        Ok(())
    }
}

fn process_record_geno(
    record: &noodles_vcf::Record,
    gw: &mut GenotypeWriter,
    ref_buf: &mut String,
    alt_buf: &mut String,
    output: &dyn Output,
) -> Result<(), CohortError> {
    let chrom = match normalize_chrom(record.reference_sequence_name()) {
        Some(c) => c,
        None => return Ok(()),
    };
    let pos = match record.variant_start() {
        Some(Ok(p)) => p.get() as i32,
        _ => return Ok(()),
    };

    ascii_uppercase_into(record.reference_bases(), ref_buf);
    let alt_bases = record.alternate_bases();
    let alt_str: &str = alt_bases.as_ref();
    let samples_raw = record.samples();
    let samples_str: &str = samples_raw.as_ref();

    for (alt_idx, alt) in alt_str.split(',').enumerate() {
        ascii_uppercase_into(alt.trim(), alt_buf);
        if alt_buf == "*" || alt_buf == "." || alt_buf.is_empty() || alt_buf.starts_with('<') {
            continue;
        }

        let (nr, na, np) = parsimony_normalize(ref_buf, alt_buf, pos);
        gw.push(chrom, np, nr, na, samples_str, (alt_idx + 1) as u8, output)?;
    }
    Ok(())
}

fn ascii_uppercase_into(src: &str, buf: &mut String) {
    buf.clear();
    for &b in src.as_bytes() {
        buf.push(b.to_ascii_uppercase() as char);
    }
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
