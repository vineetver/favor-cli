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
use crate::ingest::vcf::open_vcf;
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

    let sample_names = read_sample_names(&vcf_paths[0])?;
    let n_samples = sample_names.len();
    if n_samples == 0 {
        return Err(CohortError::Input(
            "VCF has no samples. STAAR requires a multi-sample VCF.".into(),
        ));
    }

    for vcf_path in &vcf_paths[1..] {
        let file_samples = read_sample_names(vcf_path)?;
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

    output.status(&format!(
        "Extracting genotypes: {} samples, {} file(s), {} threads, {:.1}G memory",
        n_samples,
        vcf_paths.len(),
        threads,
        available_memory as f64 / (1024.0 * 1024.0 * 1024.0)
    ));

    let geno_dir = output_dir.join("genotypes");
    let mut gw = GenotypeWriter::new(n_samples, &geno_dir, available_memory)?;

    crate::ingest::vcf::stream_genotypes(
        vcf_paths,
        &mut gw,
        available_memory,
        threads,
        geno_dir.clone(),
        output,
    )?;

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
    part_id: Option<usize>,
}

impl GenotypeWriter {
    pub fn new(
        n_samples: usize,
        output_dir: &Path,
        available_memory: u64,
    ) -> Result<Self, CohortError> {
        Self::with_part_id(n_samples, output_dir, available_memory, None)
    }

    pub fn with_part_id(
        n_samples: usize,
        output_dir: &Path,
        available_memory: u64,
        part_id: Option<usize>,
    ) -> Result<Self, CohortError> {
        let batch_size =
            raw_batch_size(n_samples, available_memory).clamp(1000, 100_000) as usize;

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
            part_id,
        })
    }

    /// Chromosome currently being written, or `None` before the first push.
    pub fn current_chrom(&self) -> Option<&str> {
        self.current_chrom.as_deref()
    }

    /// Every chromosome touched so far in insertion order, including the one
    /// currently open.
    pub fn chromosomes(&self) -> &[String] {
        &self.chromosomes
    }

    /// Push one biallelic variant. `chrom` and `pos/ref/alt` must already be
    /// normalized. `alt_idx` is the 1-based index of this alt allele in the
    /// original VCF record.
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

        // noodles exposes `Samples` as everything after INFO, which in VCF
        // includes the FORMAT column (e.g. "GT" or "GT:AD:DP"). Skip past
        // the first tab so the memchr loop sees only per-sample fields.
        let body = if samples_str.is_empty() || samples_str == "." {
            ""
        } else {
            match memchr::memchr(b'\t', samples_str.as_bytes()) {
                Some(tab) => &samples_str[tab + 1..],
                None => "",
            }
        };

        if !body.is_empty() {
            let bytes = body.as_bytes();
            let n = self.n_samples;
            let mut sample_idx: usize = 0;
            let mut start: usize = 0;

            // memchr uses AVX2/SSE4.2 to find tabs at 32 bytes/cycle. For
            // 200K-sample lines (~7MB), this replaces the dominant cost of
            // byte-by-byte split('\t').
            for tab_pos in memchr::memchr_iter(b'\t', bytes) {
                if sample_idx >= n { break; }
                let gt_len = gt_prefix_len(&bytes[start..tab_pos]);
                let dose = gt_to_dosage(&bytes[start..start + gt_len], alt_idx);
                unsafe { *self.dosages.get_unchecked_mut(sample_idx) = dose; }
                let finite = dose.is_finite();
                ac += (finite as u8 as f64) * (dose as f64);
                an += (finite as u8 as f64) * 2.0;
                sample_idx += 1;
                start = tab_pos + 1;
            }
            if sample_idx < n && start < bytes.len() {
                let gt_len = gt_prefix_len(&bytes[start..]);
                let dose = gt_to_dosage(&bytes[start..start + gt_len], alt_idx);
                self.dosages[sample_idx] = dose;
                if dose.is_finite() { ac += dose as f64; an += 2.0; }
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

        // `Drop` is implemented on `Self`, so fields can't be partial-moved;
        // swap `geno_dir` out with a default instead.
        Ok(GenotypeResult {
            sample_names: sample_names.to_vec(),
            output_dir: std::mem::take(&mut self.geno_dir),
        })
    }

    /// Flush remaining batches and close the writer. Does NOT write sidecar
    /// files (samples.txt, genotypes.json). Used by parallel workers where the
    /// orchestrator handles sidecars after all workers join.
    pub fn flush_all(&mut self) -> Result<(), CohortError> {
        self.flush()?;
        if let Some(w) = self.writer.take() {
            w.close()
                .map_err(|e| CohortError::Resource(format!("Parquet close: {e}")))?;
        }
        Ok(())
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
        let parquet_path = match self.part_id {
            Some(id) => chr_dir.join(format!("part_{id}.parquet")),
            None => chr_dir.join("data.parquet"),
        };
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

// Close parquet footer on unwind so a panicking rayon worker can't leave a
// truncated, unreadable parquet behind. Normal paths (`flush_all`, `finish`,
// `switch_chrom`) already `take` the writer, so this only fires during panic.
impl Drop for GenotypeWriter {
    fn drop(&mut self) {
        if let Some(w) = self.writer.take() {
            let _ = w.close();
        }
    }
}

/// Row-group capacity the GenotypeWriter would allocate under `available_memory`,
/// before the 1k..100k clamp. Shared with the ingest preflight so it can reject
/// configurations whose raw capacity would force flush-per-variant throughput.
///
/// The `/ 4` factor matches the writer's internal reservation: the
/// FixedSizeListBuilder holds `batch_size * n_samples * 4` bytes, Arrow's
/// `finish()` briefly doubles that during copy-out.
pub fn raw_batch_size(n_samples: usize, available_memory: u64) -> u64 {
    let bytes_per_variant = (n_samples as u64) * 4 + 200;
    (available_memory / 4) / bytes_per_variant
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

/// Length of the GT prefix in a sample field byte slice.
/// GT is FORMAT field 0, so it ends at the first ':' (or the whole slice).
/// For the common 3-byte case (0/0, 0/1, ./.) with a colon at byte 3,
/// the branch hits immediately and memchr is never called.
#[inline]
fn gt_prefix_len(seg: &[u8]) -> usize {
    if seg.len() >= 4 && seg[3] == b':' { return 3; }
    memchr::memchr(b':', seg).unwrap_or(seg.len())
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

#[cfg(test)]
mod tests {
    use super::*;
    use std::panic::{self, AssertUnwindSafe};

    struct SilentOutput;
    impl crate::output::Output for SilentOutput {
        fn status(&self, _msg: &str) {}
        fn success(&self, _msg: &str) {}
        fn warn(&self, _msg: &str) {}
        fn error(&self, _err: &CohortError) {}
        fn result_json(&self, _data: &serde_json::Value) {}
        fn table(&self, _headers: &[&str], _rows: &[Vec<String>]) {}
        fn progress(&self, _total: u64, _label: &str) -> crate::output::Progress {
            crate::output::Progress::noop()
        }
    }

    #[test]
    fn raw_batch_size_matches_expected_formula() {
        // 1 MiB budget, 256 samples -> (1 MiB / 4) / (256*4+200) = ~217
        let got = raw_batch_size(256, 1 << 20);
        assert!(got > 200 && got < 230, "got {got}");
        // Zero samples => overhead-only denominator (200 bytes).
        assert_eq!(raw_batch_size(0, 800), 1);
    }

    #[test]
    fn drop_closes_parquet_footer_on_panic() {
        // Panicking inside a rayon worker must not leave a truncated, unreadable
        // parquet: the Drop impl writes the footer during unwind.
        let tmp = tempfile::tempdir().unwrap();
        let out = SilentOutput;
        let caught = panic::catch_unwind(AssertUnwindSafe(|| {
            let mut gw = GenotypeWriter::new(2, tmp.path(), 64 << 20).unwrap();
            gw.push("22", 15_000_000, "A", "G", "GT\t0/0\t0/1", 1, &out).unwrap();
            panic!("simulated worker panic before flush_all");
        }));
        assert!(caught.is_err(), "expected panic to propagate");

        let parquet_path = tmp.path().join("chromosome=22").join("data.parquet");
        assert!(parquet_path.exists(), "parquet missing at {}", parquet_path.display());
        let f = std::fs::File::open(&parquet_path).unwrap();
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(f);
        assert!(reader.is_ok(), "parquet footer missing: {:?}", reader.err());
    }

    /// noodles hands us raw samples text that starts with the FORMAT column
    /// ("GT" or "GT:AD:DP"). `push` must strip it before parsing per-sample
    /// genotypes; otherwise sample[0]'s slot gets the FORMAT bytes, everything
    /// shifts by one, and the last sample is dropped.
    #[test]
    fn push_strips_format_before_mapping_samples() {
        let tmp = tempfile::tempdir().unwrap();
        let out = SilentOutput;
        let mut gw = GenotypeWriter::new(3, tmp.path(), 64 << 20).unwrap();
        gw.push("22", 100, "A", "G", "GT\t0/0\t0/1\t1/1", 1, &out)
            .unwrap();
        assert_eq!(gw.dosages, vec![0.0, 1.0, 2.0]);
    }

    #[test]
    fn push_handles_multi_field_format() {
        // FORMAT="GT:AD:DP" is one tab-delimited field; sample values carry
        // their own colon separators inside each tab cell.
        let tmp = tempfile::tempdir().unwrap();
        let out = SilentOutput;
        let mut gw = GenotypeWriter::new(2, tmp.path(), 64 << 20).unwrap();
        gw.push(
            "22",
            200,
            "A",
            "G",
            "GT:AD:DP\t0/1:12,3:15\t1/1:0,20:20",
            1,
            &out,
        )
        .unwrap();
        assert_eq!(gw.dosages, vec![1.0, 2.0]);
    }

    #[test]
    fn push_treats_missing_samples_field_as_all_nan() {
        let tmp = tempfile::tempdir().unwrap();
        let out = SilentOutput;
        let mut gw = GenotypeWriter::new(2, tmp.path(), 64 << 20).unwrap();
        gw.push("22", 300, "A", "G", ".", 1, &out).unwrap();
        assert!(gw.dosages.iter().all(|d| d.is_nan()));
    }
}
