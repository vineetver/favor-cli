//! VCF → Parquet: streaming, columnar, bounded memory.
//!
//! Uses noodles-vcf for lazy field parsing + Arrow/Parquet for columnar writes.
//! Memory bound: one adaptive-sized batch of variants in Arrow columnar format at a time.
//! Batch size is derived from the caller's memory budget (see `derive_batch_size`).
//! Multi-allelic sites are split into biallelic records.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{Int32Builder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

use super::{RawRecord, VariantReader};
use crate::staar::genotype::GenotypeWriter;
use crate::store::list::VariantSetWriter;
use crate::error::CohortError;
use crate::output::Output;

/// Open a VCF (plain or BGZF-compressed) for buffered reading.
/// For BGZF files, `threads` decompression workers run in parallel.
pub fn open_vcf(path: &Path, threads: usize) -> Result<Box<dyn BufRead + Send>, CohortError> {
    let is_bgzf = path
        .extension()
        .map(|e| e == "gz" || e == "bgz")
        .unwrap_or(false);
    let file = File::open(path)
        .map_err(|e| CohortError::Input(format!("Cannot open '{}': {e}", path.display())))?;

    if is_bgzf {
        let workers = NonZeroUsize::new(threads.max(1))
            .expect("threads.max(1) >= 1");
        let bgzf = noodles_bgzf::reader::Builder::default()
            .set_worker_count(workers)
            .build_from_reader(file);
        Ok(Box::new(bgzf))
    } else {
        Ok(Box::new(BufReader::with_capacity(256 * 1024, file)))
    }
}

/// Bytes per variant in a batch: 7 columns × ~20 bytes average.
const BYTES_PER_VARIANT: u64 = 140;

/// Derive batch size from memory budget. VCF ingest is streaming — only one
/// batch + the Parquet writer live in memory at a time.
/// Floored at 10K (functional minimum), capped at 1M (diminishing returns).
fn derive_batch_size(memory_budget: u64) -> usize {
    let raw = memory_budget * 9 / 10 / BYTES_PER_VARIANT;
    (raw as usize).clamp(10_000, 1_000_000)
}

/// Canonical output schema for VCF ingest.
fn vcf_schema() -> Schema {
    Schema::new(vec![
        Field::new("chromosome", DataType::Utf8, false),
        Field::new("position", DataType::Int32, false),
        Field::new("ref", DataType::Utf8, false),
        Field::new("alt", DataType::Utf8, false),
        Field::new("rsid", DataType::Utf8, true),
        Field::new("qual", DataType::Utf8, true),
        Field::new("filter", DataType::Utf8, true),
    ])
}

/// Columnar batch builder with one contiguous array per field.
struct BatchBuilder {
    chromosome: StringBuilder,
    position: Int32Builder,
    ref_allele: StringBuilder,
    alt_allele: StringBuilder,
    rsid: StringBuilder,
    qual: StringBuilder,
    filter: StringBuilder,
    count: usize,
    batch_size: usize,
}

impl BatchBuilder {
    fn new(batch_size: usize) -> Self {
        Self {
            chromosome: StringBuilder::with_capacity(batch_size, batch_size * 4),
            position: Int32Builder::with_capacity(batch_size),
            ref_allele: StringBuilder::with_capacity(batch_size, batch_size * 4),
            alt_allele: StringBuilder::with_capacity(batch_size, batch_size * 4),
            rsid: StringBuilder::with_capacity(batch_size, batch_size * 12),
            qual: StringBuilder::with_capacity(batch_size, batch_size * 4),
            filter: StringBuilder::with_capacity(batch_size, batch_size * 8),
            count: 0,
            batch_size,
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn push(
        &mut self,
        chrom: &str,
        pos: i32,
        ref_allele: &str,
        alt_allele: &str,
        rsid: Option<&str>,
        qual: Option<&str>,
        filter: Option<&str>,
    ) {
        self.chromosome.append_value(chrom);
        self.position.append_value(pos);
        self.ref_allele.append_value(ref_allele);
        self.alt_allele.append_value(alt_allele);
        match rsid {
            Some(id) if id != "." => self.rsid.append_value(id),
            _ => self.rsid.append_null(),
        }
        match qual {
            Some(q) if q != "." => self.qual.append_value(q),
            _ => self.qual.append_null(),
        }
        match filter {
            Some(f) => self.filter.append_value(f),
            None => self.filter.append_null(),
        }
        self.count += 1;
    }

    fn is_full(&self) -> bool {
        self.count >= self.batch_size
    }

    fn len(&self) -> usize {
        self.count
    }

    fn finish(&mut self) -> Result<RecordBatch, CohortError> {
        let batch = RecordBatch::try_new(
            Arc::new(vcf_schema()),
            vec![
                Arc::new(self.chromosome.finish()),
                Arc::new(self.position.finish()),
                Arc::new(self.ref_allele.finish()),
                Arc::new(self.alt_allele.finish()),
                Arc::new(self.rsid.finish()),
                Arc::new(self.qual.finish()),
                Arc::new(self.filter.finish()),
            ],
        )
        .map_err(|e| CohortError::Analysis(format!("Arrow batch error: {e}")))?;

        self.count = 0;
        Ok(batch)
    }
}

/// Normalize a chromosome name: strip "chr" prefix, M→MT, 23→X, etc.
/// Returns None for non-standard contigs (filtered out).
pub fn normalize_chrom(raw: &str) -> Option<&'static str> {
    let stripped = if raw.len() >= 3 && raw[..3].eq_ignore_ascii_case("chr") {
        &raw[3..]
    } else {
        raw
    };
    match stripped {
        "1" => Some("1"),
        "2" => Some("2"),
        "3" => Some("3"),
        "4" => Some("4"),
        "5" => Some("5"),
        "6" => Some("6"),
        "7" => Some("7"),
        "8" => Some("8"),
        "9" => Some("9"),
        "10" => Some("10"),
        "11" => Some("11"),
        "12" => Some("12"),
        "13" => Some("13"),
        "14" => Some("14"),
        "15" => Some("15"),
        "16" => Some("16"),
        "17" => Some("17"),
        "18" => Some("18"),
        "19" => Some("19"),
        "20" => Some("20"),
        "21" => Some("21"),
        "22" => Some("22"),
        "X" | "x" => Some("X"),
        "Y" | "y" => Some("Y"),
        "M" | "MT" | "m" | "mt" => Some("MT"),
        "23" => Some("X"),
        "24" => Some("Y"),
        "25" => Some("MT"),
        _ => None,
    }
}

/// Parsimony normalize a biallelic variant to minimal VCF representation.
///
/// Implements the Tan et al. 2015 parsimony algorithm ("Unified representation
/// of genetic variants", Bioinformatics 31(13):2202-2204):
///
///   1. Right-trim: remove identical trailing bases from both alleles
///   2. Left-trim: remove identical leading bases from both alleles, adjusting POS
///
/// Both steps preserve at least 1 base in each allele (the VCF anchor base).
/// This is NOT full left-alignment (which requires a reference FASTA to shift
/// through tandem repeats) — it is the parsimony/minimal representation step.
///
/// The FAVOR annotation parquets (ref_vcf/alt_vcf) use this same representation,
/// so matching on the output of this function is correct for joins.
///
/// Examples:
///   SNV:       (G, A, 100)     → (G, A, 100)       — no change
///   Insertion: (A, ACGT, 100)  → (A, ACGT, 100)     — anchor base kept
///   Deletion:  (ACGT, A, 100)  → (ACGT, A, 100)     — anchor base kept
///   Padded:    (ACGT, AGT, 100)→ (AC, A, 100)       — right-trimmed GT, left C→del
///   MNV:       (ACG, TCA, 100) → (ACG, TCA, 100)    — no common prefix/suffix
///   Redundant: (AACG, AATG, 100) → (C, T, 102)     — trimmed AA prefix + G suffix
/// Zero allocation. Returns slices borrowing directly from the inputs.
/// VCF alleles are ASCII so byte-index slicing is always valid UTF-8.
pub fn parsimony_normalize<'a>(ref_allele: &'a str, alt_allele: &'a str, pos: i32) -> (&'a str, &'a str, i32) {
    let r = ref_allele.as_bytes();
    let a = alt_allele.as_bytes();
    let rlen = r.len();
    let alen = a.len();

    if rlen == 0 || alen == 0 || ref_allele == alt_allele || (rlen == 1 && alen == 1) {
        return (ref_allele, alt_allele, pos);
    }

    // Step 1: right-trim matching suffix, keep at least 1 base
    let max_suffix = rlen.min(alen) - 1;
    let mut suffix = 0;
    while suffix < max_suffix && r[rlen - 1 - suffix] == a[alen - 1 - suffix] {
        suffix += 1;
    }

    // Step 2: left-trim matching prefix after suffix removal
    let max_prefix = (rlen - suffix).min(alen - suffix) - 1;
    let mut prefix = 0;
    while prefix < max_prefix && r[prefix] == a[prefix] {
        prefix += 1;
    }

    (&ref_allele[prefix..rlen - suffix], &alt_allele[prefix..alen - suffix], pos + prefix as i32)
}

/// Uppercase ASCII into a reusable buffer. clear + push reuses the
/// existing heap allocation so the hot loop never calls the allocator.
#[inline]
pub(crate) fn ascii_uppercase_into(src: &str, buf: &mut String) {
    buf.clear();
    for &b in src.as_bytes() {
        buf.push(b.to_ascii_uppercase() as char);
    }
}

/// Owns a reusable `noodles_vcf::Record` so iteration does not allocate per
/// row. Field wrappers (`AlternateBases`, `Ids`, `Filters`, `Samples`) live
/// as locals inside `for_each` so their `&str` borrows stay alive across the
/// closure call.
pub struct VcfVariantReader {
    inner: noodles_vcf::io::Reader<Box<dyn BufRead + Send>>,
    record: noodles_vcf::Record,
    qual_buf: String,
    path: PathBuf,
}

impl VcfVariantReader {
    pub fn open(path: &Path, threads: usize) -> Result<Self, CohortError> {
        let buf = open_vcf(path, threads)?;
        let mut inner = noodles_vcf::io::Reader::new(buf);
        inner.read_header().map_err(|e| {
            CohortError::Input(format!("VCF header in '{}': {e}", path.display()))
        })?;
        Ok(Self {
            inner,
            record: noodles_vcf::Record::default(),
            qual_buf: String::with_capacity(32),
            path: path.to_path_buf(),
        })
    }
}

impl VariantReader for VcfVariantReader {
    fn sample_names(&mut self) -> Result<Vec<String>, CohortError> {
        // Re-reads the header from disk. Acceptable because the trait's
        // sample_names() is called at most a handful of times per file
        // (cohort-build sample lookup, sidecar write, header-consistency check)
        // and the current open Reader has already advanced past the header.
        crate::staar::genotype::read_sample_names(&self.path)
    }

    fn for_each(
        &mut self,
        f: &mut dyn for<'a> FnMut(RawRecord<'a>) -> Result<(), CohortError>,
    ) -> Result<(), CohortError> {
        loop {
            let n = self.inner.read_record(&mut self.record).map_err(|e| {
                CohortError::Analysis(format!(
                    "VCF parse error in {}: {e}",
                    self.path.display()
                ))
            })?;
            if n == 0 {
                return Ok(());
            }

            let position = match self.record.variant_start() {
                Some(Ok(p)) => p.get() as i32,
                _ => continue,
            };

            self.qual_buf.clear();
            let qual: Option<&str> = match self.record.quality_score() {
                Some(Ok(q)) => {
                    use std::fmt::Write;
                    let _ = write!(self.qual_buf, "{q}");
                    Some(self.qual_buf.as_str())
                }
                _ => None,
            };

            let alt_bases = self.record.alternate_bases();
            let alt_alleles: &str = alt_bases.as_ref();

            let ids = self.record.ids();
            let ids_raw: &str = ids.as_ref();
            let rsid = if ids_raw.is_empty() || ids_raw == "." {
                None
            } else {
                ids_raw.split(';').next()
            };

            let filters = self.record.filters();
            let filter_raw: &str = filters.as_ref();
            let filter = if filter_raw.is_empty() || filter_raw == "." {
                None
            } else {
                Some(filter_raw)
            };

            let samples = self.record.samples();
            let samples_text: &str = samples.as_ref();

            let rec = RawRecord {
                chromosome: self.record.reference_sequence_name(),
                position,
                ref_allele: self.record.reference_bases(),
                alt_alleles,
                rsid,
                qual,
                filter,
                samples_text,
            };
            f(rec)?;
        }
    }
}

/// Ingest result stats.
pub struct VcfIngestResult {
    pub variant_count: u64,
    pub filtered_contigs: u64,
    pub multiallelic_split: u64,
    pub genotype_variants: u64,
}

/// Per-chromosome writer state: batch + parquet writer + variant count.
struct ChromWriter {
    batch: BatchBuilder,
    writer: ArrowWriter<File>,
    count: u64,
}

/// Mutable state for the record-processing loop. Owns the reusable buffers,
/// per-chromosome writers, and counters. Used by both sequential and parallel paths.
/// `skip_variants` switches the genotype-only mode used by `extract_genotypes`:
/// variants.parquet writers are not created and cross-file chromosome-order
/// violations are surfaced as errors (because `GenotypeWriter` is single-chrom
/// streaming and cannot reopen a closed part).
struct RecordContext<'a> {
    ref_buf: String,
    alt_buf: String,
    writers: HashMap<&'static str, ChromWriter>,
    geno_writer: Option<&'a mut GenotypeWriter>,
    schema: Arc<Schema>,
    props: WriterProperties,
    batch_size: usize,
    output_dir: PathBuf,
    part_id: Option<usize>,
    chromosome_filter: Option<crate::types::ChromosomeSet>,
    variant_count: u64,
    filtered_contigs: u64,
    multiallelic_split: u64,
    skip_variants: bool,
    finished_chroms: HashSet<&'static str>,
}

impl<'a> RecordContext<'a> {
    fn new(
        geno_writer: Option<&'a mut GenotypeWriter>,
        memory_budget: u64,
        output_dir: PathBuf,
        part_id: Option<usize>,
        chromosome_filter: Option<crate::types::ChromosomeSet>,
    ) -> Self {
        Self::build(
            geno_writer,
            memory_budget,
            output_dir,
            part_id,
            chromosome_filter,
            false,
        )
    }

    /// Genotype-only mode: `variants.parquet` is not written. Used by
    /// `extract_genotypes`. Cross-file chromosome-order violations error out
    /// because the `GenotypeWriter` is single-chrom streaming.
    fn new_genotype_only(
        geno_writer: &'a mut GenotypeWriter,
        memory_budget: u64,
        output_dir: PathBuf,
        chromosome_filter: Option<crate::types::ChromosomeSet>,
    ) -> Self {
        Self::build(
            Some(geno_writer),
            memory_budget,
            output_dir,
            None,
            chromosome_filter,
            true,
        )
    }

    fn build(
        geno_writer: Option<&'a mut GenotypeWriter>,
        memory_budget: u64,
        output_dir: PathBuf,
        part_id: Option<usize>,
        chromosome_filter: Option<crate::types::ChromosomeSet>,
        skip_variants: bool,
    ) -> Self {
        let batch_size = derive_batch_size(memory_budget);
        let schema = Arc::new(vcf_schema());
        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(Default::default()))
            .set_max_row_group_row_count(Some(batch_size))
            .build();
        Self {
            ref_buf: String::with_capacity(256),
            alt_buf: String::with_capacity(256),
            writers: HashMap::new(),
            geno_writer,
            schema,
            props,
            batch_size,
            output_dir,
            part_id,
            chromosome_filter,
            variant_count: 0,
            filtered_contigs: 0,
            multiallelic_split: 0,
            skip_variants,
            finished_chroms: HashSet::new(),
        }
    }

    fn chrom_parquet_path(&self, chrom: &str) -> Result<PathBuf, CohortError> {
        let dir = self.output_dir.join(format!("chromosome={chrom}"));
        std::fs::create_dir_all(&dir)
            .map_err(|e| CohortError::Resource(format!("Cannot create '{}': {e}", dir.display())))?;
        let name = match self.part_id {
            Some(id) => format!("part_{id}.parquet"),
            None => "data.parquet".into(),
        };
        Ok(dir.join(name))
    }

    fn process(
        &mut self,
        rec: RawRecord<'_>,
        output: &dyn Output,
    ) -> Result<(), CohortError> {
        let chrom = match normalize_chrom(rec.chromosome) {
            Some(c) => c,
            None => {
                self.filtered_contigs += 1;
                return Ok(());
            }
        };
        if let Some(filt) = &self.chromosome_filter {
            if !filt.contains_canonical(chrom) {
                self.filtered_contigs += 1;
                return Ok(());
            }
        }

        // Genotype-only mode: the `GenotypeWriter` streams one chromosome
        // at a time. Re-opening a chromosome whose part file was already
        // closed would overwrite data, so surface it as an input error.
        if self.skip_variants {
            if let Some(ref gw) = self.geno_writer {
                if gw.current_chrom() != Some(chrom)
                    && self.finished_chroms.contains(chrom)
                {
                    return Err(CohortError::Input(format!(
                        "Chromosome {chrom} appears in the current input but \
                         its writer was already closed after processing an \
                         earlier file. Sort VCF files by chromosome or use \
                         one file per chromosome."
                    )));
                }
            }
        }

        let pos = rec.position;
        ascii_uppercase_into(rec.ref_allele, &mut self.ref_buf);

        let alt_str = rec.alt_alleles;
        let alt_count = if alt_str.is_empty() {
            0
        } else {
            alt_str.matches(',').count() + 1
        };
        if alt_count > 1 {
            self.multiallelic_split += alt_count as u64 - 1;
        }

        for (alt_idx, alt) in alt_str.split(',').enumerate() {
            ascii_uppercase_into(alt.trim(), &mut self.alt_buf);
            if self.alt_buf == "*"
                || self.alt_buf == "."
                || self.alt_buf.is_empty()
                || self.alt_buf.starts_with('<')
            {
                continue;
            }

            let (nr, na, np) = parsimony_normalize(&self.ref_buf, &self.alt_buf, pos);

            if !self.skip_variants {
                let cw = Self::get_or_create_writer(
                    chrom,
                    &mut self.writers,
                    &self.output_dir,
                    self.part_id,
                    &self.schema,
                    &self.props,
                    self.batch_size,
                    output,
                )?;
                cw.batch.push(chrom, np, nr, na, rec.rsid, rec.qual, rec.filter);
                cw.count += 1;

                if cw.batch.is_full() {
                    let rb = cw.batch.finish()?;
                    cw.writer
                        .write(&rb)
                        .map_err(|e| CohortError::Resource(format!("Parquet write error: {e}")))?;
                }
            }
            self.variant_count += 1;

            if let Some(ref mut gw) = self.geno_writer {
                gw.push(chrom, np, nr, na, rec.samples_text, (alt_idx + 1) as u8, output)?;
            }
        }
        Ok(())
    }

    /// End-of-file hook: in genotype-only mode, mark every chromosome the
    /// `GenotypeWriter` has touched except the current one as finished, so
    /// the next file reopening those triggers the `skip_variants` guard in
    /// `process`.
    fn finish_file(&mut self) {
        if !self.skip_variants {
            return;
        }
        let Some(ref gw) = self.geno_writer else { return };
        let current = gw.current_chrom();
        for c in gw.chromosomes() {
            if Some(c.as_str()) != current {
                if let Some(canon) = normalize_chrom(c) {
                    self.finished_chroms.insert(canon);
                }
            }
        }
    }

    fn drive(
        &mut self,
        reader: &mut dyn VariantReader,
        label: &str,
        output: &dyn Output,
    ) -> Result<(), CohortError> {
        let pb = output.progress(0, label);
        reader.for_each(&mut |rec| {
            self.process(rec, output)?;
            pb.inc(1);
            Ok(())
        })?;
        pb.finish(&format!("{} variants ingested", self.variant_count));
        Ok(())
    }

    #[allow(clippy::too_many_arguments)]
    fn get_or_create_writer<'w>(
        chrom: &'static str,
        writers: &'w mut HashMap<&'static str, ChromWriter>,
        output_dir: &Path,
        part_id: Option<usize>,
        schema: &Arc<Schema>,
        props: &WriterProperties,
        batch_size: usize,
        output: &dyn Output,
    ) -> Result<&'w mut ChromWriter, CohortError> {
        use std::collections::hash_map::Entry;
        match writers.entry(chrom) {
            Entry::Occupied(e) => Ok(e.into_mut()),
            Entry::Vacant(e) => {
                output.status(&format!("  chr{chrom}..."));
                let dir = output_dir.join(format!("chromosome={chrom}"));
                std::fs::create_dir_all(&dir).map_err(|err| {
                    CohortError::Resource(format!("Cannot create '{}': {err}", dir.display()))
                })?;
                let name = match part_id {
                    Some(id) => format!("part_{id}.parquet"),
                    None => "data.parquet".into(),
                };
                let path = dir.join(name);
                let f = File::create(&path).map_err(|err| {
                    CohortError::Resource(format!("Cannot create '{}': {err}", path.display()))
                })?;
                let w = ArrowWriter::try_new(f, schema.clone(), Some(props.clone()))
                    .map_err(|err| CohortError::Resource(format!("Parquet writer init: {err}")))?;
                Ok(e.insert(ChromWriter {
                    batch: BatchBuilder::new(batch_size),
                    writer: w,
                    count: 0,
                }))
            }
        }
    }

    /// Flush writers and register metadata into a VariantSetWriter (sequential path).
    fn flush_into_vs(
        mut self,
        vs_writer: &mut VariantSetWriter,
    ) -> Result<VcfIngestResult, CohortError> {
        let genotype_variants = self.geno_writer.as_ref().map_or(0, |g| g.variant_count());

        vs_writer.set_columns(self.schema.fields().iter().map(|f| f.name().clone()).collect());
        let writers: Vec<_> = self.writers.drain().collect();
        for (chrom, mut cw) in writers {
            if cw.batch.len() > 0 {
                let rb = cw.batch.finish()?;
                cw.writer
                    .write(&rb)
                    .map_err(|e| CohortError::Resource(format!("Parquet write error: {e}")))?;
            }
            cw.writer
                .close()
                .map_err(|e| CohortError::Resource(format!("Parquet close error: {e}")))?;
            let path = self.chrom_parquet_path(chrom)?;
            let size = std::fs::metadata(&path).map_or(0, |m| m.len());
            vs_writer.register_chrom(chrom, cw.count, size);
        }

        Ok(VcfIngestResult {
            variant_count: self.variant_count,
            filtered_contigs: self.filtered_contigs,
            multiallelic_split: self.multiallelic_split,
            genotype_variants,
        })
    }

    /// Flush and close writers without metadata registration (parallel path).
    fn flush(self) -> Result<VcfIngestResult, CohortError> {
        let genotype_variants = self.geno_writer.as_ref().map_or(0, |g| g.variant_count());

        for (_chrom, mut cw) in self.writers {
            if cw.batch.len() > 0 {
                let rb = cw.batch.finish()?;
                cw.writer
                    .write(&rb)
                    .map_err(|e| CohortError::Resource(format!("Parquet write error: {e}")))?;
            }
            cw.writer
                .close()
                .map_err(|e| CohortError::Resource(format!("Parquet close error: {e}")))?;
        }

        Ok(VcfIngestResult {
            variant_count: self.variant_count,
            filtered_contigs: self.filtered_contigs,
            multiallelic_split: self.multiallelic_split,
            genotype_variants,
        })
    }

    fn ingest_files(
        &mut self,
        files: &[impl AsRef<Path>],
        handler: &dyn super::format::FormatHandler,
        threads: usize,
        output: &dyn Output,
    ) -> Result<(), CohortError> {
        for (file_idx, input_path) in files.iter().enumerate() {
            let input_path = input_path.as_ref();
            if files.len() > 1 {
                output.status(&format!(
                    "  File {}/{}: {}",
                    file_idx + 1,
                    files.len(),
                    input_path.file_name().unwrap_or_default().to_string_lossy()
                ));
            }
            let mut reader = handler.open_reader(input_path, threads)?;
            let label = format!(
                "ingesting {}",
                input_path.file_name().unwrap_or_default().to_string_lossy()
            );
            self.drive(reader.as_mut(), &label, output)?;
            self.finish_file();
        }
        Ok(())
    }
}

/// Stream genotypes into an existing `GenotypeWriter`. No `variants.parquet`
/// is written — this is the standalone extraction path used by the cohort
/// store when variants have already been ingested separately.
pub fn stream_genotypes(
    paths: &[PathBuf],
    handler: &(dyn super::format::FormatHandler + Send + Sync),
    gw: &mut GenotypeWriter,
    memory_budget: u64,
    threads: usize,
    geno_dir: PathBuf,
    output: &dyn Output,
) -> Result<(), CohortError> {
    let mut ctx = RecordContext::new_genotype_only(gw, memory_budget, geno_dir, None);
    ctx.ingest_files(paths, handler, threads, output)
}

/// Stream one or more variant files to per-chromosome parquet (sequential).
/// Format-agnostic: caller picks the handler (VCF, GDS, ...).
#[allow(clippy::too_many_arguments)]
pub fn ingest_vcfs(
    input_paths: &[PathBuf],
    handler: &(dyn super::format::FormatHandler + Send + Sync),
    vs_writer: &mut VariantSetWriter,
    geno_writer: Option<&mut GenotypeWriter>,
    memory_budget: u64,
    threads: usize,
    chromosome_filter: Option<&crate::types::ChromosomeSet>,
    output: &dyn Output,
) -> Result<VcfIngestResult, CohortError> {
    let mut ctx = RecordContext::new(
        geno_writer, memory_budget,
        vs_writer.root().to_path_buf(), None,
        chromosome_filter.cloned(),
    );
    output.status(&format!(
        "  Batch size: {} variants/chrom ({:.1}G memory)",
        ctx.batch_size,
        memory_budget as f64 / (1024.0 * 1024.0 * 1024.0)
    ));

    ctx.ingest_files(input_paths, handler, threads, output)?;
    ctx.flush_into_vs(vs_writer)
}

/// Validate that all input files report identical sample lists. Routes
/// through the format handler so this works for VCF (header parse) and
/// GDS (FFI read of /sample.id) alike.
fn validate_sample_consistency(
    paths: &[PathBuf],
    handler: &(dyn super::format::FormatHandler + Send + Sync),
    n_samples: usize,
) -> Result<(), CohortError> {
    use rayon::prelude::*;

    if paths.len() <= 1 || n_samples == 0 {
        return Ok(());
    }

    let first_samples = handler.open_reader(&paths[0], 1)?.sample_names()?;
    paths[1..].par_iter().try_for_each(|p| {
        let samples = handler.open_reader(p, 1)?.sample_names()?;
        if samples.len() != first_samples.len() {
            return Err(CohortError::Input(format!(
                "Sample count mismatch: '{}' has {} samples but '{}' has {}",
                paths[0].display(), first_samples.len(),
                p.display(), samples.len(),
            )));
        }
        if samples != first_samples {
            return Err(CohortError::Input(format!(
                "Sample order mismatch between '{}' and '{}'. \
                 All input files must have identical samples in the same order.",
                paths[0].display(), p.display(),
            )));
        }
        Ok(())
    })
}

/// Parallel multi-file ingest. Files are chunked across N workers, each
/// worker processes its chunk sequentially with a single RecordContext.
/// N is bounded by thread count. Validates sample consistency before
/// spawning. The handler is passed to each worker so this works for any
/// VariantReader-backed format (VCF, GDS, ...).
#[allow(clippy::too_many_arguments)]
pub fn ingest_vcfs_parallel(
    input_paths: &[PathBuf],
    handler: &(dyn super::format::FormatHandler + Send + Sync),
    output_dir: &Path,
    geno_dir: Option<&Path>,
    n_samples: usize,
    memory_budget: u64,
    threads: usize,
    chromosome_filter: Option<&crate::types::ChromosomeSet>,
    output: &dyn Output,
) -> Result<VcfIngestResult, CohortError> {
    use rayon::prelude::*;

    let threads = preflight(
        input_paths, output_dir, geno_dir, n_samples,
        memory_budget, threads, output,
    )?;

    if n_samples > 0 {
        output.status("  Validating sample consistency across files...");
        validate_sample_consistency(input_paths, handler, n_samples)?;
    }

    let n_workers = input_paths.len().min(threads);
    let memory_per_worker = memory_budget / n_workers as u64;

    let chunks: Vec<Vec<&PathBuf>> = (0..n_workers)
        .map(|w| input_paths.iter().skip(w).step_by(n_workers).collect())
        .collect();

    output.status(&format!(
        "  Parallel: {} files, {} workers, {:.1}G/worker",
        input_paths.len(), n_workers,
        memory_per_worker as f64 / (1024.0 * 1024.0 * 1024.0),
    ));

    std::fs::create_dir_all(output_dir).map_err(|e| {
        CohortError::Resource(format!("Cannot create '{}': {e}", output_dir.display()))
    })?;
    if let Some(gd) = geno_dir {
        std::fs::create_dir_all(gd).map_err(|e| {
            CohortError::Resource(format!("Cannot create '{}': {e}", gd.display()))
        })?;
    }

    let results: Vec<VcfIngestResult> = chunks
        .into_par_iter()
        .enumerate()
        .map(|(worker_id, file_chunk)| {
            run_worker(
                worker_id, &file_chunk, handler, output_dir, geno_dir,
                n_samples, memory_per_worker,
                chromosome_filter.cloned(),
                output,
            )
            .map_err(|e| e.with_context(format!(
                "worker {worker_id} ({})",
                file_chunk.iter()
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
            )))
        })
        .collect::<Result<Vec<_>, _>>()?;

    let mut total = VcfIngestResult {
        variant_count: 0, filtered_contigs: 0,
        multiallelic_split: 0, genotype_variants: 0,
    };
    for r in results {
        total.variant_count += r.variant_count;
        total.filtered_contigs += r.filtered_contigs;
        total.multiallelic_split += r.multiallelic_split;
        total.genotype_variants += r.genotype_variants;
    }
    Ok(total)
}

/// One rayon worker: own a `GenotypeWriter` and a `RecordContext`, ingest the
/// file chunk sequentially, close both cleanly. Extracted so the caller can
/// attach `(worker_id, file paths)` to any error via `with_context`.
#[allow(clippy::too_many_arguments)]
fn run_worker(
    worker_id: usize,
    file_chunk: &[&PathBuf],
    handler: &(dyn super::format::FormatHandler + Send + Sync),
    output_dir: &Path,
    geno_dir: Option<&Path>,
    n_samples: usize,
    memory_per_worker: u64,
    chromosome_filter: Option<crate::types::ChromosomeSet>,
    output: &dyn Output,
) -> Result<VcfIngestResult, CohortError> {
    let mut gw = match geno_dir {
        Some(gd) => Some(GenotypeWriter::with_part_id(
            n_samples, gd, memory_per_worker, Some(worker_id),
        )?),
        None => None,
    };
    let mut ctx = RecordContext::new(
        gw.as_mut(), memory_per_worker,
        output_dir.to_path_buf(), Some(worker_id),
        chromosome_filter,
    );
    ctx.ingest_files(file_chunk, handler, 1, output)?;
    let result = ctx.flush()?;
    if let Some(mut g) = gw {
        g.flush_all()?;
    }
    Ok(result)
}

/// Pre-flight checks for parallel VCF ingest. Catches duplicate input paths,
/// stale part files from a prior failed run, insufficient fd headroom, and
/// memory budgets that would starve the genotype writer. Returns the
/// (possibly reduced) worker count that is safe to spawn.
fn preflight(
    input_paths: &[PathBuf],
    output_dir: &Path,
    geno_dir: Option<&Path>,
    n_samples: usize,
    memory_budget: u64,
    requested_threads: usize,
    output: &dyn Output,
) -> Result<usize, CohortError> {
    dedup_canonical(input_paths)?;
    let mut roots: Vec<&Path> = vec![output_dir];
    if let Some(gd) = geno_dir { roots.push(gd); }
    clean_stale_parts(&roots, output)?;
    let threads = cap_threads_for_fds(requested_threads, output)?;
    let threads = cap_threads_for_batch_size(
        n_samples, memory_budget, threads, geno_dir.is_some(),
    )?;
    Ok(threads)
}

fn dedup_canonical(paths: &[PathBuf]) -> Result<(), CohortError> {
    use std::collections::HashSet;
    let mut seen = HashSet::with_capacity(paths.len());
    for p in paths {
        let canon = std::fs::canonicalize(p).map_err(|e| {
            CohortError::Input(format!("cannot open '{}': {e}", p.display()))
        })?;
        if !seen.insert(canon) {
            return Err(CohortError::Input(format!(
                "duplicate VCF input: '{}' (each file must appear at most once)",
                p.display(),
            )));
        }
    }
    Ok(())
}

fn clean_stale_parts(roots: &[&Path], output: &dyn Output) -> Result<(), CohortError> {
    let mut removed = 0usize;
    for root in roots {
        if !root.exists() { continue; }
        // Parts live at root/chromosome=<chr>/part_<id>.parquet (depth 2).
        for entry in walkdir::WalkDir::new(root).min_depth(2).max_depth(2) {
            let entry = entry.map_err(|e| CohortError::Resource(format!(
                "cannot scan '{}': {e}", root.display(),
            )))?;
            let name = entry.file_name().to_string_lossy();
            if name.starts_with("part_") && name.ends_with(".parquet") {
                std::fs::remove_file(entry.path()).map_err(|e| CohortError::Resource(format!(
                    "cannot remove stale part '{}': {e}", entry.path().display(),
                )))?;
                removed += 1;
            }
        }
    }
    if removed > 0 {
        output.status(&format!("  Preflight: cleaned {removed} stale part file(s)"));
    }
    Ok(())
}

#[cfg(unix)]
fn fd_soft_limit() -> Option<u64> {
    let mut rlim = libc::rlimit { rlim_cur: 0, rlim_max: 0 };
    // SAFETY: `rlim` is a valid local; RLIMIT_NOFILE is the standard descriptor resource.
    let rc = unsafe { libc::getrlimit(libc::RLIMIT_NOFILE, &mut rlim) };
    if rc == 0 { Some(rlim.rlim_cur) } else { None }
}
#[cfg(not(unix))]
fn fd_soft_limit() -> Option<u64> { None }

// Per-worker fd cost: ~24 chromosome parquet writers + bgzf reader +
// genotype writer + slack. Baseline covers stdio, mmap handles, etc.
// Usable fraction is 80% of the soft limit — leaves headroom for lib internals.
fn fd_safe_workers(limit: u64) -> usize {
    const FDS_PER_WORKER: u64 = 30;
    const BASELINE_FDS: u64 = 64;
    // `limit / 10 * 8` avoids overflow at the u64::MAX edge cases exercised by tests.
    let budget = (limit / 10 * 8).saturating_sub(BASELINE_FDS);
    (budget / FDS_PER_WORKER) as usize
}

fn cap_threads_for_fds(threads: usize, output: &dyn Output) -> Result<usize, CohortError> {
    let Some(limit) = fd_soft_limit() else { return Ok(threads); };
    let safe = fd_safe_workers(limit);
    if safe == 0 {
        return Err(CohortError::Resource(format!(
            "fd soft limit {limit} too low for parallel ingest; \
             raise with `ulimit -n 4096` and retry",
        )));
    }
    if threads <= safe { return Ok(threads); }
    output.status(&format!(
        "  Preflight: capping workers {threads} -> {safe} (fd soft limit {limit})",
    ));
    Ok(safe)
}

fn cap_threads_for_batch_size(
    n_samples: usize,
    memory_budget: u64,
    threads: usize,
    writes_genotypes: bool,
) -> Result<usize, CohortError> {
    // Viability floor for the GenotypeWriter row-group. Below this the writer
    // flushes per variant, killing throughput.
    const MIN_VIABLE_BATCH: u64 = 500;
    if !writes_genotypes || n_samples == 0 { return Ok(threads); }
    // Find the largest worker count where per-worker raw capacity still clears MIN.
    let mut t = threads.max(1);
    while t > 1 && crate::staar::genotype::raw_batch_size(n_samples, memory_budget / t as u64)
        < MIN_VIABLE_BATCH
    {
        t -= 1;
    }
    if crate::staar::genotype::raw_batch_size(n_samples, memory_budget / t as u64)
        < MIN_VIABLE_BATCH
    {
        return Err(CohortError::Resource(format!(
            "memory budget {:.1}G too small for {n_samples} samples \
             (even 1 worker can't batch {MIN_VIABLE_BATCH} variants); \
             raise --memory or reduce sample count",
            memory_budget as f64 / (1u64 << 30) as f64,
        )));
    }
    Ok(t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn snv() { assert_eq!(parsimony_normalize("G", "A", 100), ("G", "A", 100)); }
    #[test]
    fn insertion() { assert_eq!(parsimony_normalize("A", "ACGT", 100), ("A", "ACGT", 100)); }
    #[test]
    fn deletion() { assert_eq!(parsimony_normalize("ACGT", "A", 100), ("ACGT", "A", 100)); }
    #[test]
    fn suffix_trim() { assert_eq!(parsimony_normalize("ACGTG", "AG", 100), ("ACGT", "A", 100)); }
    #[test]
    fn prefix_and_suffix() { assert_eq!(parsimony_normalize("AACG", "AATG", 100), ("C", "T", 102)); }
    #[test]
    fn mnv_no_trim() { assert_eq!(parsimony_normalize("ACG", "TCA", 100), ("ACG", "TCA", 100)); }
    #[test]
    fn identical() { assert_eq!(parsimony_normalize("A", "A", 100), ("A", "A", 100)); }
    #[test]
    fn single_base_ins() { assert_eq!(parsimony_normalize("T", "TC", 100), ("T", "TC", 100)); }

    #[test]
    fn borrows_input_directly() {
        let r = "AACG";
        let a = "AATG";
        let (nr, na, _) = parsimony_normalize(r, a, 100);
        assert!(std::ptr::eq(nr.as_ptr(), r[2..3].as_ptr()));
        assert!(std::ptr::eq(na.as_ptr(), a[2..3].as_ptr()));
    }

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

    fn touch(path: &Path) {
        std::fs::write(path, b"").unwrap();
    }

    #[test]
    fn dedup_canonical_rejects_same_file_twice() {
        let tmp = tempfile::tempdir().unwrap();
        let p = tmp.path().join("a.vcf.gz");
        touch(&p);
        let err = dedup_canonical(&[p.clone(), p.clone()]).unwrap_err();
        assert!(matches!(err, CohortError::Input(ref m) if m.contains("duplicate")), "got {err:?}");
    }

    #[test]
    fn dedup_canonical_accepts_distinct_files() {
        let tmp = tempfile::tempdir().unwrap();
        let a = tmp.path().join("a.vcf.gz");
        let b = tmp.path().join("b.vcf.gz");
        touch(&a);
        touch(&b);
        assert!(dedup_canonical(&[a, b]).is_ok());
    }

    #[test]
    fn dedup_canonical_resolves_symlinks_to_same_target() {
        // Two distinct paths both pointing at the same file should be flagged.
        let tmp = tempfile::tempdir().unwrap();
        let real = tmp.path().join("real.vcf.gz");
        let link = tmp.path().join("link.vcf.gz");
        touch(&real);
        #[cfg(unix)]
        std::os::unix::fs::symlink(&real, &link).unwrap();
        #[cfg(not(unix))]
        touch(&link); // can't symlink without privileges; plain touch still makes two distinct paths
        let err = dedup_canonical(&[real, link]);
        #[cfg(unix)]
        assert!(matches!(err, Err(CohortError::Input(_))));
        #[cfg(not(unix))]
        assert!(err.is_ok()); // distinct files, not actually the same
    }

    #[test]
    fn clean_stale_parts_removes_parts_only() {
        let tmp = tempfile::tempdir().unwrap();
        let chr_dir = tmp.path().join("chromosome=22");
        std::fs::create_dir_all(&chr_dir).unwrap();
        touch(&chr_dir.join("part_0.parquet"));
        touch(&chr_dir.join("part_1.parquet"));
        touch(&chr_dir.join("data.parquet"));       // not a part file
        touch(&chr_dir.join("other.txt"));          // unrelated
        clean_stale_parts(&[tmp.path()], &SilentOutput).unwrap();
        assert!(!chr_dir.join("part_0.parquet").exists());
        assert!(!chr_dir.join("part_1.parquet").exists());
        assert!(chr_dir.join("data.parquet").exists());
        assert!(chr_dir.join("other.txt").exists());
    }

    #[test]
    fn clean_stale_parts_missing_root_is_ok() {
        let tmp = tempfile::tempdir().unwrap();
        let missing = tmp.path().join("does-not-exist");
        assert!(clean_stale_parts(&[&missing], &SilentOutput).is_ok());
    }

    #[test]
    fn fd_safe_workers_scales_with_limit() {
        assert_eq!(fd_safe_workers(64), 0);         // below baseline
        assert!(fd_safe_workers(1024) >= 20);       // default ulimit gives many workers
        assert!(fd_safe_workers(4096) > fd_safe_workers(1024));
        assert!(fd_safe_workers(u64::MAX) > 0);     // no overflow
    }

    #[test]
    fn cap_threads_for_batch_size_no_genotypes_passthrough() {
        // No genotype writer => no viability check, threads unchanged.
        let got = cap_threads_for_batch_size(0, 1 << 20, 16, false).unwrap();
        assert_eq!(got, 16);
    }

    #[test]
    fn cap_threads_for_batch_size_reduces_when_memory_low() {
        // 200k samples, 16 GiB budget — enough for a handful of workers, not 64.
        let capped = cap_threads_for_batch_size(200_000, 16u64 << 30, 64, true).unwrap();
        assert!((1..64).contains(&capped), "expected cap, got {capped}");
    }

    #[test]
    fn cap_threads_for_batch_size_errors_on_starved_single_worker() {
        // 1 MiB for 200k samples — even 1 worker can't batch 500 variants.
        let err = cap_threads_for_batch_size(200_000, 1 << 20, 1, true).unwrap_err();
        assert!(matches!(err, CohortError::Resource(_)), "got {err:?}");
    }

    #[test]
    fn cap_threads_for_batch_size_passes_when_memory_ample() {
        // 64 GiB, 10k samples, 8 workers — all fit comfortably.
        let got = cap_threads_for_batch_size(10_000, 64u64 << 30, 8, true).unwrap();
        assert_eq!(got, 8);
    }

    #[test]
    fn error_with_context_preserves_variant() {
        let e = CohortError::Input("bad".into()).with_context("worker 3");
        assert!(matches!(e, CohortError::Input(ref m) if m == "worker 3: bad"));
        let e = CohortError::Resource("hm".into()).with_context("worker 3");
        assert!(matches!(e, CohortError::Resource(ref m) if m == "worker 3: hm"));
    }
}
