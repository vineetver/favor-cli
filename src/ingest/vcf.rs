//! VCF → Parquet: streaming, columnar, bounded memory.
//!
//! Uses noodles-vcf for lazy field parsing + Arrow/Parquet for columnar writes.
//! Memory bound: one adaptive-sized batch of variants in Arrow columnar format at a time.
//! Batch size is derived from the caller's memory budget (see `derive_batch_size`).
//! Multi-allelic sites are split into biallelic records.

use std::collections::HashMap;
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
/// batch + the Parquet writer live in memory at a time, so 50% goes to the batch.
/// Floored at 10K (functional minimum), capped at 1M (diminishing returns).
fn derive_batch_size(memory_budget: u64) -> usize {
    let raw = memory_budget / 2 / BYTES_PER_VARIANT;
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
fn ascii_uppercase_into(src: &str, buf: &mut String) {
    buf.clear();
    for &b in src.as_bytes() {
        buf.push(b.to_ascii_uppercase() as char);
    }
}

/// Ingest result stats.
pub struct VcfIngestResult {
    pub variant_count: u64,
    pub filtered_contigs: u64,
    pub multiallelic_split: u64,
    /// When single-pass ingest is active, the number of genotype variants written.
    pub genotype_variants: u64,
}

/// Per-chromosome writer state: batch + parquet writer + variant count.
struct ChromWriter {
    batch: BatchBuilder,
    writer: ArrowWriter<File>,
    count: u64,
}

/// Stream one or more VCF files to per-chromosome parquet files.
/// All files share the same set of per-chromosome writers, so blocks from
/// different files (e.g., UKB b0-b22) merge into the correct chromosome output.
/// Bounded memory: 25 chromosome writers max, each with adaptive batch size.
///
/// When `geno_writer` is `Some`, genotype dosages are extracted in the same
/// pass and written to the GenotypeWriter (single-pass ingest).
pub fn ingest_vcfs(
    input_paths: &[PathBuf],
    vs_writer: &mut VariantSetWriter,
    geno_writer: Option<&mut GenotypeWriter>,
    memory_budget: u64,
    threads: usize,
    output: &dyn Output,
) -> Result<VcfIngestResult, CohortError> {
    let batch_size = derive_batch_size(memory_budget);
    output.status(&format!(
        "  Batch size: {} variants/chrom ({:.1}G memory)",
        batch_size,
        memory_budget as f64 / (1024.0 * 1024.0 * 1024.0)
    ));

    let schema = Arc::new(vcf_schema());
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_row_count(Some(batch_size))
        .build();

    let mut writers: HashMap<&'static str, ChromWriter> = HashMap::new();
    let mut variant_count: u64 = 0;
    let mut filtered_contigs: u64 = 0;
    let mut multiallelic_split: u64 = 0;
    let mut gw = geno_writer;
    let mut ref_buf = String::with_capacity(256);
    let mut alt_buf = String::with_capacity(256);

    for (file_idx, input_path) in input_paths.iter().enumerate() {
        if input_paths.len() > 1 {
            output.status(&format!(
                "  File {}/{}: {}",
                file_idx + 1,
                input_paths.len(),
                input_path.file_name().unwrap_or_default().to_string_lossy()
            ));
        }

        let reader = open_vcf(input_path, threads)?;
        let mut vcf_reader = noodles_vcf::io::Reader::new(reader);
        {
            let mut hr = vcf_reader.header_reader();
            std::io::copy(&mut hr, &mut std::io::sink()).map_err(|e| {
                CohortError::Input(format!(
                    "Skip VCF header in {}: {e}",
                    input_path.display()
                ))
            })?;
        }

        let pb = output.progress(
            0,
            &format!(
                "ingesting {}",
                input_path
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
            ),
        );
        for result in vcf_reader.records() {
            let record = result.map_err(|e| {
                CohortError::Analysis(format!("VCF parse error in {}: {e}", input_path.display()))
            })?;

            process_record(
                &record,
                &mut ref_buf,
                &mut alt_buf,
                &mut writers,
                vs_writer,
                gw.as_deref_mut(),
                &schema,
                &props,
                batch_size,
                &mut variant_count,
                &mut filtered_contigs,
                &mut multiallelic_split,
                output,
            )?;
            pb.inc(1);
        }
        pb.finish(&format!("{variant_count} variants ingested"));
    }

    let genotype_variants = gw.as_ref().map_or(0, |g| g.variant_count());

    // Flush and close all writers
    vs_writer.set_columns(schema.fields().iter().map(|f| f.name().clone()).collect());
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
        let path = vs_writer.chrom_path(chrom)?;
        let size = std::fs::metadata(&path).map_or(0, |m| m.len());
        vs_writer.register_chrom(chrom, cw.count, size);
    }

    Ok(VcfIngestResult {
        variant_count,
        filtered_contigs,
        multiallelic_split,
        genotype_variants,
    })
}

fn get_or_create_writer<'a>(
    chrom: &'static str,
    writers: &'a mut HashMap<&'static str, ChromWriter>,
    vs_writer: &VariantSetWriter,
    schema: &Arc<Schema>,
    props: &WriterProperties,
    batch_size: usize,
    output: &dyn Output,
) -> Result<&'a mut ChromWriter, CohortError> {
    use std::collections::hash_map::Entry;
    match writers.entry(chrom) {
        Entry::Occupied(e) => Ok(e.into_mut()),
        Entry::Vacant(e) => {
            output.status(&format!("  chr{chrom}..."));
            let path = vs_writer.chrom_path(chrom)?;
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

/// Process one VCF record. Zero per-record heap allocations: ref/alt
/// uppercased into caller-owned buffers, parsimony_normalize borrows
/// from those buffers, filter/samples borrowed from the noodles Record.
#[allow(clippy::too_many_arguments)]
fn process_record(
    record: &noodles_vcf::Record,
    ref_buf: &mut String,
    alt_buf: &mut String,
    writers: &mut HashMap<&'static str, ChromWriter>,
    vs_writer: &VariantSetWriter,
    mut geno_writer: Option<&mut GenotypeWriter>,
    schema: &Arc<Schema>,
    props: &WriterProperties,
    batch_size: usize,
    variant_count: &mut u64,
    filtered_contigs: &mut u64,
    multiallelic_split: &mut u64,
    output: &dyn Output,
) -> Result<(), CohortError> {
    let chrom = match normalize_chrom(record.reference_sequence_name()) {
        Some(c) => c,
        None => { *filtered_contigs += 1; return Ok(()); }
    };
    let pos = match record.variant_start() {
        Some(Ok(p)) => p.get() as i32,
        _ => return Ok(()),
    };

    ascii_uppercase_into(record.reference_bases(), ref_buf);

    let alt_bases = record.alternate_bases();
    let alt_str: &str = alt_bases.as_ref();

    let ids_wrapper = record.ids();
    let ids_str: &str = ids_wrapper.as_ref();
    let rsid = if ids_str.is_empty() || ids_str == "." {
        None
    } else {
        ids_str.split(';').next()
    };

    let qual_str = match record.quality_score() {
        Some(Ok(q)) => Some(format!("{q}")),
        _ => None,
    };

    let filters_wrapper = record.filters();
    let filter_raw: &str = filters_wrapper.as_ref();
    let filter_str: Option<&str> = if filter_raw.is_empty() || filter_raw == "." {
        None
    } else {
        Some(filter_raw)
    };

    let alt_count = if alt_str.is_empty() { 0 } else { alt_str.matches(',').count() + 1 };
    if alt_count > 1 {
        *multiallelic_split += alt_count as u64 - 1;
    }

    let samples_raw;
    let samples_str: &str = if geno_writer.is_some() {
        samples_raw = record.samples();
        samples_raw.as_ref()
    } else {
        ""
    };

    for (alt_idx, alt) in alt_str.split(',').enumerate() {
        ascii_uppercase_into(alt.trim(), alt_buf);
        if alt_buf == "*" || alt_buf == "." || alt_buf.is_empty() || alt_buf.starts_with('<') {
            continue;
        }

        let (nr, na, np) = parsimony_normalize(ref_buf, alt_buf, pos);

        let cw =
            get_or_create_writer(chrom, writers, vs_writer, schema, props, batch_size, output)?;

        cw.batch.push(chrom, np, nr, na, rsid, qual_str.as_deref(), filter_str);
        cw.count += 1;
        *variant_count += 1;

        if cw.batch.is_full() {
            let rb = cw.batch.finish()?;
            cw.writer
                .write(&rb)
                .map_err(|e| CohortError::Resource(format!("Parquet write error: {e}")))?;
        }

        if let Some(ref mut gw) = geno_writer {
            gw.push(chrom, np, nr, na, samples_str, (alt_idx + 1) as u8, output)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::parsimony_normalize;

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
}
