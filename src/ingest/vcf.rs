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

use crate::data::VariantSetWriter;
use crate::error::FavorError;
use crate::output::Output;

/// Open a VCF (plain or BGZF-compressed) for buffered reading.
/// For BGZF files, `threads` decompression workers run in parallel.
pub fn open_vcf(path: &Path, threads: usize) -> Result<Box<dyn BufRead + Send>, FavorError> {
    let is_bgzf = path
        .extension()
        .map(|e| e == "gz" || e == "bgz")
        .unwrap_or(false);
    let file = File::open(path)
        .map_err(|e| FavorError::Input(format!("Cannot open '{}': {e}", path.display())))?;

    if is_bgzf {
        // threads.max(1) is always >= 1, so NonZeroUsize::new cannot fail
        let workers = NonZeroUsize::new(threads.max(1)).unwrap();
        let bgzf = noodles_bgzf::reader::Builder::default()
            .set_worker_count(workers)
            .build_from_reader(file);
        // noodles_bgzf::Reader implements BufRead — do not double-buffer.
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

    fn finish(&mut self) -> Result<RecordBatch, FavorError> {
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
        .map_err(|e| FavorError::Analysis(format!("Arrow batch error: {e}")))?;

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
pub fn parsimony_normalize(ref_allele: &str, alt_allele: &str, pos: i32) -> (String, String, i32) {
    let r_bytes = ref_allele.as_bytes();
    let a_bytes = alt_allele.as_bytes();

    // Edge case: empty or identical alleles — return as-is
    if r_bytes.is_empty() || a_bytes.is_empty() || ref_allele == alt_allele {
        return (ref_allele.to_string(), alt_allele.to_string(), pos);
    }

    // SNV fast path: single base each, no trimming possible
    if r_bytes.len() == 1 && a_bytes.len() == 1 {
        return (ref_allele.to_string(), alt_allele.to_string(), pos);
    }

    let rlen = r_bytes.len();
    let alen = a_bytes.len();

    // Step 1: Right-trim — count matching suffix bases
    // Keep at least 1 base in each allele
    let max_suffix = (rlen.min(alen)) - 1; // leave at least 1 base
    let mut suffix_trim = 0;
    while suffix_trim < max_suffix
        && r_bytes[rlen - 1 - suffix_trim] == a_bytes[alen - 1 - suffix_trim]
    {
        suffix_trim += 1;
    }

    // Step 2: Left-trim — count matching prefix bases (after suffix removal)
    // Keep at least 1 base in each allele after both trims
    let r_after_suffix = rlen - suffix_trim;
    let a_after_suffix = alen - suffix_trim;
    let max_prefix = (r_after_suffix.min(a_after_suffix)) - 1;
    let mut prefix_trim = 0;
    while prefix_trim < max_prefix && r_bytes[prefix_trim] == a_bytes[prefix_trim] {
        prefix_trim += 1;
    }

    // Build result
    let norm_ref = &r_bytes[prefix_trim..rlen - suffix_trim];
    let norm_alt = &a_bytes[prefix_trim..alen - suffix_trim];
    let norm_pos = pos + prefix_trim as i32;

    (
        String::from_utf8_lossy(norm_ref).to_string(),
        String::from_utf8_lossy(norm_alt).to_string(),
        norm_pos,
    )
}

/// Ingest result stats.
pub struct VcfIngestResult {
    pub variant_count: u64,
    pub filtered_contigs: u64,
    pub multiallelic_split: u64,
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
pub fn ingest_vcfs(
    input_paths: &[PathBuf],
    vs_writer: &mut VariantSetWriter,
    memory_budget: u64,
    threads: usize,
    output: &dyn Output,
) -> Result<VcfIngestResult, FavorError> {
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
        let header = vcf_reader.read_header().map_err(|e| {
            FavorError::Input(format!(
                "Invalid VCF header in {}: {e}",
                input_path.display()
            ))
        })?;

        for result in vcf_reader.records() {
            let record = result.map_err(|e| {
                FavorError::Analysis(format!("VCF parse error in {}: {e}", input_path.display()))
            })?;

            process_record(
                &record,
                &header,
                &mut writers,
                vs_writer,
                &schema,
                &props,
                batch_size,
                &mut variant_count,
                &mut filtered_contigs,
                &mut multiallelic_split,
                output,
            )?;
        }
    }

    // Flush and close all writers
    vs_writer.set_columns(schema.fields().iter().map(|f| f.name().clone()).collect());
    for (chrom, mut cw) in writers {
        if cw.batch.len() > 0 {
            let rb = cw.batch.finish()?;
            cw.writer
                .write(&rb)
                .map_err(|e| FavorError::Resource(format!("Parquet write error: {e}")))?;
        }
        cw.writer
            .close()
            .map_err(|e| FavorError::Resource(format!("Parquet close error: {e}")))?;
        let path = vs_writer.chrom_path(chrom)?;
        let size = std::fs::metadata(&path).map_or(0, |m| m.len());
        vs_writer.register_chrom(chrom, cw.count, size);
    }

    Ok(VcfIngestResult {
        variant_count,
        filtered_contigs,
        multiallelic_split,
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
) -> Result<&'a mut ChromWriter, FavorError> {
    if !writers.contains_key(chrom) {
        output.status(&format!("  chr{chrom}..."));
        let path = vs_writer.chrom_path(chrom)?;
        let f = File::create(&path).map_err(|e| {
            FavorError::Resource(format!("Cannot create '{}': {e}", path.display()))
        })?;
        let w = ArrowWriter::try_new(f, schema.clone(), Some(props.clone()))
            .map_err(|e| FavorError::Resource(format!("Parquet writer init: {e}")))?;
        writers.insert(
            chrom,
            ChromWriter {
                batch: BatchBuilder::new(batch_size),
                writer: w,
                count: 0,
            },
        );
    }
    // Just inserted above if missing, so key is guaranteed present
    Ok(writers.get_mut(chrom).unwrap())
}

/// Process a single VCF record — split multi-allelics, normalize, route to per-chrom writer.
#[allow(clippy::too_many_arguments)]
fn process_record(
    record: &noodles_vcf::Record,
    _header: &noodles_vcf::Header,
    writers: &mut HashMap<&'static str, ChromWriter>,
    vs_writer: &VariantSetWriter,
    schema: &Arc<Schema>,
    props: &WriterProperties,
    batch_size: usize,
    variant_count: &mut u64,
    filtered_contigs: &mut u64,
    multiallelic_split: &mut u64,
    output: &dyn Output,
) -> Result<(), FavorError> {
    let raw_chrom = record.reference_sequence_name();
    let chrom = match normalize_chrom(raw_chrom) {
        Some(c) => c,
        None => {
            *filtered_contigs += 1;
            return Ok(());
        }
    };

    let pos = match record.variant_start() {
        Some(Ok(p)) => p.get() as i32,
        _ => return Ok(()),
    };

    let ref_allele = record.reference_bases().to_uppercase();

    let alt_bases = record.alternate_bases();
    let alt_str: &str = alt_bases.as_ref();
    let alts: Vec<&str> = alt_str.split(',').collect();

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
    let filter_str = if filter_raw.is_empty() || filter_raw == "." {
        None
    } else {
        Some(filter_raw.to_string())
    };

    if alts.len() > 1 {
        *multiallelic_split += alts.len() as u64 - 1;
    }

    for alt in &alts {
        let alt_upper = alt.trim().to_uppercase();
        // Skip missing, spanning deletion, and symbolic alleles (<DEL>, <INS>, etc.)
        if alt_upper == "*"
            || alt_upper == "."
            || alt_upper.is_empty()
            || alt_upper.starts_with('<')
        {
            continue;
        }

        let (norm_ref, norm_alt, norm_pos) = parsimony_normalize(&ref_allele, &alt_upper, pos);

        let cw =
            get_or_create_writer(chrom, writers, vs_writer, schema, props, batch_size, output)?;

        cw.batch.push(
            chrom,
            norm_pos,
            &norm_ref,
            &norm_alt,
            rsid,
            qual_str.as_deref(),
            filter_str.as_deref(),
        );
        cw.count += 1;
        *variant_count += 1;

        if cw.batch.is_full() {
            let rb = cw.batch.finish()?;
            cw.writer
                .write(&rb)
                .map_err(|e| FavorError::Resource(format!("Parquet write error: {e}")))?;

            if *variant_count % 1_000_000 == 0 {
                output.status(&format!("  {} variants processed...", variant_count));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::parsimony_normalize;

    #[test]
    fn snv() {
        assert_eq!(
            parsimony_normalize("G", "A", 100),
            ("G".into(), "A".into(), 100)
        );
    }

    #[test]
    fn insertion_with_anchor() {
        // VCF: REF=A, ALT=ACGT — insertion of CGT after pos 100
        assert_eq!(
            parsimony_normalize("A", "ACGT", 100),
            ("A".into(), "ACGT".into(), 100)
        );
    }

    #[test]
    fn deletion_with_anchor() {
        // VCF: REF=ACGT, ALT=A — deletion of CGT at pos 101-103
        assert_eq!(
            parsimony_normalize("ACGT", "A", 100),
            ("ACGT".into(), "A".into(), 100)
        );
    }

    #[test]
    fn redundant_padding_right() {
        // REF=ACGTG, ALT=AG — common suffix G, then common prefix A
        // After right-trim: ACGT, A; left-trim: CGT→""? No, keep 1 base.
        // Right-trim G: ACGT, A (suffix=1, but A is len 1 so max_suffix=0)
        // Actually: rlen=5, alen=2, max_suffix=min(5,2)-1=1. G==G → suffix_trim=1.
        // After: ACGT (4 bytes), A (1 byte). max_prefix=min(4,1)-1=0. No prefix trim.
        assert_eq!(
            parsimony_normalize("ACGTG", "AG", 100),
            ("ACGT".into(), "A".into(), 100)
        );
    }

    #[test]
    fn redundant_prefix_and_suffix() {
        // REF=AACG, ALT=AATG — common prefix AA, common suffix G
        // Right-trim: AAC, AAT (suffix=1). Left-trim: C, T (prefix=2, pos+2)
        assert_eq!(
            parsimony_normalize("AACG", "AATG", 100),
            ("C".into(), "T".into(), 102)
        );
    }

    #[test]
    fn complex_mnv_no_trim() {
        assert_eq!(
            parsimony_normalize("ACG", "TCA", 100),
            ("ACG".into(), "TCA".into(), 100)
        );
    }

    #[test]
    fn identical_alleles() {
        assert_eq!(
            parsimony_normalize("A", "A", 100),
            ("A".into(), "A".into(), 100)
        );
    }

    #[test]
    fn single_base_insertion() {
        // REF=T, ALT=TC
        assert_eq!(
            parsimony_normalize("T", "TC", 100),
            ("T".into(), "TC".into(), 100)
        );
    }
}
