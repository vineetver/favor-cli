//! Format-agnostic variant reader.
//!
//! The reader drives iteration and calls a processor closure on each record.
//! Callback form (not a lending iterator) because noodles VCF field wrappers
//! (`AlternateBases<'r>`, `Samples<'r>`, ...) expose their inner `&'r str`
//! only through `AsRef::as_ref(&self) -> &str`, which elides to the wrapper's
//! borrow rather than the record's. Keeping iteration inside the reader lets
//! wrapper locals stay alive across one closure call without any buffer copy.
//!
//! Per-variant normalization (`normalize_chrom`, `parsimony_normalize`,
//! `ascii_uppercase_into`) and multi-allelic split stay in the processor.

use crate::error::CohortError;

/// One variant as seen by the processor. Slices are valid only during the
/// `for_each` closure invocation that produced them.
pub struct RawRecord<'a> {
    pub chromosome: &'a str,
    pub position: i32,
    pub ref_allele: &'a str,
    /// Raw ALT field. Comma-separated for multi-allelic sites; the processor
    /// splits and parsimony-normalizes per alt.
    pub alt_alleles: &'a str,
    pub rsid: Option<&'a str>,
    pub qual: Option<&'a str>,
    pub filter: Option<&'a str>,
    /// Raw per-sample text (tab-separated, FORMAT-prefixed for VCF). Passed
    /// through to `GenotypeWriter::push` unchanged.
    pub samples_text: &'a str,
}

/// Streaming variant reader. One reader per file; returning `Err` from `f`
/// terminates iteration.
pub trait VariantReader: Send {
    fn for_each(
        &mut self,
        f: &mut dyn for<'a> FnMut(RawRecord<'a>) -> Result<(), CohortError>,
    ) -> Result<(), CohortError>;

    /// Sample names in the order they appear in the per-record `samples_text`.
    /// Reading is allowed to do work (header parse, FFI call), so callers
    /// should cache the result if they need it more than once.
    fn sample_names(&mut self) -> Result<Vec<String>, CohortError>;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ingest::vcf::VcfVariantReader;
    use std::io::Write;

    /// Opens a tiny inline VCF and asserts canonical fields including a
    /// multi-allelic site and a missing-filter row. Reader-level only: no
    /// `RecordContext`, no parquet side effects.
    #[test]
    fn vcf_reader_yields_canonical_records() {
        let mut tmp = tempfile::NamedTempFile::new().unwrap();
        writeln!(tmp, "##fileformat=VCFv4.3").unwrap();
        writeln!(tmp, "##contig=<ID=22>").unwrap();
        writeln!(tmp, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">").unwrap();
        writeln!(
            tmp,
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2"
        )
        .unwrap();
        writeln!(tmp, "22\t100\trs1\tA\tG\t30.5\tPASS\t.\tGT\t0/0\t0/1").unwrap();
        writeln!(tmp, "22\t200\t.\tA\tC,T\t.\t.\t.\tGT\t0/1\t1/2").unwrap();
        writeln!(tmp, "22\t300\trs3\tC\tT\t40\tLowQual\t.\tGT\t1/1\t0/0").unwrap();
        tmp.flush().unwrap();

        let mut reader = VcfVariantReader::open(tmp.path(), 1).unwrap();

        struct Row {
            chrom: String,
            pos: i32,
            r: String,
            a: String,
            rsid: Option<String>,
            qual: Option<String>,
            filter: Option<String>,
        }
        let mut rows: Vec<Row> = Vec::new();
        let mut samples_first: Option<String> = None;
        reader
            .for_each(&mut |rec| {
                rows.push(Row {
                    chrom: rec.chromosome.to_string(),
                    pos: rec.position,
                    r: rec.ref_allele.to_string(),
                    a: rec.alt_alleles.to_string(),
                    rsid: rec.rsid.map(str::to_string),
                    qual: rec.qual.map(str::to_string),
                    filter: rec.filter.map(str::to_string),
                });
                if samples_first.is_none() {
                    samples_first = Some(rec.samples_text.to_string());
                }
                Ok(())
            })
            .unwrap();

        assert_eq!(rows.len(), 3);

        assert_eq!(rows[0].chrom, "22");
        assert_eq!(rows[0].pos, 100);
        assert_eq!(rows[0].r, "A");
        assert_eq!(rows[0].a, "G");
        assert_eq!(rows[0].rsid.as_deref(), Some("rs1"));
        assert_eq!(rows[0].qual.as_deref(), Some("30.5"));
        assert_eq!(rows[0].filter.as_deref(), Some("PASS"));

        assert_eq!(rows[1].pos, 200);
        assert_eq!(rows[1].a, "C,T");
        assert!(rows[1].rsid.is_none());
        assert!(rows[1].qual.is_none());
        assert!(rows[1].filter.is_none());

        assert_eq!(rows[2].pos, 300);
        assert_eq!(rows[2].filter.as_deref(), Some("LowQual"));

        // noodles 0.73 `Samples` covers everything after INFO, which for VCF
        // includes the FORMAT column; `GenotypeWriter::push` parses this.
        assert_eq!(samples_first.as_deref(), Some("GT\t0/0\t0/1"));
    }
}
