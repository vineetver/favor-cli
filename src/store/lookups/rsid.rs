//! rsID → variant_vcf lookup. Worked example for the `Lookup` trait.
//!
//! Today's cohort builder does not emit an `rsid` column on
//! `variants.parquet` (FAVOR's annotation join carries `vid` only). The
//! impl below is structurally complete: it walks the variants.parquet
//! schema, errors out cleanly when no rsid column is present, and
//! demonstrates the file format the future builder will populate.
//!
//! Index format — `rsid.bin`:
//!
//! ```text
//! [0..4]   magic = "RSID"
//! [4..6]   version: u16 = 1
//! [6..10]  n_entries: u32
//! [10..]   entries: repeated (u16 chrom_len, chrom_bytes,
//!                              u16 rsid_len, rsid_bytes,
//!                              u32 variant_vcf), sorted by rsid
//! ```
//!
//! Query is binary search over the sorted block. The format is
//! deliberately simple — the next caller can swap in an FST without
//! changing the trait surface.

use std::path::Path;

use crate::error::CohortError;
use crate::store::cohort::types::VariantVcf;
use crate::store::cohort::CohortHandle;
use crate::types::Chromosome;

use super::Lookup;

pub struct RsidLookup;

impl Lookup for RsidLookup {
    const NAME: &'static str = "rsid";
    type Key = str;
    type Value = Vec<(Chromosome, VariantVcf)>;

    fn build(&self, cohort: &CohortHandle<'_>, index_dir: &Path) -> Result<(), CohortError> {
        // Walk every chromosome's variants.parquet schema. The first
        // chromosome whose schema does not include `rsid` makes this
        // lookup unbuildable on this cohort — surface the typed error
        // so the framework's `READY` check fails fast on retry.
        let chroms = cohort.chromosomes()?.to_vec();
        for chrom in &chroms {
            let view = cohort.chromosome(chrom)?;
            let pq = view.dir().join("variants.parquet");
            let cols = crate::store::list::parquet_column_names(&pq)?;
            if !cols.iter().any(|c| c == "rsid") {
                return Err(CohortError::DataMissing(format!(
                    "rsid lookup unsupported: cohort {} has no `rsid` column on {}; \
                     re-ingest with annotations carrying rsid before retrying",
                    cohort.id().as_str(),
                    pq.display(),
                )));
            }
        }
        // The variants.parquet readers in this codebase don't yet
        // surface an `rsid` field on `VariantIndexEntry`, so even with
        // the column present we can't populate the index without
        // touching the loader. Phase 7+ adds that read; today this
        // path is a placeholder so the trait shape is exercised.
        let _ = index_dir;
        Err(CohortError::DataMissing(
            "rsid lookup builder is reserved for the future variants.parquet `rsid` reader".into(),
        ))
    }

    fn query(
        &self,
        _cohort: &CohortHandle<'_>,
        index_dir: &Path,
        key: &Self::Key,
    ) -> Result<Self::Value, CohortError> {
        let path = index_dir.join("rsid.bin");
        let bytes = std::fs::read(&path)
            .map_err(|e| CohortError::Resource(format!("read {}: {e}", path.display())))?;
        if bytes.len() < 10 || &bytes[0..4] != b"RSID" {
            return Err(CohortError::Resource(format!(
                "{}: bad header",
                path.display()
            )));
        }
        let n_entries = u32::from_le_bytes(bytes[6..10].try_into().unwrap()) as usize;

        // Linear scan today because the binary file is already sorted
        // by rsid, and the only callers are tests + the README example.
        // Swapping in a binary search is a one-line change once a real
        // workload appears.
        let mut pos = 10usize;
        let mut out = Vec::new();
        for _ in 0..n_entries {
            if pos + 2 > bytes.len() {
                break;
            }
            let chrom_len = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            let chrom_str = std::str::from_utf8(&bytes[pos..pos + chrom_len]).map_err(|e| {
                CohortError::Resource(format!("{}: bad utf8: {e}", path.display()))
            })?;
            pos += chrom_len;
            let rsid_len = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            let rsid_str = std::str::from_utf8(&bytes[pos..pos + rsid_len]).map_err(|e| {
                CohortError::Resource(format!("{}: bad utf8: {e}", path.display()))
            })?;
            pos += rsid_len;
            let vcf = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap());
            pos += 4;
            if rsid_str == key {
                let chrom: Chromosome = chrom_str
                    .parse()
                    .map_err(|e: String| CohortError::Resource(e))?;
                out.push((chrom, VariantVcf(vcf)));
            }
        }
        Ok(out)
    }
}

