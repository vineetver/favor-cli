//! Sample → variant_vcf reverse index.
//!
//! Today's sparse_g.bin is variant-major: O(n_variants) to find every
//! variant a sample carries. This lookup builds a sample-major index
//! per chromosome by walking each variant's carrier list once and
//! emitting (sample_idx, variant_vcf) pairs sorted by sample_idx.
//!
//! On-disk format per chromosome — `chromosome={chrom}.bin`:
//!
//! ```text
//! [0..4]            magic = "SCAR"
//! [4..6]            version: u16 = 1
//! [6..10]           n_samples: u32
//! [10..18]          n_pairs: u64
//! [18..18+4*(n+1)]  offsets[n+1] of u32   — `entries[offsets[s]..offsets[s+1]]`
//! [...]             entries: u32 variant_vcfs sorted within each sample
//! ```
//!
//! `READY` is written last as the build sentinel the framework checks.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use crate::error::CohortError;
use crate::store::cohort::types::{SampleIdx, VariantVcf};
use crate::store::cohort::CohortHandle;
use crate::store::manifest::write_atomic;

use super::Lookup;

const MAGIC: &[u8; 4] = b"SCAR";
const VERSION: u16 = 1;
const HEADER_SIZE: usize = 18;

pub struct SampleCarriersLookup;

impl Lookup for SampleCarriersLookup {
    const NAME: &'static str = "sample_carriers";

    type Key = (String, SampleIdx);
    type Value = Vec<VariantVcf>;

    fn build(&self, cohort: &CohortHandle<'_>, index_dir: &Path) -> Result<(), CohortError> {
        let chroms = cohort.chromosomes()?.to_vec();
        for chrom in &chroms {
            let view = cohort.chromosome(chrom)?;
            let n_samples = view.n_samples()?;
            let n_variants = view.n_variants()?;
            let sg = view.sparse_g()?;

            // Two-pass build: count per sample, then scatter.
            let mut counts = vec![0u32; n_samples as usize];
            for v in 0..n_variants {
                let cl = sg.load_variant(v);
                for entry in &cl.entries {
                    counts[entry.sample_idx as usize] += 1;
                }
            }
            let mut offsets = vec![0u32; n_samples as usize + 1];
            let mut total = 0u32;
            for (i, &c) in counts.iter().enumerate() {
                offsets[i] = total;
                total += c;
            }
            offsets[n_samples as usize] = total;

            let mut cursor = offsets.clone();
            cursor.pop();
            let mut entries = vec![0u32; total as usize];
            for v in 0..n_variants {
                let cl = sg.load_variant(v);
                for e in &cl.entries {
                    let s = e.sample_idx as usize;
                    let pos = cursor[s] as usize;
                    entries[pos] = v;
                    cursor[s] += 1;
                }
            }

            // The per-sample lists come out variant-ordered because the
            // outer loop walks variants ascending.

            let path = index_dir.join(format!("chromosome={}.bin", chrom.label()));
            let f = File::create(&path).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", path.display()))
            })?;
            let mut w = BufWriter::with_capacity(1 << 20, f);
            w.write_all(MAGIC)?;
            w.write_all(&VERSION.to_le_bytes())?;
            w.write_all(&n_samples.to_le_bytes())?;
            w.write_all(&(total as u64).to_le_bytes())?;
            for off in &offsets {
                w.write_all(&off.to_le_bytes())?;
            }
            for vcf in &entries {
                w.write_all(&vcf.to_le_bytes())?;
            }
            w.flush()?;
            w.into_inner()
                .map_err(|e| CohortError::Resource(format!("flush {}: {e}", path.display())))?
                .sync_all()
                .map_err(|e| CohortError::Resource(format!("fsync {}: {e}", path.display())))?;
        }

        write_atomic(&index_dir.join(Self::ready_marker()), b"1\n")
    }

    fn query(
        &self,
        _cohort: &CohortHandle<'_>,
        index_dir: &Path,
        key: &Self::Key,
    ) -> Result<Self::Value, CohortError> {
        let (chrom, sample) = key;
        let path = index_dir.join(format!("chromosome={chrom}.bin"));
        let bytes = std::fs::read(&path)
            .map_err(|e| CohortError::Resource(format!("read {}: {e}", path.display())))?;
        if bytes.len() < HEADER_SIZE {
            return Err(CohortError::Resource(format!(
                "{}: header truncated",
                path.display()
            )));
        }
        if &bytes[0..4] != MAGIC {
            return Err(CohortError::Resource(format!(
                "{}: bad magic",
                path.display()
            )));
        }
        let version = u16::from_le_bytes(bytes[4..6].try_into().unwrap());
        if version != VERSION {
            return Err(CohortError::Resource(format!(
                "{}: version {version} != {VERSION}",
                path.display()
            )));
        }
        let n_samples = u32::from_le_bytes(bytes[6..10].try_into().unwrap());
        let _n_pairs = u64::from_le_bytes(bytes[10..18].try_into().unwrap());

        if sample.0 >= n_samples {
            return Err(CohortError::Input(format!(
                "sample idx {} out of range (n_samples={n_samples})",
                sample.0
            )));
        }

        let offsets_start = HEADER_SIZE;
        let offsets_len = (n_samples as usize + 1) * 4;
        if bytes.len() < offsets_start + offsets_len {
            return Err(CohortError::Resource(format!(
                "{}: offsets truncated",
                path.display()
            )));
        }
        let s = sample.0 as usize;
        let off_lo = u32::from_le_bytes(
            bytes[offsets_start + s * 4..offsets_start + s * 4 + 4]
                .try_into()
                .unwrap(),
        );
        let off_hi = u32::from_le_bytes(
            bytes[offsets_start + (s + 1) * 4..offsets_start + (s + 1) * 4 + 4]
                .try_into()
                .unwrap(),
        );
        let entries_start = offsets_start + offsets_len;
        let lo = entries_start + off_lo as usize * 4;
        let hi = entries_start + off_hi as usize * 4;
        if hi > bytes.len() {
            return Err(CohortError::Resource(format!(
                "{}: entries truncated",
                path.display()
            )));
        }
        let mut out = Vec::with_capacity((off_hi - off_lo) as usize);
        for chunk in bytes[lo..hi].chunks_exact(4) {
            out.push(VariantVcf(u32::from_le_bytes(chunk.try_into().unwrap())));
        }
        Ok(out)
    }
}

