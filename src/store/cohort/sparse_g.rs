//! SparseG: memory-mapped sparse genotype matrix over `(sample_id, variant_vcf)`.

use std::path::Path;

use super::encoding::*;
use super::variants::{CarrierEntry, CarrierList};
use crate::error::CohortError;
use crate::store::backend::{Backend, MappedBytes};

pub struct SparseG {
    mmap: MappedBytes,
    wide: bool,
    offsets: Vec<u64>,
}

impl SparseG {
    pub fn open(backend: &dyn Backend, chrom_dir: &Path) -> Result<Self, CohortError> {
        let path = chrom_dir.join("sparse_g.bin");
        let mmap = backend.mmap(&path)?;
        Self::from_mmap(mmap)
    }

    pub fn from_mmap(mmap: MappedBytes) -> Result<Self, CohortError> {
        let header = SparseGHeader::read_from(&mmap)
            .map_err(|e| CohortError::Resource(format!("sparse_g.bin: {e}")))?;

        let n_variants = header.n_variants as usize;
        let offsets_start = header.offsets_start as usize;
        let offsets_bytes = n_variants * 8;

        if mmap.len() < offsets_start + offsets_bytes {
            return Err(CohortError::Resource(format!(
                "sparse_g.bin truncated: {} bytes < expected {} (offsets table)",
                mmap.len(),
                offsets_start + offsets_bytes,
            )));
        }

        let offsets: Vec<u64> = mmap[offsets_start..offsets_start + offsets_bytes]
            .chunks_exact(8)
            .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
            .collect();

        // Scoring loop jumps from one variant's carrier block to another,
        // so default sequential readahead fetches pages we never touch and
        // evicts ones we just used. MADV_RANDOM turns readahead off; the
        // kernel now only pages in what we actually read. Best-effort:
        // `advise` is a hint, ignored on unsupported kernels.
        #[cfg(unix)]
        mmap.advise(memmap2::Advice::Random);

        Ok(Self {
            mmap,
            wide: header.wide_index(),
            offsets,
        })
    }

    #[inline]
    pub fn load_variant(&self, variant_vcf: u32) -> CarrierList {
        let byte_offset = self.offsets[variant_vcf as usize];
        let data_start = SPARSE_G_HEADER_SIZE + byte_offset as usize;
        self.parse_carrier_list(&self.mmap[data_start..])
    }

    pub fn load_variants(&self, variant_vcfs: &[u32]) -> Vec<CarrierList> {
        let n = variant_vcfs.len();
        if n == 0 {
            return Vec::new();
        }

        let mut indexed: Vec<(usize, u32)> = variant_vcfs
            .iter()
            .enumerate()
            .map(|(i, &v)| (i, v))
            .collect();
        indexed.sort_unstable_by_key(|&(_, v)| v);

        let mut sorted_results: Vec<(usize, CarrierList)> = Vec::with_capacity(n);
        for &(orig_idx, variant_vcf) in &indexed {
            sorted_results.push((orig_idx, self.load_variant(variant_vcf)));
        }

        sorted_results.sort_unstable_by_key(|&(i, _)| i);
        sorted_results.into_iter().map(|(_, cl)| cl).collect()
    }

    #[inline]
    fn parse_carrier_list(&self, data: &[u8]) -> CarrierList {
        let n_carriers = u16::from_le_bytes([data[0], data[1]]) as usize;
        let mut pos = CARRIER_COUNT_SIZE;
        let mut entries = Vec::with_capacity(n_carriers);

        if self.wide {
            for _ in 0..n_carriers {
                let sample_id =
                    u32::from_le_bytes([data[pos], data[pos + 1], data[pos + 2], data[pos + 3]]);
                let dosage = data[pos + 4];
                entries.push(CarrierEntry {
                    sample_idx: sample_id,
                    dosage,
                });
                pos += CARRIER_ENTRY_WIDE;
            }
        } else {
            for _ in 0..n_carriers {
                let sample_id = u16::from_le_bytes([data[pos], data[pos + 1]]) as u32;
                let dosage = data[pos + 2];
                entries.push(CarrierEntry {
                    sample_idx: sample_id,
                    dosage,
                });
                pos += CARRIER_ENTRY_NARROW;
            }
        }

        CarrierList { entries }
    }
}
