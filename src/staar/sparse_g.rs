//! SparseG: memory-mapped sparse genotype matrix over (sample_id, variant_vcf).
//!
//! The matrix stores only non-reference genotypes as carrier lists per variant.
//! Variants are ordered by variant_vcf (a dense index assigned at build time,
//! monotonic by position, ref, alt). An offsets table provides O(1) random
//! access to any variant's carrier list.
//!
//! All queries (region, gene, MAF, annotation) resolve to variant_vcf masks
//! externally. SparseG only knows about (sample_id, variant_vcf) → dosage.

use std::fs::File;
use std::path::Path;

use memmap2::Mmap;

use crate::error::FavorError;
use crate::staar::carrier::encoding::*;
use crate::staar::carrier::reader::{CarrierEntry, CarrierList};

/// Memory-mapped sparse genotype matrix for one chromosome.
///
/// Indexed by variant_vcf. No gene concept — just (sample_id, variant_vcf) → dosage.
/// variant_vcf is a dense, immutable index over the variant universe for this chromosome.
#[allow(dead_code)] // accessors used by store_validate and future diagnostics
pub struct SparseG {
    mmap: Mmap,
    n_samples: u32,
    n_variants: u32,
    wide: bool,
    /// Byte offset of each variant's carrier data relative to SPARSE_G_HEADER_SIZE.
    /// offsets[variant_vcf] = byte position of that variant's carrier block.
    offsets: Vec<u64>,
}

impl SparseG {
    /// Open a sparse genotype matrix for one chromosome.
    pub fn open(chrom_dir: &Path) -> Result<Self, FavorError> {
        let path = chrom_dir.join("sparse_g.bin");
        let file = File::open(&path)
            .map_err(|e| FavorError::Resource(format!("Open {}: {e}", path.display())))?;
        let mmap = unsafe {
            Mmap::map(&file)
                .map_err(|e| FavorError::Resource(format!("mmap {}: {e}", path.display())))?
        };

        let header = SparseGHeader::read_from(&mmap)
            .map_err(|e| FavorError::Resource(format!("sparse_g.bin: {e}")))?;

        let n_variants = header.n_variants as usize;
        let offsets_start = header.offsets_start as usize;
        let offsets_bytes = n_variants * 8;

        if mmap.len() < offsets_start + offsets_bytes {
            return Err(FavorError::Resource(format!(
                "sparse_g.bin truncated: {} bytes < expected {} (offsets table)",
                mmap.len(),
                offsets_start + offsets_bytes,
            )));
        }

        // Read the offsets table
        let mut offsets = Vec::with_capacity(n_variants);
        for i in 0..n_variants {
            let base = offsets_start + i * 8;
            let val = u64::from_le_bytes(mmap[base..base + 8].try_into().unwrap());
            offsets.push(val);
        }

        Ok(Self {
            mmap,
            n_samples: header.n_samples,
            n_variants: header.n_variants,
            wide: header.wide_index(),
            offsets,
        })
    }

    #[allow(dead_code)]
    pub fn n_samples(&self) -> u32 {
        self.n_samples
    }

    #[allow(dead_code)]
    pub fn n_variants(&self) -> u32 {
        self.n_variants
    }

    /// Raw offsets table. offsets[v] = byte offset of variant v's carrier data
    /// relative to SPARSE_G_HEADER_SIZE.
    #[allow(dead_code)]
    pub fn offsets(&self) -> &[u64] {
        &self.offsets
    }

    /// Total carrier data size in bytes (from header to offsets table start).
    #[allow(dead_code)]
    pub fn carrier_data_size(&self) -> u64 {
        if let Some(&_last) = self.offsets.last() {
            // The carrier data region ends where the offsets table begins.
            // offsets_start - HEADER_SIZE is the total carrier data bytes.
            // But we can also compute it from the mmap: offsets_start is stored
            // in the header. We just need the file-level size of carrier data.
            // The simplest bound: offsets table starts at some byte, carrier data
            // precedes it.
            let offsets_start_in_file = self.mmap.len() - self.offsets.len() * 8;
            (offsets_start_in_file - SPARSE_G_HEADER_SIZE) as u64
        } else {
            0
        }
    }

    /// O(1) access: load carrier list for one variant by variant_vcf.
    #[inline]
    pub fn load_variant(&self, variant_vcf: u32) -> CarrierList {
        let byte_offset = self.offsets[variant_vcf as usize];
        let data_start = SPARSE_G_HEADER_SIZE + byte_offset as usize;
        self.parse_carrier_list(&self.mmap[data_start..])
    }

    /// Batch access with block prefetching.
    ///
    /// Sorts variant_vcfs internally for sequential mmap access (better page
    /// cache and prefetcher behavior), reads carrier data in one sweep, then
    /// permutes results back to the caller's order.
    pub fn load_variants(&self, variant_vcfs: &[u32]) -> Vec<CarrierList> {
        let n = variant_vcfs.len();
        if n == 0 {
            return Vec::new();
        }

        // Build (original_index, variant_vcf) pairs, sorted by variant_vcf
        // for sequential mmap access
        let mut indexed: Vec<(usize, u32)> = variant_vcfs
            .iter()
            .enumerate()
            .map(|(i, &v)| (i, v))
            .collect();
        indexed.sort_unstable_by_key(|&(_, v)| v);

        // Read in sorted order (sequential mmap access)
        let mut sorted_results: Vec<(usize, CarrierList)> = Vec::with_capacity(n);
        for &(orig_idx, variant_vcf) in &indexed {
            sorted_results.push((orig_idx, self.load_variant(variant_vcf)));
        }

        // Permute back to caller's order
        sorted_results.sort_unstable_by_key(|&(i, _)| i);
        sorted_results.into_iter().map(|(_, cl)| cl).collect()
    }

    /// Parse a carrier list from a byte slice starting at the carrier count.
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
