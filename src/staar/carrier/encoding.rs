//! Binary format for sparse_g.bin — the sparse genotype matrix.
//!
//! The matrix is indexed by (sample_id, variant_vcf) → dosage.
//! Carrier data is stored per-variant in variant_vcf order.
//! An offsets table at the end enables O(1) random access by variant_vcf.

use std::io::{self, Write};

pub const SPARSE_G_MAGIC: [u8; 8] = *b"COHORT\x03\0";
pub const SPARSE_G_VERSION: u16 = 3;
pub const SPARSE_G_HEADER_SIZE: usize = 64;

/// Header flag: sample indices are u32 (n_samples > 65535).
pub const FLAG_WIDE_INDEX: u32 = 1 << 0;

/// sparse_g.bin header — 64 bytes on disk.
///
/// ```text
/// [0..8]    magic
/// [8..10]   version: u16 = 3
/// [10..14]  n_samples: u32
/// [14..18]  n_variants: u32
/// [18..22]  flags: u32
/// [22..30]  total_carriers: u64
/// [30..38]  offsets_start: u64   (byte offset of the offsets table)
/// [38..64]  reserved
/// ```
#[derive(Clone, Copy, Debug)]
pub struct SparseGHeader {
    pub n_samples: u32,
    pub n_variants: u32,
    pub flags: u32,
    pub total_carriers: u64,
    pub offsets_start: u64,
}

impl SparseGHeader {
    pub fn new(n_samples: u32, n_variants: u32, total_carriers: u64, offsets_start: u64) -> Self {
        Self {
            n_samples,
            n_variants,
            flags: if n_samples > 65535 {
                FLAG_WIDE_INDEX
            } else {
                0
            },
            total_carriers,
            offsets_start,
        }
    }

    pub fn wide_index(&self) -> bool {
        self.flags & FLAG_WIDE_INDEX != 0
    }

    pub fn write_to<W: Write>(&self, w: &mut W) -> io::Result<()> {
        w.write_all(&SPARSE_G_MAGIC)?; // 0..8
        w.write_all(&SPARSE_G_VERSION.to_le_bytes())?; // 8..10
        w.write_all(&self.n_samples.to_le_bytes())?; // 10..14
        w.write_all(&self.n_variants.to_le_bytes())?; // 14..18
        w.write_all(&self.flags.to_le_bytes())?; // 18..22
        w.write_all(&self.total_carriers.to_le_bytes())?; // 22..30
        w.write_all(&self.offsets_start.to_le_bytes())?; // 30..38
        w.write_all(&[0u8; 26])?; // 38..64
        Ok(())
    }

    pub fn read_from(bytes: &[u8]) -> Result<Self, &'static str> {
        if bytes.len() < SPARSE_G_HEADER_SIZE {
            return Err("too short for sparse_g header");
        }
        if bytes[0..8] != SPARSE_G_MAGIC {
            return Err("bad sparse_g magic");
        }
        let version = u16::from_le_bytes([bytes[8], bytes[9]]);
        if version != SPARSE_G_VERSION {
            return Err("unsupported sparse_g version");
        }
        // Header length checked above; all slice→array conversions are infallible.
        let read_u32 = |off: usize| u32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
        let read_u64 = |off: usize| u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap());
        Ok(Self {
            n_samples: read_u32(10),
            n_variants: read_u32(14),
            flags: read_u32(18),
            total_carriers: read_u64(22),
            offsets_start: read_u64(30),
        })
    }
}

/// Per-carrier entry size for narrow mode: u16 sample_id + u8 dosage.
pub const CARRIER_ENTRY_NARROW: usize = 3;

/// Per-carrier entry size for wide mode: u32 sample_id + u8 dosage.
pub const CARRIER_ENTRY_WIDE: usize = 5;

/// Carrier count prefix: u16 n_carriers.
pub const CARRIER_COUNT_SIZE: usize = 2;
