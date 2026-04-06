//! Per-phenotype score cache: pre-computed U vectors and K matrices.
//!
//! Layer 2 of the three-layer STAAR architecture:
//!   Layer 1 (build once): carrier store + variant index + annotations
//!   Layer 2 (per-phenotype): null model → U, K per gene → disk cache
//!   Layer 3 (per-mask): select IDs → slice cached U, K → test (seconds)
//!
//! Disk format per chromosome (`scores.bin`, v4 — vid-keyed):
//!   Header (64 bytes): magic, version, n_variants, n_genes, sigma2
//!   Section 1: vid-keyed U vector
//!     For each variant: u16 vid_len + vid_bytes + f64 u_value
//!   Section 2: per-gene blocks:
//!     gene_name[32] + m(u32) + has_k(u8)
//!     + For each variant: u16 vid_len + vid_bytes
//!     + K[f64 × m²] (if has_k)
//!
//! At load time, vid strings are resolved to current VariantIndex positions
//! via VariantIndex::resolve_vid(). Runtime scoring uses positional Vec<u32>
//! offsets — fast array indexing, no string hashing on the hot path.
//!
//! Atomic writes: scores.bin.tmp → fsync → rename to scores.bin.
//! A killed process never leaves a partial scores.bin visible to probe().

use std::collections::HashMap;
use std::io::{BufReader, BufWriter, Read as IoRead, Write as IoWrite};
use std::path::{Path, PathBuf};

use faer::Mat;
use sha2::{Digest, Sha256};

use crate::error::FavorError;
use crate::output::Output;
use crate::staar::carrier::sparse_score;
use crate::staar::carrier::{AnalysisVectors, VariantIndex};
use crate::staar::sparse_g::SparseG;
use crate::staar::store::StoreManifest;

const MAGIC: &[u8; 8] = b"FVSCORE2";
const VERSION: u16 = 4;
const HEADER_SIZE: usize = 64;
const GENE_NAME_LEN: usize = 32;

/// Genes with more variants than this store U only (no K).
/// Layer 3 falls back to carrier store for these genes.
const MAX_K_VARIANTS: usize = 2000;

// ═══════════════════════════════════════════════════════════════════════
// Types
// ═══════════════════════════════════════════════════════════════════════

/// K matrix for one gene, with explicit mapping to chromosome-wide indices.
pub struct GeneKBlock {
    /// Maps local index 0..m → chromosome-wide VariantIndex position.
    pub variant_offsets: Vec<u32>,
    /// K matrix [m × m] row-major flat. Empty if gene exceeded MAX_K_VARIANTS.
    pub k_flat: Vec<f64>,
}

impl GeneKBlock {
    /// Whether K is cached (vs. too large — compute on demand from carrier store).
    #[inline]
    pub fn has_k(&self) -> bool {
        !self.k_flat.is_empty()
    }

    #[inline]
    pub fn m(&self) -> usize {
        self.variant_offsets.len()
    }
}

/// All cached scores for one chromosome.
pub struct ChromScoreCache {
    /// U[i] for VariantIndex entry i. Length = n_variants.
    pub u_all: Vec<f64>,
    /// Per-gene K blocks with explicit variant linkage.
    pub gene_blocks: HashMap<String, GeneKBlock>,
}

// ═══════════════════════════════════════════════════════════════════════
// Cache key
// ═══════════════════════════════════════════════════════════════════════

/// SHA-256 of (store_key, trait_name, sorted covariates, known_loci path).
/// Invalidated by: new VCF, new annotations, different trait, different covariates,
/// different --known-loci file.
/// NOT invalidated by: MAF cutoff, mask definitions, SPA flag.
pub fn cache_key(
    store_key: &str,
    trait_name: &str,
    covariates: &[String],
    known_loci: Option<&Path>,
) -> String {
    let mut hasher = Sha256::new();
    hasher.update(store_key.as_bytes());
    hasher.update(b"|");
    hasher.update(trait_name.as_bytes());
    hasher.update(b"|");
    let mut sorted = covariates.to_vec();
    sorted.sort();
    for cov in &sorted {
        hasher.update(cov.as_bytes());
        hasher.update(b",");
    }
    hasher.update(b"|");
    if let Some(p) = known_loci {
        hasher.update(p.to_string_lossy().as_bytes());
    }
    format!("{:x}", hasher.finalize())
}

/// Directory for score cache files: `{store_dir}/score_cache/{key}/`
pub fn cache_dir(store_dir: &Path, key: &str) -> PathBuf {
    store_dir.join("score_cache").join(key)
}

// ═══════════════════════════════════════════════════════════════════════
// Probe — validate existing cache
// ═══════════════════════════════════════════════════════════════════════

/// Check if a valid score cache exists for all chromosomes in the manifest.
/// Validates file existence AND header consistency (not just existence).
pub fn probe(store_dir: &Path, manifest: &StoreManifest, key: &str) -> bool {
    let dir = cache_dir(store_dir, key);
    for ci in &manifest.chromosomes {
        let path = dir.join(format!("chromosome={}/scores.bin", ci.name));
        match validate_header(&path) {
            Ok(_) => {}
            Err(_) => return false,
        }
    }
    true
}

/// Read and validate the header of a scores.bin file.
/// Returns (n_variants, n_genes) on success.
fn validate_header(path: &Path) -> Result<(u32, u32), FavorError> {
    let file = std::fs::File::open(path)
        .map_err(|e| FavorError::Resource(format!("Open {}: {e}", path.display())))?;
    let file_len = file
        .metadata()
        .map_err(|e| FavorError::Resource(format!("Stat {}: {e}", path.display())))?
        .len();
    if file_len < HEADER_SIZE as u64 {
        return Err(FavorError::Resource(format!(
            "scores.bin truncated: {} bytes < header size {HEADER_SIZE}",
            file_len
        )));
    }

    let mut reader = BufReader::new(file);
    let mut header = [0u8; HEADER_SIZE];
    reader
        .read_exact(&mut header)
        .map_err(|e| FavorError::Resource(format!("Read header: {e}")))?;

    if &header[0..8] != MAGIC {
        return Err(FavorError::Resource("scores.bin: bad magic".into()));
    }
    let version = u16::from_le_bytes([header[8], header[9]]);
    if version != VERSION {
        return Err(FavorError::Resource(format!(
            "scores.bin: version {version} != {VERSION}"
        )));
    }
    let n_variants = u32::from_le_bytes(header[10..14].try_into().unwrap());
    let n_genes = u32::from_le_bytes(header[14..18].try_into().unwrap());

    // v4: U section is variable-length (vid-keyed), so we can only check
    // that the file is larger than just the header. Full validation at load time.
    if file_len <= HEADER_SIZE as u64 && n_variants > 0 {
        return Err(FavorError::Resource(format!(
            "scores.bin truncated: {} bytes with {n_variants} variants",
            file_len
        )));
    }

    Ok((n_variants, n_genes))
}

// ═══════════════════════════════════════════════════════════════════════
// Build — compute U/K for all genes, write to disk
// ═══════════════════════════════════════════════════════════════════════

/// Intermediate result from per-gene computation.
struct GeneComputed {
    gene_name: String,
    variant_vids: Vec<Box<str>>,
    /// Positional offsets for U scatter (build-time only, not written to disk).
    variant_offsets: Vec<u32>,
    u_values: Vec<f64>,
    k_flat: Vec<f64>, // empty if m > MAX_K_VARIANTS
}

/// Build score cache for one chromosome.
///
/// Iterates all genes via VariantIndex (membership.parquet). Carrier data
/// loaded from SparseG with block prefetching. O(MAC × k) per gene.
pub fn build_chromosome(
    sparse_g: &SparseG,
    variant_index: &VariantIndex,
    analysis: &AnalysisVectors,
    out_dir: &Path,
    chrom: &str,
    out: &dyn Output,
) -> Result<(), FavorError> {
    let n_variants = variant_index.len();

    // Compute U/K per gene, scatter into chromosome-wide U vector
    let mut u_all = vec![0.0f64; n_variants];
    let mut gene_computed: Vec<GeneComputed> = Vec::with_capacity(variant_index.n_genes());

    let gene_names: Vec<String> = variant_index.gene_names().map(|s| s.to_string()).collect();

    for gene_name in &gene_names {
        let gene_vcfs = variant_index.gene_variant_vcfs(gene_name);
        let m = gene_vcfs.len();
        if m == 0 {
            continue;
        }

        // Load carrier data once for all variants in this gene
        let carriers = sparse_g.load_variants(gene_vcfs);
        let variant_offsets: Vec<u32> = gene_vcfs.to_vec();
        let variant_vids: Vec<Box<str>> = gene_vcfs
            .iter()
            .map(|&v| variant_index.get(v).vid.clone())
            .collect();

        let gc = if m > MAX_K_VARIANTS {
            let u_values = sparse_score::compute_u_only(&carriers, analysis);
            GeneComputed {
                gene_name: gene_name.clone(),
                variant_vids,
                variant_offsets,
                u_values,
                k_flat: Vec::new(),
            }
        } else {
            let (u_mat, k_mat) = sparse_score::score_gene_sparse(&carriers, analysis);
            let u_values: Vec<f64> = (0..m).map(|i| u_mat[(i, 0)]).collect();
            let mut k_flat = vec![0.0f64; m * m];
            for r in 0..m {
                for c in 0..m {
                    k_flat[r * m + c] = k_mat[(r, c)];
                }
            }
            GeneComputed {
                gene_name: gene_name.clone(),
                variant_vids,
                variant_offsets,
                u_values,
                k_flat,
            }
        };

        // Scatter U values into chromosome-wide vector
        for (local, &global) in gc.variant_offsets.iter().enumerate() {
            let gi = global as usize;
            if gi < n_variants {
                u_all[gi] = gc.u_values[local];
            }
        }
        gene_computed.push(gc);
    }

    gene_computed.sort_by(|a, b| a.gene_name.cmp(&b.gene_name));

    // Phase 3: write binary
    let chrom_dir = out_dir.join(format!("chromosome={chrom}"));
    std::fs::create_dir_all(&chrom_dir)
        .map_err(|e| FavorError::Resource(format!("mkdir {}: {e}", chrom_dir.display())))?;

    // Atomic write: write to .tmp, fsync, rename to final path.
    // A killed process never leaves a partial scores.bin visible to probe().
    let final_path = chrom_dir.join("scores.bin");
    let tmp_path = chrom_dir.join("scores.bin.tmp");
    let file = std::fs::File::create(&tmp_path)
        .map_err(|e| FavorError::Resource(format!("Create {}: {e}", tmp_path.display())))?;
    let mut w = BufWriter::new(file);

    // Header (64 bytes)
    let mut header = [0u8; HEADER_SIZE];
    header[0..8].copy_from_slice(MAGIC);
    header[8..10].copy_from_slice(&VERSION.to_le_bytes());
    header[10..14].copy_from_slice(&(n_variants as u32).to_le_bytes());
    header[14..18].copy_from_slice(&(gene_computed.len() as u32).to_le_bytes());
    header[18..26].copy_from_slice(&analysis.sigma2.to_le_bytes());
    w.write_all(&header)
        .map_err(|e| FavorError::Resource(format!("Write header: {e}")))?;

    // Section 1: vid-keyed U vector
    // Each entry: u16 vid_len + vid_bytes + f64 u_value
    for (i, entry) in variant_index.all_entries().iter().enumerate() {
        let vid = entry.vid.as_bytes();
        w.write_all(&(vid.len() as u16).to_le_bytes())
            .map_err(|e| FavorError::Resource(format!("Write vid len: {e}")))?;
        w.write_all(vid)
            .map_err(|e| FavorError::Resource(format!("Write vid: {e}")))?;
        w.write_all(&u_all[i].to_le_bytes())
            .map_err(|e| FavorError::Resource(format!("Write U: {e}")))?;
    }

    // Section 2: per-gene blocks with vid keys
    for gc in &gene_computed {
        // Gene name (null-padded to GENE_NAME_LEN)
        let mut name_buf = [0u8; GENE_NAME_LEN];
        let name_bytes = gc.gene_name.as_bytes();
        let copy_len = name_bytes.len().min(GENE_NAME_LEN - 1);
        name_buf[..copy_len].copy_from_slice(&name_bytes[..copy_len]);
        w.write_all(&name_buf)
            .map_err(|e| FavorError::Resource(format!("Write gene name: {e}")))?;

        // n_variants (u32)
        let m = gc.variant_vids.len() as u32;
        w.write_all(&m.to_le_bytes())
            .map_err(|e| FavorError::Resource(format!("Write m: {e}")))?;

        // Explicit has_k flag (1 byte): 1 if K is cached, 0 if not
        let has_k: u8 = if gc.k_flat.is_empty() { 0 } else { 1 };
        w.write_all(&[has_k])
            .map_err(|e| FavorError::Resource(format!("Write has_k: {e}")))?;

        // Per-variant vid strings (length-prefixed)
        for vid in &gc.variant_vids {
            let vid_bytes = vid.as_bytes();
            w.write_all(&(vid_bytes.len() as u16).to_le_bytes())
                .map_err(|e| FavorError::Resource(format!("Write vid len: {e}")))?;
            w.write_all(vid_bytes)
                .map_err(|e| FavorError::Resource(format!("Write vid: {e}")))?;
        }

        // K matrix [f64 × m²] — only if has_k
        for &val in &gc.k_flat {
            w.write_all(&val.to_le_bytes())
                .map_err(|e| FavorError::Resource(format!("Write K: {e}")))?;
        }
    }

    w.flush()
        .map_err(|e| FavorError::Resource(format!("Flush: {e}")))?;
    // fsync before rename to ensure data is durable
    w.into_inner()
        .map_err(|e| FavorError::Resource(format!("BufWriter finish: {e}")))?
        .sync_all()
        .map_err(|e| FavorError::Resource(format!("fsync {}: {e}", tmp_path.display())))?;
    std::fs::rename(&tmp_path, &final_path).map_err(|e| {
        FavorError::Resource(format!(
            "Rename {} -> {}: {e}",
            tmp_path.display(),
            final_path.display()
        ))
    })?;

    let n_cached = gene_computed
        .iter()
        .filter(|g| !g.k_flat.is_empty())
        .count();
    let n_large = gene_computed.len() - n_cached;
    out.status(&format!(
        "  chr{chrom}: cached U for {} variants, K for {} genes ({} large genes U-only)",
        n_variants, n_cached, n_large,
    ));

    Ok(())
}

// ═══════════════════════════════════════════════════════════════════════
// Load — read cached scores from disk
// ═══════════════════════════════════════════════════════════════════════

/// Load all cached scores for one chromosome.
/// Vid strings on disk are resolved to current VariantIndex positions at load time.
pub fn load_chromosome(
    cache_dir: &Path,
    chrom: &str,
    variant_index: &VariantIndex,
) -> Result<ChromScoreCache, FavorError> {
    let path = cache_dir.join(format!("chromosome={chrom}/scores.bin"));
    let data = std::fs::read(&path)
        .map_err(|e| FavorError::Resource(format!("Read {}: {e}", path.display())))?;

    if data.len() < HEADER_SIZE {
        return Err(FavorError::Resource("scores.bin truncated".into()));
    }

    // Parse header
    if &data[0..8] != MAGIC {
        return Err(FavorError::Resource("scores.bin: bad magic".into()));
    }
    let version = u16::from_le_bytes(data[8..10].try_into().unwrap());
    if version != VERSION {
        return Err(FavorError::Resource(format!(
            "scores.bin: version {version} != {VERSION}. Delete cache to rebuild."
        )));
    }
    let n_variants_on_disk = u32::from_le_bytes(data[10..14].try_into().unwrap()) as usize;
    let n_genes = u32::from_le_bytes(data[14..18].try_into().unwrap()) as usize;

    // Section 1: vid-keyed U vector — resolve vids to current VariantIndex positions
    let mut u_all = vec![0.0f64; variant_index.len()];
    let mut pos = HEADER_SIZE;

    for _ in 0..n_variants_on_disk {
        if pos + 2 > data.len() {
            return Err(FavorError::Resource(
                "scores.bin: U section truncated".into(),
            ));
        }
        let vid_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
        pos += 2;
        if pos + vid_len + 8 > data.len() {
            return Err(FavorError::Resource("scores.bin: U entry truncated".into()));
        }
        let vid = std::str::from_utf8(&data[pos..pos + vid_len])
            .map_err(|_| FavorError::Resource("scores.bin: invalid UTF-8 in vid".into()))?;
        pos += vid_len;
        let u_val = f64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
        pos += 8;

        // Resolve vid to current variant_vcf — skip if variant no longer exists
        if let Some(vcf_id) = variant_index.resolve_vid(vid) {
            u_all[vcf_id as usize] = u_val;
        }
    }

    // Section 2: per-gene blocks — resolve vid keys to positional offsets
    let mut gene_blocks = HashMap::with_capacity(n_genes);

    for _ in 0..n_genes {
        // Gene name (GENE_NAME_LEN) + m (4 bytes) + has_k (1 byte)
        let gene_header_size = GENE_NAME_LEN + 4 + 1;
        if pos + gene_header_size > data.len() {
            return Err(FavorError::Resource(
                "scores.bin: gene block truncated".into(),
            ));
        }
        let name_raw = &data[pos..pos + GENE_NAME_LEN];
        let name_end = name_raw
            .iter()
            .position(|&b| b == 0)
            .unwrap_or(GENE_NAME_LEN);
        let gene_name = String::from_utf8_lossy(&name_raw[..name_end]).into_owned();
        pos += GENE_NAME_LEN;

        let m = u32::from_le_bytes(data[pos..pos + 4].try_into().unwrap()) as usize;
        pos += 4;

        let has_k = data[pos] != 0;
        pos += 1;

        // Read m vid strings, resolve each to a VariantIndex position
        let mut variant_offsets = Vec::with_capacity(m);
        let mut all_resolved = true;
        for _ in 0..m {
            if pos + 2 > data.len() {
                return Err(FavorError::Resource(format!(
                    "scores.bin: gene {gene_name} vid truncated"
                )));
            }
            let vid_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
            pos += 2;
            if pos + vid_len > data.len() {
                return Err(FavorError::Resource(format!(
                    "scores.bin: gene {gene_name} vid bytes truncated"
                )));
            }
            let vid = std::str::from_utf8(&data[pos..pos + vid_len])
                .map_err(|_| FavorError::Resource("scores.bin: invalid UTF-8 in vid".into()))?;
            pos += vid_len;

            match variant_index.resolve_vid(vid) {
                Some(vcf_id) => variant_offsets.push(vcf_id),
                None => {
                    all_resolved = false;
                    break;
                }
            }
        }

        // K matrix — present only if has_k flag is set
        let k_flat = if has_k {
            let k_size = m * m * 8;
            if pos + k_size > data.len() {
                return Err(FavorError::Resource(format!(
                    "scores.bin: gene {gene_name} K matrix truncated ({} bytes needed, {} available)",
                    k_size, data.len() - pos
                )));
            }
            let k: Vec<f64> = (0..m * m)
                .map(|i| {
                    let off = pos + i * 8;
                    f64::from_le_bytes(data[off..off + 8].try_into().unwrap())
                })
                .collect();
            pos += k_size;
            k
        } else {
            Vec::new()
        };

        // Only insert gene if all vids resolved successfully
        if all_resolved && variant_offsets.len() == m {
            gene_blocks.insert(
                gene_name,
                GeneKBlock {
                    variant_offsets,
                    k_flat,
                },
            );
        }
    }

    Ok(ChromScoreCache { u_all, gene_blocks })
}

// ═══════════════════════════════════════════════════════════════════════
// Window K assembly — for SCANG and sliding windows
// ═══════════════════════════════════════════════════════════════════════

/// Assemble K matrix for a window spanning potentially multiple genes.
///
/// Returns block-diagonal K: within-gene terms from cached K blocks,
/// cross-gene terms are zero. This is more correct than the current
/// approach (which silently drops all but one gene's variants).
///
/// `window_global_indices`: chromosome-wide VariantIndex positions for
/// the variants in this window.
pub fn assemble_window_k(cache: &ChromScoreCache, window_global_indices: &[usize]) -> Mat<f64> {
    let wm = window_global_indices.len();
    let mut k = Mat::zeros(wm, wm);

    // Map global index → window-local index
    let global_to_window: HashMap<usize, usize> = window_global_indices
        .iter()
        .enumerate()
        .map(|(local, &global)| (global, local))
        .collect();

    // Group window variants by gene, extract K submatrices
    for block in cache.gene_blocks.values() {
        if !block.has_k() {
            continue;
        }
        // Find which of this gene's variants are in the window
        let pairs: Vec<(usize, usize)> = block
            .variant_offsets
            .iter()
            .enumerate()
            .filter_map(|(gene_local, &global)| {
                global_to_window
                    .get(&(global as usize))
                    .map(|&win_local| (gene_local, win_local))
            })
            .collect();
        if pairs.is_empty() {
            continue;
        }

        // Copy K submatrix entries from gene block into window K
        let gm = block.m();
        for &(gi, wi) in &pairs {
            for &(gj, wj) in &pairs {
                k[(wi, wj)] = block.k_flat[gi * gm + gj];
            }
        }
    }

    k
}

/// Slice U values for a window from the chromosome-wide U vector.
pub fn slice_window_u(cache: &ChromScoreCache, global_indices: &[usize]) -> Mat<f64> {
    let n = cache.u_all.len();
    let m = global_indices.len();
    Mat::from_fn(m, 1, |i, _| {
        let gi = global_indices[i];
        debug_assert!(gi < n, "slice_window_u: index {gi} >= u_all.len() {n}");
        cache.u_all[gi]
    })
}

// ═══════════════════════════════════════════════════════════════════════
#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: write a length-prefixed vid string.
    fn write_vid(w: &mut BufWriter<std::fs::File>, vid: &str) {
        use std::io::Write;
        let b = vid.as_bytes();
        w.write_all(&(b.len() as u16).to_le_bytes()).unwrap();
        w.write_all(b).unwrap();
    }

    /// Build a synthetic scores.bin in the v4 (vid-keyed) format and verify
    /// that the binary layout is correct.
    #[test]
    fn round_trip_write_read_v4() {
        let dir = tempfile::tempdir().unwrap();
        let chrom_dir = dir.path().join("chromosome=test");
        std::fs::create_dir_all(&chrom_dir).unwrap();

        let vids = [
            "1-100-A-T",
            "1-200-G-C",
            "1-300-A-G",
            "1-400-T-C",
            "1-500-G-A",
        ];
        let u_vals: [f64; 5] = [0.1, 0.2, 0.3, 0.4, 0.5];
        let gene_a_vids = ["1-100-A-T", "1-200-G-C", "1-300-A-G"]; // first 3
        let gene_a_k: Vec<f64> = vec![1.0, 0.1, 0.2, 0.1, 1.0, 0.3, 0.2, 0.3, 1.0];

        let path = chrom_dir.join("scores.bin");
        {
            let file = std::fs::File::create(&path).unwrap();
            let mut w = BufWriter::new(file);

            // Header
            let mut header = [0u8; HEADER_SIZE];
            header[0..8].copy_from_slice(MAGIC);
            header[8..10].copy_from_slice(&VERSION.to_le_bytes());
            header[10..14].copy_from_slice(&(vids.len() as u32).to_le_bytes());
            header[14..18].copy_from_slice(&1u32.to_le_bytes()); // 1 gene
            header[18..26].copy_from_slice(&1.5f64.to_le_bytes());
            {
                use std::io::Write;
                w.write_all(&header).unwrap();
            }

            // Section 1: vid-keyed U
            for (vid, &u) in vids.iter().zip(u_vals.iter()) {
                write_vid(&mut w, vid);
                use std::io::Write;
                w.write_all(&u.to_le_bytes()).unwrap();
            }

            // Section 2: Gene A with K
            {
                use std::io::Write;
                let mut name_buf = [0u8; GENE_NAME_LEN];
                name_buf[..5].copy_from_slice(b"GENEA");
                w.write_all(&name_buf).unwrap();
                w.write_all(&3u32.to_le_bytes()).unwrap();
                w.write_all(&[1u8]).unwrap(); // has_k
            }
            for vid in &gene_a_vids {
                write_vid(&mut w, vid);
            }
            {
                use std::io::Write;
                for &val in &gene_a_k {
                    w.write_all(&val.to_le_bytes()).unwrap();
                }
                w.flush().unwrap();
            }
        }

        // Verify raw binary structure
        let data = std::fs::read(&path).unwrap();
        assert_eq!(&data[0..8], MAGIC);
        assert_eq!(u16::from_le_bytes(data[8..10].try_into().unwrap()), VERSION);
        assert_eq!(u32::from_le_bytes(data[10..14].try_into().unwrap()), 5); // 5 variants
        assert_eq!(u32::from_le_bytes(data[14..18].try_into().unwrap()), 1); // 1 gene

        // Parse Section 1: check first vid-keyed U entry
        let mut pos = HEADER_SIZE;
        let vid_len = u16::from_le_bytes(data[pos..pos + 2].try_into().unwrap()) as usize;
        pos += 2;
        let vid = std::str::from_utf8(&data[pos..pos + vid_len]).unwrap();
        assert_eq!(vid, "1-100-A-T");
        pos += vid_len;
        let u = f64::from_le_bytes(data[pos..pos + 8].try_into().unwrap());
        assert!((u - 0.1).abs() < 1e-15);
    }

    /// Verify that cache_key changes when known_loci changes.
    #[test]
    fn cache_key_invalidation() {
        let covs = vec!["age".to_string(), "sex".to_string()];

        let k1 = cache_key("store1", "BMI", &covs, None);
        let k2 = cache_key("store1", "BMI", &covs, None);
        assert_eq!(k1, k2, "same inputs should produce same key");

        // Different trait
        let k3 = cache_key("store1", "LDL", &covs, None);
        assert_ne!(k1, k3);

        // Different covariates
        let k4 = cache_key("store1", "BMI", &["age".to_string()], None);
        assert_ne!(k1, k4);

        // Covariate order doesn't matter
        let covs_rev = vec!["sex".to_string(), "age".to_string()];
        let k5 = cache_key("store1", "BMI", &covs_rev, None);
        assert_eq!(k1, k5, "covariate order should not affect key");

        // known_loci changes key
        let k6 = cache_key("store1", "BMI", &covs, Some(Path::new("/path/to/loci.tsv")));
        assert_ne!(k1, k6, "known_loci should change cache key");

        // Different known_loci path changes key
        let k7 = cache_key("store1", "BMI", &covs, Some(Path::new("/other/loci.tsv")));
        assert_ne!(k6, k7);
    }

    /// Verify that validate_header rejects truncated files.
    #[test]
    fn validate_header_rejects_truncated() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("scores.bin");

        // Empty file
        std::fs::write(&path, []).unwrap();
        assert!(validate_header(&path).is_err());

        // Header only (no U vector)
        let mut header = [0u8; HEADER_SIZE];
        header[0..8].copy_from_slice(MAGIC);
        header[8..10].copy_from_slice(&VERSION.to_le_bytes());
        header[10..14].copy_from_slice(&100u32.to_le_bytes()); // claims 100 variants
        std::fs::write(&path, header).unwrap();
        assert!(
            validate_header(&path).is_err(),
            "should reject: header claims 100 variants but file has no U vector"
        );

        // Valid small file
        header[10..14].copy_from_slice(&0u32.to_le_bytes()); // 0 variants
        std::fs::write(&path, header).unwrap();
        assert!(validate_header(&path).is_ok());

        // Wrong magic
        let mut bad = header;
        bad[0..8].copy_from_slice(b"NOTVALID");
        std::fs::write(&path, bad).unwrap();
        assert!(validate_header(&path).is_err());

        // Wrong version
        let mut bad = header;
        bad[8..10].copy_from_slice(&99u16.to_le_bytes());
        std::fs::write(&path, bad).unwrap();
        assert!(validate_header(&path).is_err());
    }

    /// Verify assemble_window_k produces correct block-diagonal K.
    #[test]
    fn assemble_window_k_multi_gene() {
        // Two genes: A (variants at global 0,1,2) and B (variants at global 5,6)
        let mut gene_blocks = HashMap::new();
        gene_blocks.insert(
            "GENEA".to_string(),
            GeneKBlock {
                variant_offsets: vec![0, 1, 2],
                k_flat: vec![1.0, 0.1, 0.2, 0.1, 2.0, 0.3, 0.2, 0.3, 3.0],
            },
        );
        gene_blocks.insert(
            "GENEB".to_string(),
            GeneKBlock {
                variant_offsets: vec![5, 6],
                k_flat: vec![4.0, 0.5, 0.5, 5.0],
            },
        );

        let cache = ChromScoreCache {
            u_all: vec![0.0; 10],
            gene_blocks,
        };

        // Window spanning variants from both genes: global indices [1, 2, 5]
        let window_globals = vec![1, 2, 5];
        let k = assemble_window_k(&cache, &window_globals);

        assert_eq!(k.nrows(), 3);
        assert_eq!(k.ncols(), 3);

        // Within gene A: k[0,0] = K_A[1,1] = 2.0, k[0,1] = K_A[1,2] = 0.3, k[1,1] = K_A[2,2] = 3.0
        assert!((k[(0, 0)] - 2.0).abs() < 1e-15);
        assert!((k[(0, 1)] - 0.3).abs() < 1e-15);
        assert!((k[(1, 0)] - 0.3).abs() < 1e-15);
        assert!((k[(1, 1)] - 3.0).abs() < 1e-15);

        // Within gene B: k[2,2] = K_B[0,0] = 4.0
        assert!((k[(2, 2)] - 4.0).abs() < 1e-15);

        // Cross-gene terms are zero
        assert!((k[(0, 2)]).abs() < 1e-15);
        assert!((k[(1, 2)]).abs() < 1e-15);
        assert!((k[(2, 0)]).abs() < 1e-15);
        assert!((k[(2, 1)]).abs() < 1e-15);
    }

    /// Verify slice_window_u extracts correct values.
    #[test]
    fn slice_window_u_correct() {
        let cache = ChromScoreCache {
            u_all: vec![0.0, 1.1, 2.2, 3.3, 4.4, 5.5],
            gene_blocks: HashMap::new(),
        };
        let u = slice_window_u(&cache, &[1, 3, 5]);
        assert_eq!(u.nrows(), 3);
        assert!((u[(0, 0)] - 1.1).abs() < 1e-15);
        assert!((u[(1, 0)] - 3.3).abs() < 1e-15);
        assert!((u[(2, 0)] - 5.5).abs() < 1e-15);
    }
}
