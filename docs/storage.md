# Genotype Storage

COHORT stores genotypes as a sparse matrix over `(sample_id, variant_vcf) -> dosage`. All queries resolve to aligned vectors over `variant_vcf`. All operations are vectorized transformations over these aligned index spaces. No alternative variant ordering or implicit joins anywhere.

> **Coordinate convention:** All genomic positions are **1-based** (VCF convention). The FAVOR annotation data, ingested variant sets, and all pipeline output use 1-based positions. Input in 0-based formats (e.g. BED) is converted at ingest. The `variant_vcf` dense index (0, 1, 2, ..., N-1) is a separate concept — an internal array index, not a genomic coordinate.

> [Back to README](../README.md)

---

## Contents

- [Core Model](#core-model)
- [The variant_vcf Invariant](#the-variant_vcf-invariant)
- [Storage Layout](#storage-layout)
- [Sparse Carrier Lists](#sparse-carrier-lists)
- [Score Cache](#score-cache)
- [Build and Load](#build-and-load)

---

## Core Model

For a chromosome with **N** variants and **S** samples:

```
variant_vcf:  [0, 1, 2, ..., N-1]     dense index, immutable once built
```

Everything is aligned to `variant_vcf`:

```
# Base metadata vectors, all length N, row i = variant_vcf i
position:      [i32; N]
ref_allele:    [Utf8; N]
alt_allele:    [Utf8; N]
vid:           [Utf8; N]               "chr-pos-ref-alt", e.g. "22-51065592-G-A"
maf:           [f64; N]
region_type:   [RegionType; N]
consequence:   [Consequence; N]
cadd_phred:    [f64; N]
revel:         [f64; N]
regulatory:    [RegulatoryFlags; N]

# Named weight vectors, first-class, composable
weights: Vec<WeightVector>             each is { name: &str, values: [f64; N] }

# Sparse genotype matrix
G: sparse (S x N)                      G[s, v] = dosage for carrier s at variant v

# Masks, compiled on demand, never pre-materialized
gene("BRCA2")      -> &[u32]          variant_vcfs in that gene
maf < 0.01         -> boolean filter
is_plof()          -> boolean filter
```

The 11 STAAR annotation weight channels:

| Channel | Source |
|---------|--------|
| `w_cadd` | CADD PHRED |
| `w_linsight` | LINSIGHT |
| `w_fathmm_xf` | FATHMM-XF |
| `w_apc_epi_active` | Epigenetics (active) |
| `w_apc_epi_repressed` | Epigenetics (repressed) |
| `w_apc_epi_transcription` | Epigenetics (transcription) |
| `w_apc_conservation` | Conservation |
| `w_apc_protein_function` | Protein function |
| `w_apc_local_nd` | Local nucleotide diversity |
| `w_apc_mutation_density` | Mutation density |
| `w_apc_tf` | Transcription factor binding |

Every analysis reduces to:

```
mask     = compile(gene, maf_cutoff, consequence_predicate)
weights  = select(annotation_channels)[mask]
carriers = G[:, mask]                            O(total_MAC), not O(S x |mask|)
result   = f(carriers, weights)                  score tests on small arrays
```

---

## The variant_vcf Invariant

`variant_vcf` is a dense 0-based index over a fixed variant universe. Assigned at store build time. Sorted by `(position, ref_allele, alt_allele)`.

1. `variant_vcf` never changes within a store lifetime. Score caches, masks, cross-dataset references all depend on this.
2. If the variant universe changes (new VCF, different annotations), the entire store is rebuilt and all caches invalidated.
3. A content-hash fingerprint in `manifest.json` (SHA-256 of VCF head/tail/size + annotation fingerprint) enforces rebuild detection.
4. `variant_vcf = 42` always means the same physical variant for the lifetime of the store.

Metadata vectors, weight vectors, genotype columns, score cache vectors, and mask index sets all share this coordinate space. No translation layer.

---

## Storage Layout

```
.genotype_store/
+-- manifest.json
+-- samples.txt
+-- chromosome={chr}/
    +-- sparse_g.bin
    +-- variants.parquet
    +-- membership.parquet
```

### manifest.json

```json
{
  "version": 3,
  "key": "a1b2c3...",
  "n_samples": 200643,
  "n_variants": 17178,
  "chromosomes": [
    { "name": "22", "n_variants": 17178 }
  ],
  "created_at": "1743879600",
  "cohort_version": "0.1.0"
}
```

`key` = SHA-256 of the VCF fingerprint (first 1MB + last 1MB + file size) combined with the annotation fingerprint. Different input data = different key = full rebuild.

### samples.txt

One sample ID per line, in VCF column order. Line number = sample index used in `sparse_g.bin`.

### sparse_g.bin

Binary sparse genotype matrix. O(1) random access per variant via offset table.

```
+-----------------------------------------------------------+
|  Header (64 bytes)                                        |
|                                                           |
|  [0..8]   magic: "COHORT\x03\0"                           |
|  [8..10]  version: u16 = 3                                |
|  [10..14] n_samples: u32                                  |
|  [14..18] n_variants: u32                                 |
|  [18..22] flags: u32  (bit 0 = wide mode if S > 65535)   |
|  [22..30] total_carriers: u64                             |
|  [30..38] offsets_start: u64                              |
|  [38..64] reserved                                        |
+-----------------------------------------------------------+
|  Carrier data (variant_vcf 0, 1, 2, ...)                  |
|                                                           |
|  variant 0: [u16 n_carriers] [entries...]                 |
|  variant 1: [u16 n_carriers] [entries...]                 |
|  ...                                                      |
|  variant N-1: [u16 n_carriers] [entries...]               |
+-----------------------------------------------------------+
|  Offset table                                             |
|                                                           |
|  [u64; N]  byte offset of each variant's carrier data     |
+-----------------------------------------------------------+
```

Carrier entry format:

| Mode | When | Size | Layout |
|------|------|------|--------|
| Narrow | S <= 65,535 | 3 bytes | `u16 sample_idx` + `u8 dosage` |
| Wide | S > 65,535 | 5 bytes | `u32 sample_idx` + `u8 dosage` |

Only non-reference carriers stored. A variant at MAF=0.001 in 200K samples stores ~400 entries instead of 200K zeros.

To read carriers for `variant_vcf = v`:
1. Read `offsets[v]` from offset table
2. Seek to `header_size + offsets[v]`
3. Read `n_carriers` then `n_carriers` entries

### variants.parquet

Row `i` = `variant_vcf i`. ZSTD compressed. 25 columns:

| Index | Column | Type |
|-------|--------|------|
| 0 | `variant_vcf` | UInt32 |
| 1 | `position` | Int32 |
| 2 | `ref_allele` | Utf8 |
| 3 | `alt_allele` | Utf8 |
| 4 | `vid` | Utf8 |
| 5 | `maf` | Float64 |
| 6 | `region_type` | Utf8 |
| 7 | `consequence` | Utf8 |
| 8 | `cadd_phred` | Float64 |
| 9 | `revel` | Float64 |
| 10 | `is_cage_promoter` | Boolean |
| 11 | `is_cage_enhancer` | Boolean |
| 12 | `is_ccre_promoter` | Boolean |
| 13 | `is_ccre_enhancer` | Boolean |
| 14-24 | 11 weight channels | Float64 |

No `gene_name` column. Gene membership is many-to-many and stored separately.

### membership.parquet

| Column | Type |
|--------|------|
| `variant_vcf` | UInt32 |
| `gene_name` | Utf8 |

One row per (variant, gene) pair. A variant in an overlapping gene region has multiple rows. Sorted by `(gene_name, variant_vcf)`.

---

## Sparse Carrier Lists

At biobank scale with rare variants:

| | Dense matrix | Sparse carriers |
|---|---|---|
| 200K samples, 20 variants, MAC=5 | 4M entries (32 MB) | ~100 entries (1.2 KB) |
| Storage ratio | 1x | ~25,000x smaller |

This changes computational complexity for every scoring operation:

| Operation | Dense | Sparse |
|-----------|-------|--------|
| G'r (score vector) | O(S * m) | O(total_MAC) |
| G'G (kernel matrix) | O(S * m^2) | O(total_MAC * compound_hets) |
| G'X (covariate projection) | O(S * m * k) | O(total_MAC * k) |

Concrete numbers for 200K samples, 20 variants at MAC=5, 6 covariates: dense does 24M multiply-adds, sparse does 3,300. The gap grows with cohort size because total_MAC stays constant while S increases.

For a typical gene the entire working set (carrier entries + metadata + weights) fits in L1 cache (~600 bytes of carrier data for 20 variants at MAC=5).

---

## Score Cache

Pre-computed score statistics (U vectors and K matrices) cached per phenotype.

```
.genotype_store/
+-- score_cache/
    +-- {cache_key}/
        +-- chromosome={chr}/
            +-- scores.bin
```

Three layers:

**Layer 1** (build once): `sparse_g.bin` + `variants.parquet` + `membership.parquet`. Hours.

**Layer 2** (per phenotype): fit null model, compute U and K for all genes, write `scores.bin`. Minutes.

**Layer 3** (per mask): select variant_vcfs, slice cached U and K, run test battery. Seconds.

Cache key = `SHA-256(store_key, trait_name, covariates, known_loci)`. Changing MAF cutoff or mask definition only repeats Layer 3. No genotype I/O.

### scores.bin (v4, vid-keyed)

```
Header (64 bytes):
  [0..8]   magic: "FVSCORE2"
  [8..10]  version: u16 = 4
  [10..14] n_variants: u32
  [14..18] n_genes: u32
  [18..26] sigma2: f64

Section 1: vid-keyed U vector
  For each variant:
    u16 vid_len, [u8] vid_bytes, f64 u_value

Section 2: per-gene K blocks
  For each gene:
    [u8; 32] gene_name (null-padded)
    u32 m (variant count)
    u8 has_k (1 = cached, 0 = compute on demand)
    For each variant: u16 vid_len, [u8] vid_bytes
    If has_k: [f64; m * m] K matrix (row-major)
```

Variants identified by vid strings on disk, resolved to current `variant_vcf` positions at load time. Runtime scoring uses positional array indexing, not string lookups.

Genes with more than 2,000 variants store U only. K computed on demand for those.

Writes are atomic: `scores.bin.tmp` then fsync then rename. A killed process never leaves a partial file visible.

---

## Build and Load

### Build

```
VCF + FAVOR annotations
  |
  +- 1. Extract genotypes (VCF parse, 8 threads)
  |     200K samples x ~18K variants -> dosage lists
  |
  +- 2. Join with annotations (DataFusion HashJoin)
  |     genotypes x FAVOR full-tier -> annotated rare variants
  |
  +- 3. Deduplicate by (position, ref, alt)
  |
  +- 4. Write sparse_g.bin
  |     Carrier lists in variant_vcf order + offset table
  |
  +- 5. Write membership.parquet
  |     (variant_vcf, gene_name) pairs
  |
  +- 6. Write variants.parquet
        All metadata aligned to variant_vcf
```

### Load VariantIndex

1. Read `variants.parquet` into `Vec<VariantIndexEntry>` indexed by `variant_vcf`
2. Read `membership.parquet` into `HashMap<gene_name, Vec<variant_vcf>>`
3. Build `vid_to_vcf: HashMap<vid_string, u32>` for score cache resolution
4. Build 11 named weight channels as `Vec<f64>`, each length N

### Load SparseG

1. Memory-map `sparse_g.bin`
2. Validate header (magic, version, flags)
3. Load offset table (N * 8 bytes)
4. Carrier access is O(1) per variant via offset lookup

Batch loading sorts requested `variant_vcf` values for sequential mmap access, then permutes results back to the caller's order.
