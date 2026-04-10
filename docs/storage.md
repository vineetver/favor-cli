# Genotype Store

Takes VCF genotypes and FAVOR annotations, writes them sparse and aligned
to a single variant index so STAAR/burden/SKAT runs are fast.

```text
  VCF / GDS / TileDB-VCF
         │  ingest + annotate
         ▼
    genotype store            <- this doc
         │  mask + score
         ▼
    STAAR / burden / SKAT
```

If the input VCF changes, the store is rebuilt.

> [Back to README](../README.md) · [Performance](performance.md)

## Why Sparse?

Rare variants. Most cells are zero:

```text
               variants →
samples ↓      v0   v1   v2   v3   v4   v5
            +----+----+----+----+----+----+
sample 0    | 0  | 0  | 1  | 0  | 0  | 2  |
sample 1    | 0  | 0  | 0  | 0  | 0  | 0  |
sample 2    | 1  | 0  | 0  | 0  | 1  | 0  |
sample 3    | 0  | 0  | 0  | 0  | 0  | 0  |
sample 4    | 0  | 2  | 0  | 0  | 0  | 0  |
            +----+----+----+----+----+----+

Dense storage pays for every cell.
```

We only keep non-zero entries:

```text
1) samples.txt             2) variants.parquet
--------------             -------------------
0 -> sample 0              variant_vcf  vid            pos   maf  ...
1 -> sample 1              0            22-...-A-G     ...
2 -> sample 2              1            22-...-C-T     ...
3 -> sample 3              2            22-...-G-A     ...
4 -> sample 4              3            22-...-T-C     ...
                           4            22-...-A-AT    ...
                           5            22-...-G-C     ...

3) sparse_g.bin
---------------
offsets:                   carrier blocks:
v0 -> byte 100             v0: [(sample 2, dosage 1)]
v1 -> byte 112             v1: [(sample 4, dosage 2)]
v2 -> byte 121             v2: [(sample 0, dosage 1)]
v3 -> byte 127             v3: []
v4 -> byte 130             v4: [(sample 2, dosage 1)]
v5 -> byte 139             v5: [(sample 0, dosage 2)]
```

Zeros are implied, not stored. At 200K samples, 20 variants, MAC=5:

```text
dense matrix:    4,000,000 entries
sparse carriers:       100 entries
```

Kernel costs follow:

```text
operation    dense            sparse
─────────    ────────────     ─────────────────────────
G'r          O(S × m)         O(total_MAC)
G'G          O(S × m²)        O(total_MAC × compound_hets)
G'X          O(S × m × k)     O(total_MAC × k)
```

## Core Model

One chromosome, `N` variants, `S` samples. Everything hangs off `variant_vcf`:

```text
variant_vcf: [0, 1, 2, ..., N-1]

        ┌─ position[v]
        ├─ ref_allele[v]
        ├─ alt_allele[v]
        ├─ vid[v]              "chr-pos-ref-alt"
  v ────├─ maf[v]
        ├─ region_type[v]
        ├─ consequence[v]
        ├─ weights[channel][v]
        └─ G[:, v]             sparse carrier list
```

Queries are index arithmetic:

```text
gene("BRCA2")      -> [variant_vcf]       gene membership lookup
maf < 0.01         -> mask[variant_vcf]    column filter
is_plof()          -> mask[variant_vcf]    column filter
carriers(v)        -> sparse list          direct offset jump

mask     = compile(gene, maf_cutoff, annotation predicate)
weights  = weights[mask]
carriers = G[:, mask]
result   = score(carriers, weights)
```

## Invariants

`variant_vcf` is dense, 0-based, fixed for the lifetime of a store. 1-based genomic positions are stored alongside but the index itself is 0-based.

1. Assigned at build time, sorted by `(position, ref_allele, alt_allele)`.
2. Variant universe changes = full rebuild + cache invalidation.
3. `variant_vcf = 42` means the same physical variant for the life of that build.
4. Metadata, masks, score-cache all share this coordinate space.

## On-Disk Layout

```text
.genotype_store/
├── manifest.json          schema version, content fingerprint, counts
├── samples.txt            sample ID -> index
├── score_cache/
│   └── {cache_key}/
│       └── chromosome={chr}/
│           └── scores.bin
└── chromosome={chr}/
    ├── sparse_g.bin       sparse genotype matrix
    ├── variants.parquet   per-variant metadata + weights
    └── membership.parquet gene -> variant mapping
```

### manifest.json

Content key = hash of VCF fingerprint + annotation fingerprint. If either input changes, store is stale.

### samples.txt

One sample ID per line, VCF column order. Line number = sample index in `sparse_g.bin`.

### sparse_g.bin

```text
┌──────────────────────────┐
│ header                   │  counts, flags, offsets location
├──────────────────────────┤
│ carrier block for v=0    │  variable length
│ carrier block for v=1    │
│ ...                      │
│ carrier block for v=N-1  │
├──────────────────────────┤
│ offset table             │  byte offset for every variant
└──────────────────────────┘
```

Carrier entry format:

```text
mode     when            size      layout
──────   ─────────────   ───────   ──────────────────────────
narrow   S ≤ 65,535      3 bytes   u16 sample_idx + u8 dosage
wide     S > 65,535      5 bytes   u32 sample_idx + u8 dosage
```

Reading one variant:

```text
offsets[v] ──jump──> carrier block ──read──> n_carriers entries
     O(1) lookup         sequential read
```

### variants.parquet

One row per `variant_vcf`:

- identity: `variant_vcf`, `position`, `ref_allele`, `alt_allele`, `vid`
- annotations: `maf`, `region_type`, `consequence`, `cadd_phred`, `revel`
- regulatory flags: CAGE and cCRE booleans
- 11 STAAR weight channels (still `Float64`, probably too wide)

### membership.parquet

Gene-to-variant mapping:

```text
variant_vcf (UInt32)  │  gene_name (Utf8)
──────────────────────┼──────────────────
 12                   │  BRCA2
 13                   │  BRCA2
 14                   │  BRCA2
 12                   │  ZNF331
```

Sorted by `(gene_name, variant_vcf)` so gene lookups are a range scan, no need
to put `gene_name` in `variants.parquet`.

## Query Patterns

### Fast

```text
query                               why
────────────────────────────────     ─────────────────────────────────────
variant by variant_vcf              O(1) offset lookup
all variants in one gene            membership.parquet pre-sorted
gene + MAF + consequence mask       metadata aligned to variant_vcf
burden / SKAT inputs for one gene   sparse carriers + aligned weights
rescore with a different mask       score cache slices by variant_vcf
```

### Workable

```text
query                               limitation
────────────────────────────────     ────────────────────────────────────────
variant by vid                      in-memory vid -> variant_vcf map at load
region query (SNVs)                 binary search, no interval index
multi-region query                  repeated range resolution, not batched
PRS-style scan over selected vars   scan + sum, no specialized primitive
```

### Gaps

```text
query                               gap
────────────────────────────────     ────────────────────────────────────────
sample -> carried variants           no sample-side index
trio / family lookups                no sample-side index
long indel / SV overlap query        no end_position, no interval sidecar
add samples or variants              full rebuild
VCF/BCF round-trip export            raw VCF structure not retained
object-store-backed reads            assumes local mmap
```

## Score Cache

Caches per-phenotype U and K so only the mask step (seconds) repeats,
not the null model fit (minutes) or the store build (hours):

```text
layer   what                             cost
─────   ──────────────────────────────   ────────
  1     build store from VCF + FAVOR     hours
  2     fit null model, compute U and K  minutes
  3     apply masks, run test battery    seconds
```

`scores.bin` layout:

```text
┌─────────────────────────────┐
│ header                      │  magic, version, counts, sigma²
├─────────────────────────────┤
│ section 1: U values         │  keyed by vid
├─────────────────────────────┤
│ section 2: K blocks         │  per gene
└─────────────────────────────┘
```

Still a custom binary format. Tracked separately.

## Build and Load

### Build

```text
VCF + FAVOR annotations
  │
  ├─ extract genotypes ──────────> sparse_g.bin
  ├─ join with annotations ──────> variants.parquet
  ├─ deduplicate by (pos,ref,alt)
  └─ extract gene membership ───> membership.parquet
```

### Load

```text
VariantIndex:                       SparseG:
  1. read variants.parquet            1. mmap sparse_g.bin
  2. read membership.parquet          2. validate header
  3. build vid -> variant_vcf         3. load offset table
  4. materialize weight vectors       4. carrier lookups by offset jump
```

Batch carrier loads sort requested `variant_vcf` values so reads stay sequential.
