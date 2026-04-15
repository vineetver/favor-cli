# Ingest

Streams VCF (plain or BGZF) into columnar parquet. For multi-sample VCFs
(biobank dosages) it also writes the sparse genotype store the downstream
STAAR path needs.

```text
  *.vcf / *.vcf.gz
         │
         │  streaming parse (noodles-vcf) + memchr tab split
         ▼
  +----------------------+       variants.parquet    (columnar, per-chr)
  |  parallel ingest     │──┬──► membership.parquet  (gene ↔ variant_vcf)
  |  (rayon)             │  │
  +----------------------+  └──► sparse_g.bin        (carrier list, mmap-ready)
```

Genotypes stay dense in memory only long enough to identify non-zero
carriers. See [Genotype store](storage.md) for the on-disk shape.

> [Back to README](../README.md) · [Genotype store](storage.md) · [STAAR](staar.md)

## Data Shapes

```text
INPUT                              TRANSFORMATION              OUTPUT
──────────────────────────         ──────────────────────      ──────────────────────────
path(s) to *.vcf[.gz]         ─┐                               variants.parquet  (per chr)
  samples in columns           ├─ parse + normalize            membership.parquet
  variants in rows             │  split multi-allelic          sparse_g.bin      (if genotypes)
                               │  parsimony (ref, alt)          samples.txt       (sidecar)
annotated variant set         ─┤  left-join on (chr, pos,      genotypes.json    (sidecar)
  (from favor annotate)        │  ref, alt)
                               └─ carrier-only genotype cast
```

## Running It

### Workstation (single box, many cores)

One command, all files. Let rayon saturate cores:

```bash
favor ingest --annotations variants.annotated/ \
  data/ukb23157_c1_b1.vcf.gz data/ukb23157_c1_b2.vcf.gz ...
```

Threads and memory are auto-detected from cgroup + SLURM; override with
`--threads` / `SLURM_MEM_PER_NODE` if needed. Files are chunked across
workers; each worker processes its chunk sequentially with one
`RecordContext`.

### HPC (per-chromosome SLURM jobs)

One job per chromosome. Each job uses the full node's cores for its
chromosome's files. Good when the cluster is bursty and a single
multi-chromosome job would starve:

```bash
# submit.sh
for chr in {1..22} X Y; do
  srun -p xlin -t 0-04:00 --mem=64G -c 16 -J "ingest_chr${chr}" \
    favor ingest --annotations variants.annotated/ \
      data/ukb23157_c${chr}_b*.vcf.gz &
done
wait
```

The parallel ingest path does not split a single file across workers —
block-level parallelism is the next milestone. If your files span one
chromosome each, the per-file parallel ingest already does the right
thing on a single node.

### File ordering

The chunker interleaves files round-robin across workers
(`files[w], files[w + N], files[w + 2N], …`). For multi-chromosome inputs
this spreads chromosomes across workers — fine for even load, but each
worker pays a ~100 ms writer-switch cost per chromosome change. If all
files for a chromosome live in a contiguous span of `--input` args, the
round-robin still works; if mixed, sort by chromosome first.

## Preflight

Before any worker starts, `ingest_vcfs_parallel` runs four cheap checks.
Each fails fast with a typed error.

```text
check                       why                                            behavior
──────────────────────      ──────────────────────────────────────────     ──────────────────
duplicate inputs            canonicalize paths; passing the same file      errors (exit 1)
                            twice double-counts carriers
stale part files            a previous run that died mid-way leaves        removes them
                            part_*.parquet under chromosome=*/; picking
                            both up corrupts variant counts
fd soft limit               each worker opens ~24 chrom writers + bgzf     caps --threads
                            reader + geno writer; default ulimit is 1024
batch-size viability        row-group must clear 500 variants or the       caps --threads, or
                            writer flushes per variant                      errors (exit 3)
```

Preflight messages are prefixed with `  Preflight:` in human mode. In
machine mode (`--format json`) they surface as structured status lines.

## Throughput Tips

```text
bottleneck                   fix                                            where it lives
──────────────────────       ──────────────────────────────────────         ──────────────────
bgzf decompression           more BGZF worker threads (auto-scales          noodles-bgzf
                             with `--threads`)                               reader::Builder
VCF line scan                SIMD tab search via memchr AVX2/SSE4.2         memchr_iter in
                             (already enabled; no tuning)                    GenotypeWriter
sample-column parse          SIMD tab find, then gt_prefix_len on
                             each GT block
row-group churn              raise per-worker memory so raw_batch_size      MIN_VIABLE_BATCH
                             clears 500; preflight will cap workers          in preflight
                             if too tight
fd exhaustion                `ulimit -n 4096` on the shell before ingest    fd soft limit
                             (preflight will cap otherwise)
```

For a 200K-sample exome chr22 on 16 cores and 64 GB the warm path sits
in the 10–15 minute range. The hot cost is BGZF decompression + sample
column parse — both SIMD already.

## Machine Mode

`--format json` emits one line per status event plus a final JSON
summary. Orchestrators consume this without parsing human text:

```bash
favor ingest --format json ... | jq 'select(.event == "variant_count")'
```

`--dry-run` validates inputs (preflight runs) and prints the derived
plan without opening writers.

## Related

- [Genotype store](storage.md) — on-disk shape of `sparse_g.bin`, `variants.parquet`, `membership.parquet`
- [STAAR](staar.md) — the consumer of what ingest produces
- [Performance](performance.md) — memory pool, hot paths
