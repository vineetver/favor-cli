<p align="center">
  <h1 align="center">FAVOR CLI</h1>
  <p align="center">
    Raw variants in. Rare-variant results out.
    <br />
    <strong>Annotate. Enrich. Analyze. Interpret.</strong>
    <br />
    <br />
    <a href="#install">Install</a> &middot; <a href="#quick-start">Quick Start</a> &middot; <a href="#commands">Commands</a> &middot; <a href="#roadmap">Roadmap</a> &middot; <a href="#citation">Citation</a>
  </p>
</p>

<p align="center">
  <a href="https://github.com/vineetver/favor-cli/actions/workflows/ci.yml"><img src="https://github.com/vineetver/favor-cli/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://github.com/vineetver/favor-cli/releases/latest"><img src="https://img.shields.io/github/v/release/vineetver/favor-cli?color=blue" alt="Release"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue" alt="License"></a>
  <img src="https://img.shields.io/badge/rust-stable-orange" alt="Rust">
  <img src="https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey" alt="Platform">
</p>

---

> **Pre-1.0 software.** Interfaces, resource profiles, and on-disk layouts may change between releases.

## Install

```bash
curl -fsSL https://raw.githubusercontent.com/vineetver/favor-cli/master/install.sh | sh
```

## Quick Start

```bash
# 1. configure: point at a data directory + choose annotation tier
favor setup --root /data/favor --tier base

# 2. pull annotation data (~200 GB base, ~508 GB full)
favor data pull

# 3. ingest and annotate variants
favor ingest variants.vcf.gz
favor annotate variants.ingested

# 4. run STAAR rare-variant association
favor staar --genotypes cohort.vcf.gz --phenotype pheno.tsv \
  --trait-name LDL --covariates age,sex,PC1,PC2 \
  --annotations variants.annotated
```

## Commands

| Command | What it does |
|---------|-------------|
| `favor setup` | Configure data root, annotation tier, environment |
| `favor data pull` | Download annotation parquets and optional packs |
| `favor ingest` | Normalize VCF/TSV/CSV into canonical parquet variant sets |
| `favor annotate` | Join variants against FAVOR base or full annotations |
| `favor enrich` | Overlay tissue-specific eQTL, regulatory, enhancer-gene data |
| `favor staar` | STAAR rare-variant association testing |
| `favor meta-staar` | Cross-study meta-analysis from summary statistics |
| `favor schema` | Inspect annotation table columns and types |
| `favor manifest` | Show installed data and available commands |

Use `--format json` for machine-readable output. Use `--dry-run` before heavy computation.

## Data layout

FAVOR CLI uses two separate storage areas:

**Data root** (`--root` during setup) holds annotation parquets shared across projects:

```
/data/favor/
  base/chromosome=*/sorted.parquet      # base tier (~200 GB)
  full/chromosome=*/sorted.parquet      # full tier (~508 GB)
  tissue/                               # optional enrichment packs
    reference/                          #   gene index, cCRE regions (40 MB, always installed)
    rollups/                            #   gene-level summaries (49 MB, always installed)
    variant_in_region/                  #   variant-region junction (155 GB, always installed)
    variant_eqtl/                       #   GTEx eQTL (3 GB, optional)
    region_ccre_tissue_signals/         #   ENCODE regulatory (18 GB, optional)
    ...
```

**Project store** (`.cohort/` in your working directory) holds per-project cohort data:

```
my_study/
  .cohort/
    cohorts/<id>/                       # built by favor ingest or favor staar
      manifest.json
      samples.txt
      chromosome=*/
        sparse_g.bin                    # sparse genotype matrix (mmap'd)
        variants.parquet                # variant metadata + STAAR weights
        membership.parquet              # gene-variant assignments
    cache/score_cache/                  # reused across mask/MAF reruns
    annotations/refs.toml               # attached annotation databases
```

The store root is resolved as: `--store-path` flag > `FAVOR_STORE` env > walk up for `.cohort/` > `<cwd>/.cohort/`.

See [Setup guide](docs/setup.md) for detailed configuration, pack selection, HPC tips, and working directory organization.

## Resource requirements

Tested on UKB exome chr22 (~200K samples, ~400K variants, ~17K rare) with 64 GB. Full genome not yet tested. Memory scales with sample count, not variant count.

```text
samples    RAM       notes
───────    ──────    ─────────────────────────────
 10K       32 GB     comfortable
200K       64 GB     tested (UKB exome chr22)
```

Memory, threads, and temp directory are auto-detected from SLURM and cgroup. Override with:

```text
SLURM_MEM_PER_NODE     memory pool
FAVOR_KINSHIP_MEM_GB   kinship budget (default 16 GB)
TMPDIR                 scratch space
```

## Docs

- **[Setup guide](docs/setup.md)** - installation, configuration, data management, HPC best practices
- [Genotype store](docs/storage.md) - sparse genotype store for rare-variant analysis
- [Validation](docs/validation.md) - statistical accuracy vs R reference
- [Performance](docs/performance.md) - benchmarks and optimization roadmap
- [Agent reference](AGENTS.md) - machine interface for LLM agents

## Roadmap

| Milestone | Focus |
|-----------|-------|
| [v0.2.0 - STAAR hardening](https://github.com/vineetver/favor-cli/milestone/1) | GRM, score validation, multi-VCF input, performance profiling |
| [v0.3.0 - MetaSTAAR](https://github.com/vineetver/favor-cli/milestone/2) | cross-biobank meta-analysis, allele flip, conditional, effect sizes |
| [v0.4.0 - Interpret](https://github.com/vineetver/favor-cli/milestone/3) | variant interpretation, fine-mapping, colocalization, V2G, tiers |
| [v0.5.0 - memory and thread pool overhaul](https://github.com/vineetver/favor-cli/milestone/5) | one compute handle, bounded scratch, machine-visible resource control |
| [v0.6.0 - storage and query engine](https://github.com/vineetver/favor-cli/milestone/6) | store format, query paths, incremental ingest, cloud I/O, agent-friendly queries |
| [v1.0.0 - Production](https://github.com/vineetver/favor-cli/milestone/4) | orchestration, provenance, QC, full test suite |

## Citation

FAVOR CLI implements the [STAAR](https://github.com/xihaoli/STAARpipeline) framework and the [FAVOR](https://favor.genohub.org) annotation database. If you use this tool, please cite:

> Li Z\*, Li X\*, Zhou H, et al. **A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies.** *Nature Methods*, 19(12), 1599-1611 (2022). [DOI: 10.1038/s41592-022-01640-x](https://doi.org/10.1038/s41592-022-01640-x)

> Li X\*, Li Z\*, Zhou H, et al. **Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis of large whole-genome sequencing studies at scale.** *Nature Genetics*, 52(9), 969-983 (2020). [DOI: 10.1038/s41588-020-0676-4](https://doi.org/10.1038/s41588-020-0676-4)

> Zhou H, Verma V, Li X, et al. **FAVOR 2.0: A reengineered functional annotation of variants online resource for interpreting genomic variation.** *Nucleic Acids Research*, 54(D1), D1405-D1414 (2026). [DOI: 10.1093/nar/gkaf1217](https://doi.org/10.1093/nar/gkaf1217)

> Li TC, Zhou H, Verma V, et al. **FAVOR-GPT: a generative natural language interface to whole genome variant functional annotations.** *Bioinformatics Advances*, 4(1), vbae143 (2024). [DOI: 10.1093/bioadv/vbae143](https://doi.org/10.1093/bioadv/vbae143)

## License

GPL-3.0
