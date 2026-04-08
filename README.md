<p align="center">
  <h1 align="center">COHORT CLI</h1>
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
  <a href="https://github.com/vineetver/cohort-cli/actions/workflows/ci.yml"><img src="https://github.com/vineetver/cohort-cli/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://github.com/vineetver/cohort-cli/releases/latest"><img src="https://img.shields.io/github/v/release/vineetver/cohort-cli?color=blue" alt="Release"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-GPL--3.0-blue" alt="License"></a>
  <img src="https://img.shields.io/badge/rust-stable-orange" alt="Rust">
  <img src="https://img.shields.io/badge/platform-linux%20%7C%20macos%20%7C%20windows-lightgrey" alt="Platform">
</p>

---

## Install

```bash
curl -fsSL https://raw.githubusercontent.com/vineetver/cohort-cli/master/install.sh | sh
```

## Quick Start

```bash
cohort setup
cohort data pull
cohort ingest input.vcf.gz
cohort annotate input.ingested
cohort staar --genotypes cohort.vcf.gz --phenotype pheno.tsv \
  --trait-name LDL --covariates age,sex,PC1,PC2 \
  --annotations input.annotated
```

## Commands

| Command | Purpose |
|---------|---------|
| `cohort ingest` | Normalize variant inputs into a canonical variant set |
| `cohort annotate` | Join variants against FAVOR base or full annotations |
| `cohort enrich` | Join tissue-specific tables onto annotated variants |
| `cohort staar` | Run single-study rare-variant association tests |
| `cohort meta-staar` | Run cross-study meta-analysis from summary-stat outputs |
| `cohort data` | Pull, inspect, and verify annotation packs |
| `cohort schema` | Inspect installed table schemas |
| `cohort manifest` | Show installed data and command availability |
| `cohort setup` | Configure paths, tier, and environment defaults |

Machine mode: use `--format json`. For heavy runs, use `--dry-run` first.

Docs: [Storage](docs/storage.md), [Validation](docs/validation.md), [Performance](docs/performance.md), [Agent reference](AGENTS.md).

## Roadmap

| Milestone | Focus |
|-----------|-------|
| [v0.2.0 - STAAR hardening](https://github.com/vineetver/cohort-cli/milestone/1) | GRM, score validation, multi-VCF input, performance profiling |
| [v0.3.0 - MetaSTAAR](https://github.com/vineetver/cohort-cli/milestone/2) | cross-biobank meta-analysis, allele flip, conditional, effect sizes |
| [v0.4.0 - Interpret](https://github.com/vineetver/cohort-cli/milestone/3) | variant interpretation, fine-mapping, colocalization, V2G, tiers |
| [v0.5.0 - memory and thread pool overhaul](https://github.com/vineetver/cohort-cli/milestone/5) | one compute handle, bounded scratch, machine-visible resource control |
| [v0.6.0 - storage and query engine](https://github.com/vineetver/cohort-cli/milestone/6) | store format, query paths, incremental ingest, cloud I/O, agent-friendly queries |
| [v1.0.0 - Production](https://github.com/vineetver/cohort-cli/milestone/4) | orchestration, provenance, QC, full test suite |

## Citation

COHORT CLI implements the [STAAR](https://github.com/xihaoli/STAARpipeline) framework and the [FAVOR](https://favor.genohub.org) annotation database. If you use this tool, please cite:

> Li Z\*, Li X\*, Zhou H, et al. **A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies.** *Nature Methods*, 19(12), 1599-1611 (2022). [DOI: 10.1038/s41592-022-01640-x](https://doi.org/10.1038/s41592-022-01640-x)

> Li X\*, Li Z\*, Zhou H, et al. **Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis of large whole-genome sequencing studies at scale.** *Nature Genetics*, 52(9), 969-983 (2020). [DOI: 10.1038/s41588-020-0676-4](https://doi.org/10.1038/s41588-020-0676-4)

> Zhou H, Verma V, Li X, et al. **FAVOR 2.0: A reengineered functional annotation of variants online resource for interpreting genomic variation.** *Nucleic Acids Research*, 54(D1), D1405-D1414 (2026). [DOI: 10.1093/nar/gkaf1217](https://doi.org/10.1093/nar/gkaf1217)

> Li TC, Zhou H, Verma V, et al. **FAVOR-GPT: a generative natural language interface to whole genome variant functional annotations.** *Bioinformatics Advances*, 4(1), vbae143 (2024). [DOI: 10.1093/bioadv/vbae143](https://doi.org/10.1093/bioadv/vbae143)

## License

GPL-3.0
