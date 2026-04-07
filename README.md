<p align="center">
  <h1 align="center">COHORT CLI</h1>
  <p align="center">
    From raw variants to biological mechanisms in one tool.
    <br />
    <strong>Annotate. Enrich. Analyze. Interpret.</strong>
    <br />
    <br />
    <a href="#install">Install</a> &middot; <a href="#quick-start">Quick Start</a> &middot; <a href="docs/storage.md">Storage Format</a> &middot; <a href="docs/validation.md">Validation</a> &middot; <a href="AGENTS.md">Agent Reference</a> &middot; <a href="#roadmap">Roadmap</a>
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

## Quick start

```bash
cohort setup
cohort ingest input.vcf.gz
cohort annotate input.ingested.parquet
cohort enrich input.annotated.parquet --tissue brain
cohort staar --genotypes cohort.vcf.gz --phenotype pheno.tsv \
  --trait-name LDL --covariates age,sex,PC1,PC2 --annotations annotated.parquet
```

## Commands

| Command | Description |
|---------|-------------|
| `cohort setup` | Configure tier, data paths, HPC environment, memory budget |
| `cohort ingest` | Ingest any variant format: WGS VCF, variant lists, credible sets, TSV/CSV/parquet |
| `cohort annotate` | Annotate with FAVOR functional annotations (200–508 GB, 24 chromosomes) |
| `cohort enrich` | Tissue-specific overlay: eQTL, sQTL, ChromBPNet, enhancer-gene links |
| `cohort staar` | STAAR rare-variant association testing (single-study) |
| `cohort meta-staar` | MetaSTAAR cross-biobank rare-variant meta-analysis (summary-stat based) |
| `cohort schema` | Inspect annotation table schemas |

All commands support `--format json` and `--dry-run`. See [AGENTS.md](AGENTS.md) for the machine interface.

## STAAR architecture

Genotypes are stored as a canonical sparse matrix **G** over `(sample_id, variant_vcf)`. All queries — region, gene, MAF, annotation — resolve to aligned vectors over `variant_vcf`. Scoring is carrier-indexed at **O(total_MAC)**, not O(n_samples x n_variants).

| Layer | What | Cache key |
|-------|------|-----------|
| **Build** | VCF x annotations -> `sparse_g.bin` + `variants.parquet` + `membership.parquet` | SHA-256(VCF content, annotation content) |
| **Score cache** | Null model -> U, K per gene | (store key, trait, covariates) |
| **Test** | Slice cached U/K per mask -> Burden, SKAT, ACAT-V -> omnibus | (mask predicate, MAF cutoff) |

Interactive results: Plotly.js summary with Manhattan, QQ, and volcano plots.

See [docs/storage.md](docs/storage.md) for the storage format and [docs/validation.md](docs/validation.md) for validation against the R STAARpipeline.

## Data packs

| Pack | Size | Description |
|------|------|-------------|
| **favor-base** | 200 GB | 40 curated annotation columns including pathogenicity, frequency, clinical, conservation, regulatory, aPC STAAR channels |
| **favor-full** | 508 GB | All 54 annotation columns including dbnsfp, ENCODE, MaveDB, COSMIC |
| eqtl | 3 GB | GTEx v10 eQTL/sQTL/apaQTL, 50 tissues, SuSiE fine-mapped |
| sc-eqtl | 48 GB | Single-cell eQTL: OneK1K, DICE, PsychENCODE |
| regulatory | 18 GB | cCRE tissue signals, chromatin states, accessibility |
| enhancer-gene | 12 GB | ABC, EPIraction, rE2G, EpiMap, CRISPRi |
| tissue-scores | 5 GB | ChromBPNet, allelic imbalance |

## Roadmap

Tracked in [GitHub Issues](https://github.com/vineetver/cohort-cli/issues) with milestones:

| Milestone | Focus |
|-----------|-------|
| [v0.2.0](https://github.com/vineetver/cohort-cli/milestone/1) | STAAR hardening: GRM, multi-VCF, AI-STAAR, MultiSTAAR, performance |
| [v0.3.0](https://github.com/vineetver/cohort-cli/milestone/2) | MetaSTAAR: cross-biobank meta-analysis |
| [v0.4.0](https://github.com/vineetver/cohort-cli/milestone/3) | Variant interpretation: scoring, fine-mapping, colocalization, V2G |
| [v1.0.0](https://github.com/vineetver/cohort-cli/milestone/4) | Nextflow orchestration, provenance, QC |

## Citation

COHORT CLI implements the [STAAR](https://github.com/xihaoli/STAARpipeline) framework and the [FAVOR](https://favor.genohub.org) annotation database. If you use this tool, please cite:

> Li Z\*, Li X\*, Zhou H, et al. **A framework for detecting noncoding rare variant associations of large-scale whole-genome sequencing studies.** *Nature Methods*, 19(12), 1599-1611 (2022). [DOI: 10.1038/s41592-022-01640-x](https://doi.org/10.1038/s41592-022-01640-x)

> Li X\*, Li Z\*, Zhou H, et al. **Dynamic incorporation of multiple in silico functional annotations empowers rare variant association analysis of large whole-genome sequencing studies at scale.** *Nature Genetics*, 52(9), 969-983 (2020). [DOI: 10.1038/s41588-020-0676-4](https://doi.org/10.1038/s41588-020-0676-4)

> Zhou H, Verma V, Li X, et al. **FAVOR 2.0: A reengineered functional annotation of variants online resource for interpreting genomic variation.** *Nucleic Acids Research*, 54(D1), D1405-D1414 (2026). [DOI: 10.1093/nar/gkaf1217](https://doi.org/10.1093/nar/gkaf1217)

> Li TC, Zhou H, Verma V, et al. **FAVOR-GPT: a generative natural language interface to whole genome variant functional annotations.** *Bioinformatics Advances*, 4(1), vbae143 (2024). [DOI: 10.1093/bioadv/vbae143](https://doi.org/10.1093/bioadv/vbae143)

## License

GPL-3.0
