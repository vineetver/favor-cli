# COHORT CLI — Agent Guide

Machine interface: always use `--format json`. Use `--dry-run` before committing resources.

## What this tool does

COHORT is an end-to-end genomic variant analysis toolkit. It takes called variants,
annotates them with functional predictions, overlays tissue-specific regulatory data,
runs association testing (rare-variant and common-variant), and interprets results
down to target gene, mechanism, and cell type.

```
                         ingest → annotate → enrich
                              /                \
               rare-variant lane            common-variant lane
               staar (association)          gwas (summary stats)
               prs (polygenic scores)       coloc (eQTL overlap)
                              \                /
                            interpret
                 (variant → gene → mechanism → cell type)
```

## Flags

| Flag | Description |
|------|-------------|
| `--format json` | Always use. stdout=result, stderr=status. |
| `--dry-run` | Estimate memory, validate inputs. No real work. |
| `-o <path>` | Override output path |

## Memory

COHORT adapts to available memory. DataFusion commands (annotate, enrich, ingest) work with any
allocation. STAAR needs the genotype matrix in RAM per-chromosome — use `--dry-run` first:

```bash
cohort staar --dry-run --format json ...
# → {"memory": {"recommended": "64G"}}
srun --mem=64G -c 8 cohort staar --format json ...
```

## Output conventions

- **Parquet**: main result data
- **`.meta.json`**: parameters, counts, data versions
- **stdout**: structured JSON result summary
- **stderr**: progress/status (JSON lines in machine mode)
- **Exit codes**: 0=ok, 1=input, 2=data missing, 3=resource, 4=analysis, 5=internal

## Querying results

```sql
-- Annotated variants
SELECT vid, gencode.genes[1] AS gene, gencode.consequence, main.cadd.phred
FROM 'annotated.parquet' WHERE main.cadd.phred > 20;

-- STAAR significant genes
SELECT * FROM 'staar_results/coding_pLoF_missense.parquet'
WHERE "STAAR-O" < 2.5e-6 ORDER BY "STAAR-O";

-- Enriched: join tissue data by vid
SELECT a.*, e.* FROM 'enriched/annotated.parquet' a
JOIN 'enriched/eqtl.parquet' e ON a.vid = e.vid
WHERE e.tissue_name LIKE '%Liver%';
```
