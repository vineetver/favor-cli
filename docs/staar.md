# STAAR

Takes a fitted null model and a sparse genotype slice, runs a score-test
battery per gene or window, writes per-mask parquet results.

```text
  phenotype + covariates (+ kinship if related samples)
              │
              │  fit null model   (expensive, once per trait)
              ▼
        null model
        ├─ residuals  r     what the covariates don't explain in y
        └─ projection P̂    how to measure a genotype's effect under the null
              │
              │  +  sparse genotype store (from favor ingest)
              │
              │  score test      (cheap, once per gene × mask)
              ▼
        per-gene parquet results
```

Fit the null once. Reuse it for thousands of gene × mask score tests.
Genotypes only enter at the score step.

> [Back to README](../README.md) · [Genotype store](storage.md) · [Validation](validation.md) · [Divergences](statistical-divergences.md)

## How The Pieces Tie Together

Three independent choices, picked once per run. They compose:

```text
null model kind       ×   score mode        ×   masks
───────────────           ──────────            ──────────────
Glm           (no K)      Standard              coding
Logistic      (no K)      Spa       (binary)    noncoding
KinshipReml   (with K)    AiStaar   (ancestry)  sliding-window
KinshipPql    (with K)                          scang
```

- **null model kind** — trait numeric or yes/no × samples unrelated or related.
- **score mode** — how to turn `U` into a p-value.
- **masks** — which variants enter each test.

## Pipeline Stages

`src/staar/pipeline.rs`. One resumable `run.json` tracks each stage.

```text
Validate          check annotation schema, 11 STAAR weight channels
EnsureStore       probe / build / load sparse genotype store
LoadPhenotype     read pheno.tsv, align to samples.txt, infer trait type
FitNullModel      dispatch on kind → r, P̂ pieces + disk cache
EmitSumstats      (RunMode=EmitSumstats) write U/K segments, early exit
EnsureScoreCache  fill chrom-wide U vector + per-gene K blocks
RunScoring        per chromosome × gene × mask → StaarResult
WriteResults      per-mask parquet + staar.meta.json
```

## Null Model Fit

```text
INPUT                      TRANSFORMATION                 OUTPUT
────────────────────       ─────────────────────────      ───────────────────
y    (n,)         ─┐                                       r     (n,)     residuals
X    (n, p+1)     ─┼─ fit (OLS / IRLS / REML / PQL) ───►   P̂    pieces   tuning for score test
K    (n, n) opt   ─┘                                       h²    scalar   heritability (if K)
```

`y` = phenotype per sample. `X` = intercept column + covariates. `K` =
kinship matrix when samples are related.

Worked example (fake LDL). Left = input, middle = what fit does, right = output `r`:

```text
sample age sex  LDL         fit finds Ŷ = f(age, sex)      r
──────── ─── ── ────        ──────────────────────────     ───
  0     45  M   130         expected 125 ··················+5
  1     32  F    90         expected  95 ··················-5
  2     61  M   165         expected 140 ··················+25
  3     29  F   105         expected 100 ··················+5
  4     54  M   148         expected 135 ··················+13
```

Genotypes never enter. `r` is the only per-sample output the score test reads.

Which fitter runs:

```text
trait numeric?   samples related?   fitter         under-the-hood
──────────────   ──────────────     ───────────    ─────────────────────
yes              no                 Glm            OLS via QR
yes              yes (kinship)      KinshipReml    REML, sparse Cholesky
yes/no           no                 Logistic       IRLS logistic
yes/no           yes (kinship)      KinshipPql     penalized quasi-likelihood
```

### Speed

Linear algebra runs on [`faer`](https://github.com/sarah-quinones/faer-rs),
a SIMD-vectorized BLAS-style backend. IRLS converges in a handful of
iterations on exome-scale cohorts. Kinship fits use sparse Cholesky
(`faer::sparse::Llt`) when `K` is block-diagonal, Hutchinson stochastic
trace for log-det, and Takahashi's formula for the partial inverse the
score test needs. `Glm` / `Logistic` fits cache to disk keyed by hash of
(phenotype, covariates, kinship). Kinship fits don't cache in v0.2.

## Score Test

```text
INPUT                         TRANSFORMATION              OUTPUT
────────────                  ───────────────────         ──────────────────
r     (n,)       from null  ─┐                            U    (m,)      signal
G_S   (n, m)     mask slice  ├─ U = G_Sᵀ r ────────►      K    (m, m)   noise floor
P̂   pieces      from null  ─┤  K = G_Sᵀ P̂ G_S ────►
w     (11, m)    annot wts   ─┘  → 3 × 2 × 11 tests       66 weighted p-values
                                 → Cauchy combine          STAAR-O (1 per gene)
```

Reading `U`, `K`: think signal and noise. `U[v]` = how much variant `v`
lines up with the residuals. Big `|U|` = suspicious. `K` = how much `U`
could wiggle by chance given the covariates already used.

### Annotation channel catalog

The eleven channels are fixed. Changing names, dropping one, or reordering
silently changes STAAR-O; `favor staar` validates the annotation parquet
against this catalog at startup and errors loud if any channel is missing.

```text
order   column name              source
─────   ──────────────────────   ────────────────────────────────
  1     w_cadd                   CADD-Phred (main annotation)
  2     w_linsight               LINSIGHT
  3     w_fathmm_xf              fathmm-XF
  4     w_apc_epi_active         aPC epigenetics, active mark
  5     w_apc_epi_repressed      aPC epigenetics, repressed mark
  6     w_apc_epi_transcription  aPC epigenetics, transcription
  7     w_apc_conservation       aPC conservation (v2)
  8     w_apc_protein_function   aPC protein function (v3)
  9     w_apc_local_nd           aPC local nucleotide diversity (v3)
 10     w_apc_mutation_density   aPC mutation density
 11     w_apc_tf                 aPC transcription-factor binding
```

Extra `w_*` columns are allowed — reads happen by name, so added channels
don't shift any existing index. Missing columns fail the preflight with
the column names called out.

### Test battery

Three test shapes × two β-weight shapes = six base tests. Each runs against
eleven annotation channels = sixty-six weighted tests. All Cauchy-combined
into `STAAR-O`.

```text
test      asks                                         good at
────────  ────────────────────────────────────────     ──────────────────────
Burden    do carriers shift the trait one direction?   all variants push same way
SKAT      is there more spread than chance?            mix of protective & harmful
ACAT-V    does any single variant stand out?           one bad actor in the set
```

### Score modes

```text
mode       when                                         what it changes
─────────  ──────────────────────────────────────       ──────────────────────────
Standard   default                                       reads cached U / K
Spa        rare binary trait, heavy case imbalance       saddlepoint p-value tails
AiStaar    ancestry-mixed cohort                         per-ancestry MAFs, ensemble
```

### Speed

`U = G_Sᵀ r` on sparse carrier lists is `O(total MAC)`, not `O(n × m)`.
For 200K samples and MAC = 5 per variant, that is a ~40,000× win over
dense. `K` blocks run through `faer` (SIMD). Chromosome-wide `U` and
per-gene `K` blocks are precomputed once into a score cache — rerunning
with a different mask reslices existing blocks instead of recomputing.
`rayon` parallelizes across genes within a chromosome.

## Masks

A mask is a predicate over annotated variants.

```text
INPUT                       TRANSFORMATION              OUTPUT
──────────────────────      ─────────────────────       ──────────────────
variants.parquet  ─┐        for each variant v:         variant_vcf list
  (gene's variants) ├──►    pred(v) AND maf(v) < cut ─► (feeds G_S slice)
mask predicate    ─┘
```

Categories:

```text
category         predicates expanded                                    how variants chosen
──────────────   ──────────────────────────────────────────────────    ──────────────────────
coding           pLoF, missense, disruptive_missense, pLoF_missense,   consequence + CADD + REVEL
                 synonymous, ptv, ptv_ds
noncoding        upstream, downstream, UTR, promoter_CAGE, promoter_   region_type + regulatory
                 DHS, enhancer_CAGE, enhancer_DHS, ncRNA                flags
sliding-window   fixed-width chunks (default 2 kb, step 2 kb)          positional
scang            variable-width chunks, L ∈ [lmin, lmax] variants      positional (data-adaptive)
custom           user BED (not yet wired)                              user-supplied
```

### Speed

`variants.parquet` is columnar, so a predicate filter is a column scan
over `consequence` / `cadd_phred` / regulatory flags. `membership.parquet`
is pre-sorted by `(gene, variant_vcf)` so "all variants in gene G" is a
range scan.

## Output Shape

One row per gene, one parquet per mask.

```text
INPUT                TRANSFORMATION              OUTPUT ROW (per gene)
──────────────────   ─────────────────────       ────────────────────────────────
U (m,)           ─┐                               ensembl_id, gene_symbol, chr,
K (m, m)          ├─ test battery ─────────────► start, end, n_variants, cMAC
w (11, m)        ─┘  + Cauchy combine            6 base tests       (f64, nullable)
                                                 66 weighted tests  (test × channel)
                                                 6 per-test omnibus (STAAR-B/S/A)
                                                 ACAT-O, STAAR-O
```

A gene row comes out null on a mask if `< 2` variants qualify or the
kernel degenerates.

## EmitSumstats and MetaSTAAR

`EmitSumstats` stops after the null fit and writes per-variant stats in a
form another study can pool later. Shape of what gets written per segment
(~500 kb genomic slice):

```text
INPUT                           TRANSFORMATION              OUTPUT (segment row)
────────────────────────        ─────────────────────       ──────────────────────────
r (n,)  from null          ─┐                               segment_id   u32
genotype slice (n, m)       ├─ U = G'r                      positions    [i32; m]
  (m variants in segment)   │  K = G' P̂ G  (lower tri)      refs         [utf8; m]
per-variant MAC, n_obs      ─┘                              alts         [utf8; m]
                                                            u_stat       [f64; m]
                                                            cov_lower    [f64; m(m+1)/2]
                                                            mac          [u32; m]
                                                            n_obs        [i32; m]
```

Lower triangle only — `K` is symmetric, so storing `m(m+1)/2` halves the
disk. Alleles are stored in canonical (lex-min) orientation so later
UNION-ALL across studies works without per-row flip logic at query time.

`favor meta-staar` pools these across studies:

```text
INPUT                              TRANSFORMATION                  OUTPUT
────────────────────────────       ─────────────────────────       ───────────────
study 1: (U₁, K₁) per segment  ─┐                                   U_meta  (M,)
study 2: (U₂, K₂) per segment  ─┼─ per canonical variant v:          K_meta  (M, M)
...                             │    U_meta[v] = Σ sign_s · U_s[v]    → same test
                                └─   K_meta[v,w] = Σ s_v·s_w · K_s    battery as
                                                                      single-study
                                                                      per-mask parquet
```

`sign_s ∈ {+1, −1}` flips when a study encoded the variant with the minor
allele on the other side. `M` is the union size over studies. The
canonical orientation in sumstats pushes all flip bookkeeping to one
lookup at merge time.

## MultiTrait

Joint null for `k` continuous traits sharing one covariate set.

```text
INPUT                       TRANSFORMATION              OUTPUT
──────────────────────      ─────────────────────       ──────────────────────
Y     (n, k)   traits    ─┐                              R       (n, k)    residuals
X     (n, p+1) covars    ─┼─ joint OLS / REML ─────►     Σ_res   (k, k)    residual cov
                                                         Σ_res⁻¹           (cached for SKAT)
```

`Σ_res` replaces `σ²` in the score-test kernel. Each gene × mask still
yields one joint p-value. v0.2: continuous traits only, no kinship.

## Caching

Three layers (see [storage](storage.md) for details):

```text
layer   artifact                         invalidated by
─────   ──────────────────────────────   ──────────────────────────────────
  1     sparse genotype store            VCF or annotation fingerprint change
  2     null model (Glm / Logistic)      phenotype or covariate change
  3     score cache (U, K blocks)        store rebuild, mask / MAF change
```

Swap the mask → layer 3 only. Swap the phenotype → 2 and 3. Rebuild store
→ all three.

## See Also

- [Genotype store](storage.md) — sparse G, variant index, score cache
- [Validation](validation.md) — R STAAR reference and tolerances
- [Statistical divergences](statistical-divergences.md) — known differences from R STAAR / SKAT
- [Performance](performance.md) — hot paths and memory picture
