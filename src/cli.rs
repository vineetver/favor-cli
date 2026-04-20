use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

use crate::output::Format;

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum ConditionalModel {
    Homogeneous,
    Heterogeneous,
}

#[derive(Debug, Clone, ValueEnum)]
pub enum GenomeBuild {
    Hg38,
    Hg19,
}

#[derive(Parser)]
#[command(
    name = "favor",
    version,
    about = "FAVOR CLI - genomic variant annotation and rare-variant association testing",
    long_about = "Genomic variant annotation, enrichment, and analysis toolkit.\n\n\
        Human mode (default in terminal): colored output, progress bars, interactive prompts.\n\
        Machine mode (--format json or FAVOR_MACHINE=1): structured JSON for LLM agents."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Option<Command>,

    /// Output format: auto (detect tty), json (machine), human (forced)
    #[arg(long, global = true, default_value = "auto")]
    pub format: Format,

    /// Validate inputs and emit a JSON execution plan without doing real work
    #[arg(long, global = true, default_value_t = false)]
    pub dry_run: bool,

    /// Override auto-detected thread count (otherwise SLURM/cgroup/nproc)
    #[arg(long, global = true)]
    pub threads: Option<usize>,

    /// Path to the cohort store root. Falls back to FAVOR_STORE, then a
    /// parent-walk for `.cohort/`, then `<cwd>/.cohort/`.
    #[arg(long, global = true)]
    pub store_path: Option<PathBuf>,
}

#[derive(Subcommand)]
#[allow(clippy::large_enum_variant)]
pub enum Command {
    /// Initialize a project directory with agent context files
    Init {
        /// Project directory [default: current directory]
        path: Option<PathBuf>,

        /// Overwrite existing files without prompting
        #[arg(long)]
        force: bool,
    },

    /// Configure data tier, paths, resources. Flag-driven; the previous
    /// interactive TUI wizard has been removed.
    Setup {
        /// Annotation root directory (where FAVOR / tissue data lives)
        #[arg(long)]
        root: PathBuf,

        /// Annotation tier: base or full
        #[arg(long)]
        tier: String,

        /// Comma-separated optional pack ids to install
        #[arg(long, value_delimiter = ',')]
        packs: Vec<String>,

        /// Environment type: hpc or workstation
        #[arg(long)]
        environment: Option<String>,

        /// Default memory budget: e.g. 16GB, 64GB
        #[arg(long)]
        memory_budget: Option<String>,
    },

    /// Manage annotation data (pull, status)
    Data {
        #[command(subcommand)]
        action: DataAction,
    },

    /// Manage attached annotation databases (refs.toml)
    Annotation {
        #[command(subcommand)]
        action: AnnotationAction,
    },

    /// Cohort store maintenance
    Store {
        #[command(subcommand)]
        action: StoreAction,
    },

    /// Ingest variants: auto-detect format, normalize, output parquet.
    /// For a multi-sample VCF (biobank dosages, e.g. UKB/TOPMed/All of Us)
    /// pass `--annotations <annotated-set>` and the cohort genotype store is
    /// built under `.cohort/cohorts/<id>/` instead of a variant-only set.
    Ingest {
        /// Input file(s) or directory (VCF, TSV, CSV, parquet)
        #[arg(required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,

        /// Output directory [default: <input_stem>.ingested/].
        /// Ignored for multi-sample VCFs (cohort path uses `.cohort/cohorts/<id>/`).
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Emit a SQL script instead of running it (for editing)
        #[arg(long)]
        emit_sql: bool,

        /// Genome build of input: hg38 or hg19 [default: auto-detect]
        #[arg(long)]
        build: Option<GenomeBuild>,

        /// Annotated variant set from `favor annotate` (required for
        /// multi-sample VCF input — provides STAAR weights + gene assignments)
        #[arg(long)]
        annotations: Option<PathBuf>,

        /// Cohort id under `.cohort/cohorts/` (multi-sample VCF only).
        /// Defaults to the VCF stem.
        #[arg(long)]
        cohort_id: Option<String>,

        /// Force rebuild of the cohort store even if a valid one exists
        #[arg(long)]
        rebuild: bool,

        /// Restrict ingest to these chromosomes. Accepts `22`, `1,2,3`,
        /// `1-22`, or combinations like `1-22,X,Y,MT`. Default: all standard
        /// chromosomes.
        #[arg(long, value_name = "SPEC")]
        chromosome: Option<crate::types::ChromosomeSet>,
    },

    /// Annotate variants against FAVOR base or full tier
    Annotate {
        /// Input variant set directory (from ingest)
        input: PathBuf,

        /// Output directory for annotated variants [default: <input_stem>.annotated/]
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Use FAVOR full tier instead of base
        #[arg(long)]
        full: bool,
    },

    /// Enrich annotated variants with tissue-specific data
    Enrich {
        /// Input annotated variant set directory
        input: PathBuf,

        /// Tissue name or group
        #[arg(long)]
        tissue: String,

        /// Output path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Variant interpretation - ACMG (coding) + regulatory framework (noncoding)
    Interpret {
        /// Input annotated variant set directory
        input: PathBuf,

        /// Target tissue for context-aware scoring
        #[arg(long)]
        tissue: Option<String>,

        /// Disease context for re-ranking
        #[arg(long)]
        disease: Option<String>,

        /// Output path
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// STAAR rare variant association testing
    Staar {
        /// Multi-sample VCF(s) with genotypes (.vcf or .vcf.gz). Accepts
        /// multiple files, globs, or a directory of per-chromosome VCFs.
        /// Required unless `--cohort <id>` is given to load a pre-built
        /// cohort store.
        #[arg(long, num_args = 1..)]
        genotypes: Vec<PathBuf>,

        /// Pre-built cohort id (from `favor ingest <vcf> --annotations ... --cohort-id <id>`).
        /// Skips probe + rebuild — trusts the manifest at `.cohort/cohorts/<id>/`.
        #[arg(long)]
        cohort: Option<String>,

        /// Phenotype file (TSV with sample_id as first column)
        #[arg(long)]
        phenotype: PathBuf,

        /// Trait column name(s) in phenotype file (comma-separated for MultiSTAAR)
        #[arg(long, value_delimiter = ',')]
        trait_name: Vec<String>,

        /// Covariate columns (comma-separated, e.g. age,sex,PC1,PC2)
        #[arg(long, value_delimiter = ',')]
        covariates: Vec<String>,

        /// Annotated variant set from `favor annotate` (provides aPC weights + gene assignments)
        #[arg(long)]
        annotations: Option<PathBuf>,

        /// Mask types to test: coding, noncoding, sliding-window, scang [default: coding]
        #[arg(long, value_delimiter = ',', default_value = "coding")]
        masks: Vec<String>,

        /// MAF cutoff for rare variants [default: 0.01]
        #[arg(long, default_value = "0.01")]
        maf_cutoff: f64,

        /// Sliding window size in bp [default: 2000]
        #[arg(long, default_value = "2000")]
        window_size: u32,

        /// Run per-variant individual score tests. Gaussian + unrelated samples only.
        #[arg(long)]
        individual: bool,

        /// Use saddlepoint approximation for imbalanced binary traits
        #[arg(long)]
        spa: bool,

        /// AI-STAAR: population group column in phenotype file (activates ancestry-informed weights)
        #[arg(long)]
        ancestry_col: Option<String>,

        /// AI-STAAR ensemble base tests B
        #[arg(long, default_value = "5")]
        ai_base_tests: usize,

        /// AI-STAAR ensemble weight RNG seed
        #[arg(long, default_value = "7590")]
        ai_seed: u64,

        /// SCANG minimum variants per window [default: 40]
        #[arg(long, default_value = "40")]
        scang_lmin: usize,

        /// SCANG maximum variants per window [default: 300]
        #[arg(long, default_value = "300")]
        scang_lmax: usize,

        /// SCANG step between window sizes [default: 10]
        #[arg(long, default_value = "10")]
        scang_step: usize,

        /// Kinship matrix file(s) for mixed-model analysis (TSV: sample_i sample_j kinship).
        /// Pass multiple files comma-separated to combine matrices (e.g. population GRM + pedigree).
        #[arg(long, value_delimiter = ',')]
        kinship: Vec<PathBuf>,

        /// Phenotype column naming a categorical group for heteroscedastic residual variance.
        /// Each distinct level becomes a separate τ_e component in the AI-REML fit.
        #[arg(long)]
        kinship_groups: Option<String>,

        /// Phenotype column carrying per-sample time values for the
        /// longitudinal random-slope LMM (GMMAT `glmmkin.R:174-183`).
        /// Each input --kinship matrix expands into three variance
        /// components (intercept, intercept-slope covariance, slope).
        /// Requires --kinship; gaussian family only; single-trait only.
        #[arg(long)]
        random_slope: Option<String>,

        /// Known loci file for conditional analysis (one chr:pos:ref:alt per line)
        #[arg(long)]
        known_loci: Option<PathBuf>,

        /// Load a pre-fit null model from FVNULLM1 binary (skips the fit stage).
        /// Unrelated samples only; kinship-aware import lands once the
        /// sparse Cholesky factor is serializable.
        #[arg(long, value_name = "PATH")]
        null_model: Option<PathBuf>,

        /// Export summary statistics for MetaSTAAR instead of running tests
        #[arg(long)]
        emit_sumstats: bool,

        /// Force rebuild the genotype store even if a valid cache exists
        #[arg(long)]
        rebuild_store: bool,

        /// Cohort id under the store root. Defaults to the genotype VCF stem.
        #[arg(long)]
        cohort_id: Option<String>,

        /// Skip the genotype store (no-op; reserved for the streaming path).
        #[arg(long)]
        no_store: bool,

        /// Load priors from prior runs and record this run's results
        /// (no-op; reserved for the adaptive-scoring loop).
        #[arg(long)]
        adaptive: bool,

        /// Column name mapping for phenotype file (key=value pairs, e.g. id=IID,trait=BMI)
        #[arg(long, value_delimiter = ',')]
        column_map: Vec<String>,

        /// Output directory
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Compute a sparse GRM and PCA scores from a cohort + KING .seg output.
    /// Outputs are cached under `.cohort/cache/grm/<cohort>/` and consumed
    /// by `favor staar --kinship <grm.tsv> --pca-covariates <pca_scores.tsv>`.
    Grm {
        /// Pre-built cohort id (under the store root).
        #[arg(long)]
        cohort: String,

        /// KING .seg IBD output file (external tool, not reimplemented).
        #[arg(long)]
        king_seg: PathBuf,

        /// Maximum relatedness degree (default 4 = up to 3rd cousins).
        #[arg(long, default_value = "4")]
        degree: u8,

        /// Number of PCA components (default 20).
        #[arg(long, default_value = "20")]
        n_pcs: usize,

        /// SNP block size for memory control (default 5000).
        #[arg(long, default_value = "5000")]
        block_size: usize,

        /// Output directory (default: cohort GRM cache path)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Forward-selection LD pruning on conditional score-test p-values
    #[command(name = "ld-prune")]
    LdPrune {
        /// Pre-built cohort id (under the store root).
        #[arg(long)]
        cohort: String,

        /// Phenotype file (TSV with sample_id as first column)
        #[arg(long)]
        phenotype: PathBuf,

        /// Trait column name in the phenotype file
        #[arg(long)]
        trait_name: String,

        /// Covariate columns (comma-separated, e.g. age,sex,PC1,PC2)
        #[arg(long, value_delimiter = ',')]
        covariates: Vec<String>,

        /// Candidate variants file. Tab-delimited or colon-delimited with
        /// four fields per row: CHR POS REF ALT. `#`-prefixed lines skipped.
        #[arg(long)]
        variants: PathBuf,

        /// Minor allele frequency floor for candidates (default 0.01)
        #[arg(long, default_value = "0.01")]
        maf_cutoff: f64,

        /// Conditional p-value threshold at which forward selection stops
        /// (default 1e-4, matches STAARpipeline LD_pruning).
        #[arg(long, default_value = "1e-4")]
        cond_p_thresh: f64,

        /// Column name mapping for phenotype file (key=value pairs)
        #[arg(long, value_delimiter = ',')]
        column_map: Vec<String>,

        /// Output TSV path (default: <cohort>.ld_pruned.tsv)
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Meta-analysis of STAAR across studies (MetaSTAAR)
    #[command(name = "meta-staar")]
    MetaStaar {
        /// Directories containing per-study summary statistics (from --emit-sumstats)
        #[arg(long, value_delimiter = ',')]
        studies: Vec<PathBuf>,

        /// Mask types to test: coding, noncoding, sliding-window [default: coding]
        #[arg(long, value_delimiter = ',', default_value = "coding")]
        masks: Vec<String>,

        /// MAF cutoff for meta-analysis [default: 0.01]
        #[arg(long, default_value = "0.01")]
        maf_cutoff: f64,

        /// Sliding window size in bp [default: 2000]
        #[arg(long, default_value = "2000")]
        window_size: u32,

        /// Known loci for conditional analysis (one chr:pos:ref:alt per line)
        #[arg(long)]
        known_loci: Option<PathBuf>,

        /// Conditional model: homogeneous (shared effects across studies) or
        /// heterogeneous (study-specific effects)
        #[arg(long, default_value = "homogeneous")]
        conditional_model: ConditionalModel,

        /// Output directory
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Show schema of annotation tables (for LLM agents)
    Schema {
        /// Table name: "base", "full", or a tissue pack name
        table: Option<String>,
    },

    /// List available analyses and data status (for LLM agents)
    Manifest,

    /// Remove favor binary and config
    Uninstall,
}

#[derive(Subcommand)]
pub enum DataAction {
    /// Download annotation data or add-on packs
    Pull {
        /// Download FAVOR full tier (508GB) instead of configured tier
        #[arg(long)]
        full: bool,

        /// Show what would be downloaded without downloading
        #[arg(long)]
        dry_run: bool,

        /// Skip download confirmation prompt
        #[arg(short, long)]
        yes: bool,

        /// Download a specific tissue add-on pack (e.g., eqtl, regulatory)
        #[arg(long)]
        pack: Option<String>,
    },

    /// Show installed data status
    Status,

    /// Verify local data against remote manifests
    Verify {
        /// Verify a specific pack (default: all installed)
        #[arg(long)]
        pack: Option<String>,

        /// Also verify SHA-256 checksums (slow — reads every byte)
        #[arg(long)]
        checksums: bool,
    },

    /// Generate v2 manifests and upload packs to MinIO (admin)
    Publish {
        /// Source data root directory containing annotation parquets
        #[arg(long)]
        source_root: PathBuf,

        /// Publish only this pack (default: all)
        #[arg(long)]
        pack: Option<String>,

        /// Generate manifests without uploading
        #[arg(long)]
        dry_run: bool,
    },
}

#[derive(Subcommand)]
pub enum AnnotationAction {
    /// Attach an external annotation database to refs.toml
    Attach {
        /// Alias to register (e.g. "ukb-base", "favor-full-2024")
        name: String,

        /// Annotation kind: favor-base, favor-full, or tissue
        #[arg(long)]
        kind: AttachKind,

        /// Path to the annotation directory
        #[arg(long)]
        path: PathBuf,
    },
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum AttachKind {
    FavorBase,
    FavorFull,
    Tissue,
}

#[derive(Subcommand)]
pub enum StoreAction {
    /// Garbage-collect orphaned cohort caches
    Gc,
}
