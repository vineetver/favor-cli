use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

use crate::output::Format;

#[derive(Debug, Clone, ValueEnum)]
pub enum GenomeBuild {
    Hg38,
    Hg19,
}

#[derive(Parser)]
#[command(
    name = "favor",
    version,
    about = "FAVOR - Functional Annotation of Variants Online Resource",
    long_about = "Genomic variant annotation, enrichment, and analysis toolkit.\n\n\
        Human mode (default in terminal): colored output, progress bars, interactive prompts.\n\
        Machine mode (--format json or FAVOR_MACHINE=1): structured JSON for LLM agents."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    /// Output format: auto (detect tty), json (machine), human (forced)
    #[arg(long, global = true, default_value = "auto")]
    pub format: Format,

    /// Validate inputs and emit a JSON execution plan without doing real work
    #[arg(long, global = true, default_value_t = false)]
    pub dry_run: bool,

    /// Override auto-detected thread count (otherwise SLURM/cgroup/nproc)
    #[arg(long, global = true)]
    pub threads: Option<usize>,
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

    /// Interactive setup wizard - configure data tier, paths, resources
    Setup {
        /// Environment type (for machine mode): hpc or workstation
        #[arg(long)]
        environment: Option<String>,

        /// Default memory budget (for machine mode): e.g. 16GB, 64GB
        #[arg(long)]
        memory_budget: Option<String>,
    },

    /// Manage annotation data (pull, status)
    Data {
        #[command(subcommand)]
        action: DataAction,
    },

    /// Ingest variants: auto-detect format, normalize, output parquet
    Ingest {
        /// Input file(s) or directory (VCF, TSV, CSV, parquet)
        #[arg(required = true, num_args = 1..)]
        inputs: Vec<PathBuf>,

        /// Output directory [default: <input_stem>.ingested/]
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Emit a SQL script instead of running it (for editing)
        #[arg(long)]
        emit_sql: bool,

        /// Genome build of input: hg38 or hg19 [default: auto-detect]
        #[arg(long)]
        build: Option<GenomeBuild>,
    },

    /// Annotate variants against favor-base or favor-full
    Annotate {
        /// Input variant set directory (from ingest)
        input: PathBuf,

        /// Output directory for annotated variants [default: <input_stem>.annotated/]
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Use favor-full instead of favor-base
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
        /// Multi-sample VCF with genotypes (.vcf or .vcf.gz)
        #[arg(long)]
        genotypes: PathBuf,

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

        /// Run per-variant individual score tests
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

        /// Known loci file for conditional analysis (one chr:pos:ref:alt per line)
        #[arg(long)]
        known_loci: Option<PathBuf>,

        /// Export summary statistics for MetaSTAAR instead of running tests
        #[arg(long)]
        emit_sumstats: bool,

        /// Force rebuild the genotype store even if a valid cache exists
        #[arg(long)]
        rebuild_store: bool,

        /// Skip genotype store entirely (no caching, one-off analysis)
        #[arg(long)]
        no_store: bool,

        /// Use feedback loop: load priors from prior runs, record results for future runs
        #[arg(long)]
        adaptive: bool,

        /// Column name mapping for phenotype file (key=value pairs, e.g. id=IID,trait=BMI)
        #[arg(long, value_delimiter = ',')]
        column_map: Vec<String>,

        /// Output directory
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
        /// Download favor-full (508GB) instead of configured tier
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
