//! Genotype store: build, cache, probe, and manifest management.
//!
//! On first run the store is built from the VCF x annotation join. Subsequent
//! runs with the same VCF and annotations skip both steps entirely and load
//! pre-joined data via memory-mapped sparse_g.bin.
//!
//! Layout:
//! ```text
//! .genotype_store/
//!   manifest.json
//!   samples.txt
//!   chromosome={chr}/
//!     sparse_g.bin         # sparse genotype matrix (sample_id, variant_vcf) → dosage
//!     variants.parquet     # aligned metadata vectors, row i = variant_vcf i
//!     membership.parquet   # (variant_vcf, gene_name) many-to-many
//! ```

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::column;
use crate::data::{VariantSet, VariantSetKind};
use crate::engine::DfEngine;
use crate::error::FavorError;
use crate::ingest::{ColumnContract, ColumnRequirement};
use crate::output::Output;
use crate::resource::Resources;

#[derive(Serialize, Deserialize)]
pub struct StoreManifest {
    pub version: u32,
    pub key: String,
    pub n_samples: usize,
    pub n_variants: usize,
    pub chromosomes: Vec<ChromInfo>,
    pub created_at: String,
    pub favor_version: String,
}

#[derive(Serialize, Deserialize)]
pub struct ChromInfo {
    pub name: String,
    pub n_variants: usize,
}

/// Content-based file fingerprint: SHA-256 of (first 1MB + last 1MB + file size).
/// Survives path renames and HPC scratch migrations — only invalidates when
/// actual file content changes.
fn file_content_fingerprint(path: &Path) -> Result<Vec<u8>, FavorError> {
    use std::io::{Read as IoRead, Seek, SeekFrom};
    const CHUNK: u64 = 1024 * 1024;

    let mut f = File::open(path)
        .map_err(|e| FavorError::Resource(format!("open {}: {e}", path.display())))?;
    let size = f
        .metadata()
        .map_err(|e| FavorError::Resource(format!("stat {}: {e}", path.display())))?
        .len();

    let mut hasher = Sha256::new();
    hasher.update(size.to_le_bytes());

    // First 1MB
    let head = CHUNK.min(size) as usize;
    let mut buf = vec![0u8; head];
    f.read_exact(&mut buf)
        .map_err(|e| FavorError::Resource(format!("read {}: {e}", path.display())))?;
    hasher.update(&buf);

    // Last 1MB (skip if file <= 1MB since it's already fully read)
    if size > CHUNK {
        let tail_start = size.saturating_sub(CHUNK);
        f.seek(SeekFrom::Start(tail_start))
            .map_err(|e| FavorError::Resource(format!("seek {}: {e}", path.display())))?;
        let tail_len = (size - tail_start) as usize;
        buf.resize(tail_len, 0);
        f.read_exact(&mut buf)
            .map_err(|e| FavorError::Resource(format!("read {}: {e}", path.display())))?;
        hasher.update(&buf);
    }

    Ok(hasher.finalize().to_vec())
}

/// For annotation directories, fingerprint the meta.json file.
fn dir_fingerprint(path: &Path) -> Result<Vec<u8>, FavorError> {
    let meta_path = path.join("meta.json");
    if meta_path.exists() {
        return file_content_fingerprint(&meta_path);
    }
    // Fallback: hash the directory's own mtime+size (stat only)
    let meta = fs::metadata(path)
        .map_err(|e| FavorError::Resource(format!("stat {}: {e}", path.display())))?;
    let mtime = meta
        .modified()
        .unwrap_or(SystemTime::UNIX_EPOCH)
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();
    let mut hasher = Sha256::new();
    hasher.update(mtime.to_le_bytes());
    Ok(hasher.finalize().to_vec())
}

pub fn compute_key(vcf_path: &Path, annotations_path: &Path) -> Result<String, FavorError> {
    let vcf_fp = file_content_fingerprint(vcf_path)?;
    let ann_fp = if annotations_path.is_dir() {
        dir_fingerprint(annotations_path)?
    } else {
        file_content_fingerprint(annotations_path)?
    };
    let mut hasher = Sha256::new();
    hasher.update(&vcf_fp);
    hasher.update(&ann_fp);
    Ok(format!("{:x}", hasher.finalize()))
}

pub struct StoreProbe {
    pub store_dir: PathBuf,
    /// `Some` only when a valid, key-matching, fully-materialized store exists.
    pub manifest: Option<StoreManifest>,
}

pub fn probe(store_dir: &Path, vcf_path: &Path, annotations_path: &Path) -> StoreProbe {
    let store_dir = store_dir.to_path_buf();
    let miss = || StoreProbe { store_dir: store_dir.clone(), manifest: None };

    let Ok(s) = fs::read_to_string(store_dir.join("manifest.json")) else { return miss(); };
    let Ok(manifest) = serde_json::from_str::<StoreManifest>(&s) else { return miss(); };
    let Ok(key) = compute_key(vcf_path, annotations_path) else { return miss(); };

    if manifest.key != key || manifest.version != 3 {
        return miss();
    }

    for ci in &manifest.chromosomes {
        let sparse_g = store_dir.join(format!("chromosome={}/sparse_g.bin", ci.name));
        let variants = store_dir.join(format!("chromosome={}/variants.parquet", ci.name));
        if !sparse_g.exists() || !variants.exists() {
            return miss();
        }
    }

    StoreProbe { store_dir, manifest: Some(manifest) }
}

pub const STAAR_ANNOTATION_COLUMNS: &[ColumnRequirement] = &[
    ColumnRequirement {
        name: "gencode",
        source: "FAVOR full annotations",
        used_by: "gene/region/consequence extraction",
    },
    ColumnRequirement {
        name: "main",
        source: "FAVOR full annotations",
        used_by: "CADD score extraction",
    },
    ColumnRequirement {
        name: "cage",
        source: "FAVOR full annotations",
        used_by: "regulatory mask predicates",
    },
    ColumnRequirement {
        name: "apc",
        source: "FAVOR full annotations",
        used_by: "11 annotation weight channels",
    },
    ColumnRequirement {
        name: "dbnsfp",
        source: "FAVOR full annotations",
        used_by: "REVEL score for masks",
    },
    ColumnRequirement {
        name: "ccre",
        source: "FAVOR full annotations",
        used_by: "cCRE regulatory masks",
    },
    ColumnRequirement {
        name: "linsight",
        source: "FAVOR full annotations",
        used_by: "LINSIGHT annotation weight",
    },
    ColumnRequirement {
        name: "fathmm_xf",
        source: "FAVOR full annotations",
        used_by: "FATHMM-XF annotation weight",
    },
];

pub struct GenoStoreResult {
    pub sample_names: Vec<String>,
    pub store_dir: PathBuf,
    pub manifest: StoreManifest,
    pub engine: DfEngine,
    pub geno_output_dir: Option<PathBuf>,
}

pub fn setup_resources(out: &dyn Output) -> Result<Resources, FavorError> {
    let resources = Resources::detect_configured();

    out.status(&format!(
        "STAAR: {} memory, {} threads ({})",
        resources.memory_human(),
        resources.threads,
        resources.environment()
    ));

    Ok(resources)
}

pub fn build_or_load_store(
    genotypes: &Path,
    annotations: &Path,
    store_dir: &Path,
    geno_staging_dir: &Path,
    rebuild_store: bool,
    res: &Resources,
    out: &dyn Output,
) -> Result<GenoStoreResult, FavorError> {
    let probe_result = if rebuild_store {
        StoreProbe { store_dir: store_dir.to_path_buf(), manifest: None }
    } else {
        probe(store_dir, genotypes, annotations)
    };

    if let Some(manifest) = probe_result.manifest {
        out.status(&format!(
            "Step 1/4: Using cached genotype store ({} variants x {} samples)",
            manifest.n_variants, manifest.n_samples,
        ));

        let sample_names = read_sample_names(&probe_result.store_dir)?;
        let engine = DfEngine::new(res)?;

        return Ok(GenoStoreResult {
            sample_names,
            store_dir: probe_result.store_dir,
            manifest,
            engine,
            geno_output_dir: None,
        });
    }

    out.status("Step 1/4: Building genotype store...");

    out.status("  Extracting genotypes from VCF...");
    let geno = crate::staar::genotype::extract_genotypes(
        genotypes,
        geno_staging_dir,
        res.memory_bytes,
        res.threads,
        out,
    )?;

    out.status("  Joining genotypes with annotations...");
    let engine = run_annotation_join(annotations, &geno, res)?;

    // Deduplicate: same (chromosome, position, ref, alt) from multi-allelic splits or overlapping VCFs.
    let dup_count = engine.query_scalar(&column::dedup_count_sql())?;
    if dup_count > 0 {
        out.warn(&format!(
            "  {dup_count} duplicate variants found — keeping first occurrence."
        ));
        engine.execute(&column::dedup_sql())?;
        engine.execute("DROP TABLE _rare_all")?;
        engine.execute("ALTER TABLE _rare_dedup RENAME TO _rare_all")?;
    }

    let n_all = engine.query_scalar("SELECT COUNT(*) FROM _rare_all")?;
    let n_genes = engine.query_scalar(&column::gene_count_sql())?;
    out.status(&format!(
        "  {} annotated variants, {} genes",
        n_all, n_genes
    ));

    if n_all == 0 {
        return Err(FavorError::Analysis(format!(
            "No variants found after joining genotypes ({}) with annotations ({}). \
             Check that both use the same genome build and allele normalization.",
            genotypes.display(),
            annotations.display(),
        )));
    }

    let manifest = build(&engine, &geno, store_dir, genotypes, annotations, out)?;

    Ok(GenoStoreResult {
        sample_names: geno.sample_names.clone(),
        store_dir: store_dir.to_path_buf(),
        manifest,
        engine,
        geno_output_dir: Some(geno.output_dir),
    })
}

fn read_sample_names(store_dir: &Path) -> Result<Vec<String>, FavorError> {
    let path = store_dir.join("samples.txt");
    let content = std::fs::read_to_string(&path)
        .map_err(|e| FavorError::Resource(format!("Read {}: {e}", path.display())))?;
    Ok(content.lines().map(|s| s.to_string()).collect())
}

pub fn run_annotation_join(
    annotations: &Path,
    geno: &crate::staar::genotype::GenotypeResult,
    res: &Resources,
) -> Result<DfEngine, FavorError> {
    let engine = DfEngine::new(res)?;

    let ann_vs = VariantSet::open(annotations)?;

    // Check required annotation columns via table schema (not a query —
    // LIMIT 0 can return zero batches in DataFusion, losing the schema).
    let ann_cols: Vec<String> = {
        engine.register_parquet_dir("_ann_check", ann_vs.root())?;
        let cols = engine.table_columns("_ann_check")?;
        let _ = engine.execute("DROP TABLE IF EXISTS _ann_check");
        cols
    };

    let contract = ColumnContract {
        command: "staar (annotation join)",
        required: STAAR_ANNOTATION_COLUMNS,
    };
    let missing = contract.check(&ann_cols);
    if !missing.is_empty() {
        let tier_hint = match ann_vs.kind() {
            Some(VariantSetKind::Annotated {
                tier: crate::config::Tier::Base,
            }) => " Your data was annotated with base tier. Re-run: `favor annotate --full`.",
            _ => " Re-run: `favor annotate --full`.",
        };
        return Err(FavorError::DataMissing(format!(
            "Missing annotation columns in {}:\n{}\n\
             STAAR requires favor-full annotations.{}",
            annotations.display(),
            ColumnContract::format_missing(&missing),
            tier_hint,
        )));
    }

    let geno_dir = &geno.output_dir;
    engine.register_parquet_sorted("_genotypes", geno_dir, &["position", "ref", "alt"])?;
    engine.register_parquet_sorted(
        "_annotations",
        ann_vs.root(),
        &["position", "ref_vcf", "alt_vcf"],
    )?;

    // The store keeps all variants; MAF filtering happens at read time.
    engine.execute(&column::annotation_join_sql())?;

    Ok(engine)
}

pub fn build(
    engine: &DfEngine,
    geno_result: &super::genotype::GenotypeResult,
    store_dir: &Path,
    vcf_path: &Path,
    annotations_path: &Path,
    out: &dyn Output,
) -> Result<StoreManifest, FavorError> {
    if store_dir.exists() {
        fs::remove_dir_all(store_dir)
            .map_err(|e| FavorError::Resource(format!("Clean store dir: {e}")))?;
    }
    fs::create_dir_all(store_dir).map_err(|e| {
        FavorError::Resource(format!("Cannot create '{}': {e}", store_dir.display()))
    })?;

    let key = compute_key(vcf_path, annotations_path)?;
    let n_samples = geno_result.sample_names.len();

    let samples_path = store_dir.join("samples.txt");
    fs::write(&samples_path, geno_result.sample_names.join("\n")).map_err(|e| {
        FavorError::Resource(format!("Cannot write '{}': {e}", samples_path.display()))
    })?;

    let chroms = engine.query_strings(&column::distinct_chroms_sql())?;

    let mut chrom_infos: Vec<ChromInfo> = Vec::new();
    let mut total_variants = 0usize;
    let mut total_genes = 0usize;
    let mut total_carriers = 0u64;

    for chrom in &chroms {
        out.status(&format!("  Building sparse genotype store chr{chrom}..."));

        let chrom_dir = store_dir.join(format!("chromosome={chrom}"));
        fs::create_dir_all(&chrom_dir).map_err(|e| {
            FavorError::Resource(format!("Cannot create '{}': {e}", chrom_dir.display()))
        })?;

        let geno_path = format!(
            "{}/chromosome={chrom}/data.parquet",
            geno_result.output_dir.display()
        );
        let stats = crate::staar::sparse_g_writer::build_chromosome(
            engine, &geno_path, chrom, &chrom_dir, n_samples, out,
        )?;

        if stats.n_variants == 0 {
            continue;
        }

        chrom_infos.push(ChromInfo {
            name: chrom.clone(),
            n_variants: stats.n_variants,
        });
        total_variants += stats.n_variants;
        total_genes += stats.n_genes;
        total_carriers += stats.total_carriers;
    }

    let manifest = StoreManifest {
        version: 3,
        key,
        n_samples,
        n_variants: total_variants,
        chromosomes: chrom_infos,
        created_at: now_string(),
        favor_version: env!("CARGO_PKG_VERSION").to_string(),
    };

    let manifest_json = serde_json::to_string_pretty(&manifest)
        .map_err(|e| FavorError::Resource(format!("Manifest serialize: {e}")))?;
    fs::write(store_dir.join("manifest.json"), manifest_json)?;

    out.success(&format!(
        "Genotype store: {} variants x {} samples, {} genes, {} carriers",
        total_variants, n_samples, total_genes, total_carriers,
    ));

    Ok(manifest)
}

fn now_string() -> String {
    let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default();
    format!("{}", d.as_secs())
}
