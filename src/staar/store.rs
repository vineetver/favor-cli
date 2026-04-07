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

/// Why a `probe()` returned a miss. Surface this in `RunManifest` so
/// operators can tell "no store on disk" apart from "store exists but
/// content fingerprint changed" without grepping logs.
#[derive(Debug, Clone, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum ProbeReason {
    NoManifest,
    UnreadableManifest,
    ContentKeyChanged,
    SchemaVersionMismatch { found: u32, expected: u32 },
    MissingChromosomeArtifact { chromosome: String, file: String },
    FingerprintFailed,
}

pub struct StoreProbe {
    pub store_dir: PathBuf,
    /// `Some` only when a valid, key-matching, fully-materialized store exists.
    pub manifest: Option<StoreManifest>,
    /// Populated only when `manifest` is `None`.
    pub miss_reason: Option<ProbeReason>,
}

pub fn probe(store_dir: &Path, vcf_path: &Path, annotations_path: &Path) -> StoreProbe {
    let store_dir = store_dir.to_path_buf();
    let miss = |reason: ProbeReason| StoreProbe {
        store_dir: store_dir.clone(),
        manifest: None,
        miss_reason: Some(reason),
    };

    let Ok(s) = fs::read_to_string(store_dir.join("manifest.json")) else {
        return miss(ProbeReason::NoManifest);
    };
    let Ok(manifest) = serde_json::from_str::<StoreManifest>(&s) else {
        return miss(ProbeReason::UnreadableManifest);
    };
    let Ok(key) = compute_key(vcf_path, annotations_path) else {
        return miss(ProbeReason::FingerprintFailed);
    };

    if manifest.version != 3 {
        return miss(ProbeReason::SchemaVersionMismatch {
            found: manifest.version,
            expected: 3,
        });
    }
    if manifest.key != key {
        return miss(ProbeReason::ContentKeyChanged);
    }

    for ci in &manifest.chromosomes {
        let sparse_g = store_dir.join(format!("chromosome={}/sparse_g.bin", ci.name));
        let variants = store_dir.join(format!("chromosome={}/variants.parquet", ci.name));
        if !sparse_g.exists() {
            return miss(ProbeReason::MissingChromosomeArtifact {
                chromosome: ci.name.clone(),
                file: "sparse_g.bin".into(),
            });
        }
        if !variants.exists() {
            return miss(ProbeReason::MissingChromosomeArtifact {
                chromosome: ci.name.clone(),
                file: "variants.parquet".into(),
            });
        }
    }

    StoreProbe {
        store_dir,
        manifest: Some(manifest),
        miss_reason: None,
    }
}

/// Human-friendly summary of a probe miss for `out.status()` and run.json.
pub fn describe_miss(reason: &ProbeReason) -> String {
    match reason {
        ProbeReason::NoManifest => "no manifest.json on disk".into(),
        ProbeReason::UnreadableManifest => "manifest.json present but unparseable".into(),
        ProbeReason::ContentKeyChanged => {
            "VCF or annotations content fingerprint changed since last build".into()
        }
        ProbeReason::SchemaVersionMismatch { found, expected } => {
            format!("manifest schema v{found} on disk, this build expects v{expected}")
        }
        ProbeReason::MissingChromosomeArtifact { chromosome, file } => {
            format!("chromosome={chromosome} is missing {file}")
        }
        ProbeReason::FingerprintFailed => "could not fingerprint VCF or annotations".into(),
    }
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

impl GenoStoreResult {
    /// Build a `GenotypeResult` view for downstream stages that need a
    /// sample-list + parquet directory pair (phenotype loading, kinship
    /// group loading, known-loci joins). Replaces hand-built clones at
    /// every call site.
    pub fn to_genotype_result(&self) -> super::genotype::GenotypeResult {
        super::genotype::GenotypeResult {
            sample_names: self.sample_names.clone(),
            output_dir: self
                .geno_output_dir
                .clone()
                .unwrap_or_else(|| self.store_dir.clone()),
        }
    }
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
        StoreProbe {
            store_dir: store_dir.to_path_buf(),
            manifest: None,
            miss_reason: None,
        }
    } else {
        probe(store_dir, genotypes, annotations)
    };

    if let Some(manifest) = probe_result.manifest {
        out.status(&format!(
            "Using cached genotype store ({} variants x {} samples)",
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

    let why = match (&probe_result.miss_reason, rebuild_store) {
        (_, true) => "  Rebuild requested by --rebuild-store".to_string(),
        (Some(r), _) => format!("  Cache miss: {}", describe_miss(r)),
        (None, _) => "  Cache miss: probe returned no reason".to_string(),
    };
    out.status(&why);
    out.status("Building genotype store...");

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
    write_atomic(&store_dir.join("manifest.json"), manifest_json.as_bytes())?;

    out.success(&format!(
        "Genotype store: {} variants x {} samples, {} genes, {} carriers",
        total_variants, n_samples, total_genes, total_carriers,
    ));

    Ok(manifest)
}

/// Write `bytes` to `path` atomically: write to `path.tmp`, fsync the file,
/// rename, then fsync the parent directory so the rename survives a crash.
///
/// Parent-dir fsync is best-effort: on filesystems where opening a directory
/// is not allowed (rare on Linux but possible on some network mounts), the
/// rename has already happened so we accept a degraded durability guarantee
/// rather than failing the build.
pub fn write_atomic(path: &Path, bytes: &[u8]) -> Result<(), FavorError> {
    use std::io::Write;
    let parent = path
        .parent()
        .ok_or_else(|| FavorError::Resource(format!("path '{}' has no parent", path.display())))?;
    let tmp = path.with_extension(format!(
        "{}.tmp",
        path.extension().and_then(|s| s.to_str()).unwrap_or("")
    ));
    {
        let mut f = File::create(&tmp)
            .map_err(|e| FavorError::Resource(format!("create {}: {e}", tmp.display())))?;
        f.write_all(bytes)
            .map_err(|e| FavorError::Resource(format!("write {}: {e}", tmp.display())))?;
        f.sync_all()
            .map_err(|e| FavorError::Resource(format!("fsync {}: {e}", tmp.display())))?;
    }
    fs::rename(&tmp, path).map_err(|e| {
        FavorError::Resource(format!(
            "rename {} -> {}: {e}",
            tmp.display(),
            path.display()
        ))
    })?;
    if let Ok(dir) = File::open(parent) {
        let _ = dir.sync_all();
    }
    Ok(())
}

fn now_string() -> String {
    let d = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap_or_default();
    format!("{}", d.as_secs())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn probe_no_manifest_returns_typed_reason() {
        let dir = tempfile::tempdir().unwrap();
        let bogus_vcf = dir.path().join("does-not-exist.vcf.gz");
        let bogus_ann = dir.path().join("does-not-exist.annotated");
        let p = probe(dir.path(), &bogus_vcf, &bogus_ann);
        assert!(p.manifest.is_none());
        assert_eq!(p.miss_reason, Some(ProbeReason::NoManifest));
    }

    #[test]
    fn probe_unreadable_manifest_returns_typed_reason() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("manifest.json"), "{ not valid json").unwrap();
        let bogus_vcf = dir.path().join("does-not-exist.vcf.gz");
        let bogus_ann = dir.path().join("does-not-exist.annotated");
        let p = probe(dir.path(), &bogus_vcf, &bogus_ann);
        assert!(p.manifest.is_none());
        assert_eq!(p.miss_reason, Some(ProbeReason::UnreadableManifest));
    }

    #[test]
    fn write_atomic_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("manifest.json");
        write_atomic(&path, b"hello world").unwrap();
        assert_eq!(std::fs::read(&path).unwrap(), b"hello world");
        // No leftover .tmp.
        let entries: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(Result::ok)
            .collect();
        assert_eq!(entries.len(), 1);
    }

    #[test]
    fn describe_miss_humanises_each_variant() {
        assert!(describe_miss(&ProbeReason::NoManifest).contains("manifest"));
        assert!(describe_miss(&ProbeReason::ContentKeyChanged).contains("fingerprint"));
        assert!(describe_miss(&ProbeReason::SchemaVersionMismatch {
            found: 2,
            expected: 3
        })
        .contains("v2"));
    }
}
