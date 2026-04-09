//! Genotype cohort store: build, cache, probe, and manifest management.

pub mod builder;
pub mod encoding;
pub mod handle;
pub mod membership;
pub mod sparse_g;
pub mod types;
pub mod variants;

pub use handle::{BuildOpts, ChromosomeView, CohortHandle, CohortSources};

pub use crate::store::ids::CohortId;

use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::column;
use crate::engine::DfEngine;
use crate::error::CohortError;
use crate::ingest::{ColumnContract, ColumnRequirement};
use crate::output::Output;
use crate::store::list::VariantSet;

#[derive(Serialize, Deserialize)]
pub struct CohortManifest {
    pub version: u32,
    pub key: String,
    pub n_samples: usize,
    pub n_variants: usize,
    pub chromosomes: Vec<ChromInfo>,
    pub created_at: String,
    pub cohort_version: String,
}

#[derive(Serialize, Deserialize)]
pub struct ChromInfo {
    pub name: String,
    pub n_variants: usize,
}

/// Content-based fingerprint so cache keys survive path renames.
fn file_content_fingerprint(path: &Path) -> Result<Vec<u8>, CohortError> {
    use std::io::{Read as IoRead, Seek, SeekFrom};
    const CHUNK: u64 = 1024 * 1024;

    let mut f = File::open(path)
        .map_err(|e| CohortError::Resource(format!("open {}: {e}", path.display())))?;
    let size = f
        .metadata()
        .map_err(|e| CohortError::Resource(format!("stat {}: {e}", path.display())))?
        .len();

    let mut hasher = Sha256::new();
    hasher.update(size.to_le_bytes());

    let head = CHUNK.min(size) as usize;
    let mut buf = vec![0u8; head];
    f.read_exact(&mut buf)
        .map_err(|e| CohortError::Resource(format!("read {}: {e}", path.display())))?;
    hasher.update(&buf);

    if size > CHUNK {
        let tail_start = size.saturating_sub(CHUNK);
        f.seek(SeekFrom::Start(tail_start))
            .map_err(|e| CohortError::Resource(format!("seek {}: {e}", path.display())))?;
        let tail_len = (size - tail_start) as usize;
        buf.resize(tail_len, 0);
        f.read_exact(&mut buf)
            .map_err(|e| CohortError::Resource(format!("read {}: {e}", path.display())))?;
        hasher.update(&buf);
    }

    Ok(hasher.finalize().to_vec())
}

fn dir_fingerprint(path: &Path) -> Result<Vec<u8>, CohortError> {
    let meta_path = path.join("meta.json");
    if meta_path.exists() {
        return file_content_fingerprint(&meta_path);
    }
    let meta = fs::metadata(path)
        .map_err(|e| CohortError::Resource(format!("stat {}: {e}", path.display())))?;
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

pub fn compute_key(vcf_path: &Path, annotations_path: &Path) -> Result<String, CohortError> {
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
    pub manifest: Option<CohortManifest>,
    pub miss_reason: Option<ProbeReason>,
}

pub fn probe(store_dir: &Path, vcf_path: &Path, annotations_path: &Path) -> StoreProbe {
    let store_dir = store_dir.to_path_buf();
    let miss = |reason: ProbeReason| StoreProbe {
        store_dir: store_dir.clone(),
        manifest: None,
        miss_reason: Some(reason),
    };

    if let Ok(staging) = staging_path(&store_dir) {
        let _ = finish_interrupted_swap(&staging, &store_dir);
    }

    let Ok(s) = fs::read_to_string(store_dir.join("manifest.json")) else {
        return miss(ProbeReason::NoManifest);
    };
    let Ok(manifest) = serde_json::from_str::<CohortManifest>(&s) else {
        return miss(ProbeReason::UnreadableManifest);
    };
    let Ok(key) = compute_key(vcf_path, annotations_path) else {
        return miss(ProbeReason::FingerprintFailed);
    };

    if manifest.version != 4 {
        return miss(ProbeReason::SchemaVersionMismatch {
            found: manifest.version,
            expected: 4,
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
    // For Existing cohorts geno.output_dir is the cohort dir itself; for Fresh
    // builds it is the staging dir from extract_genotypes. store_dir always
    // points at the cohort dir under `.cohort/cohorts/<id>/`.
    pub geno: crate::staar::genotype::GenotypeResult,
    pub manifest: CohortManifest,
    pub store_dir: PathBuf,
}

pub(crate) fn read_sample_names_at(store_dir: &Path) -> Result<Vec<String>, CohortError> {
    let path = store_dir.join("samples.txt");
    let content = std::fs::read_to_string(&path)
        .map_err(|e| CohortError::Resource(format!("Read {}: {e}", path.display())))?;
    Ok(content.lines().map(|s| s.to_string()).collect())
}

/// Register the annotation + genotype tables on `df` and execute the join,
/// leaving `_rare_all` registered for downstream queries. Caller owns the
/// engine — typically `engine.df()` from the composition root. Drops any
/// pre-existing copies of the temp tables so a second build inside the
/// same process is idempotent.
pub fn run_annotation_join(
    df: &DfEngine,
    annotations: &Path,
    geno: &crate::staar::genotype::GenotypeResult,
) -> Result<(), CohortError> {
    for tbl in ["_rare_all", "_ann_check", "_genotypes", "_annotations"] {
        let _ = df.execute(&format!("DROP TABLE IF EXISTS {tbl}"));
    }

    let ann_vs = VariantSet::open(annotations)?;

    // Check the table schema directly; LIMIT 0 may return no batches in DataFusion.
    let ann_cols: Vec<String> = {
        df.register_parquet_dir("_ann_check", ann_vs.root())?;
        let cols = df.table_columns("_ann_check")?;
        let _ = df.execute("DROP TABLE IF EXISTS _ann_check");
        cols
    };

    let contract = ColumnContract {
        command: "staar (annotation join)",
        required: STAAR_ANNOTATION_COLUMNS,
    };
    let missing = contract.check(&ann_cols);
    if !missing.is_empty() {
        let tier_hint = if matches!(ann_vs.tier(), Some(crate::config::Tier::Base)) {
            " Your data was annotated with base tier. Re-run: `cohort annotate --full`."
        } else {
            " Re-run: `cohort annotate --full`."
        };
        return Err(CohortError::DataMissing(format!(
            "Missing annotation columns in {}:\n{}\n\
             STAAR requires FAVOR full-tier annotations.{}",
            annotations.display(),
            ColumnContract::format_missing(&missing),
            tier_hint,
        )));
    }

    let geno_dir = &geno.output_dir;
    df.register_parquet_sorted("_genotypes", geno_dir, &["position", "ref", "alt"])?;
    df.register_parquet_sorted(
        "_annotations",
        ann_vs.root(),
        &["position", "ref_vcf", "alt_vcf"],
    )?;

    df.execute(&column::annotation_join_sql())?;

    Ok(())
}

pub fn build(
    engine: &DfEngine,
    geno_result: &crate::staar::genotype::GenotypeResult,
    store_dir: &Path,
    vcf_path: &Path,
    annotations_path: &Path,
    out: &dyn Output,
) -> Result<CohortManifest, CohortError> {
    let staging = staging_path(store_dir)?;
    finish_interrupted_swap(&staging, store_dir)?;
    fs::create_dir_all(&staging).map_err(|e| {
        CohortError::Resource(format!("Cannot create '{}': {e}", staging.display()))
    })?;

    // Hold a writer lock for the duration of the build so two parallel
    // `cohort ingest` runs cannot race the staging swap. Lock file lives
    // next to staging so it survives the rename of `store_dir`.
    let _lock = BuildLock::acquire(&staging.with_extension("lock"))?;

    let key = compute_key(vcf_path, annotations_path)?;
    let n_samples = geno_result.sample_names.len();

    let samples_path = staging.join("samples.txt");
    fs::write(&samples_path, geno_result.sample_names.join("\n")).map_err(|e| {
        CohortError::Resource(format!("Cannot write '{}': {e}", samples_path.display()))
    })?;

    let chroms = engine.query_strings(&column::distinct_chroms_sql())?;

    let mut chrom_infos: Vec<ChromInfo> = Vec::new();
    let mut total_variants = 0usize;
    let mut total_genes = 0usize;
    let mut total_carriers = 0u64;

    for chrom in &chroms {
        out.status(&format!("  Building sparse genotype store chr{chrom}..."));

        let chrom_dir = staging.join(format!("chromosome={chrom}"));
        fs::create_dir_all(&chrom_dir).map_err(|e| {
            CohortError::Resource(format!("Cannot create '{}': {e}", chrom_dir.display()))
        })?;

        let geno_path = format!(
            "{}/chromosome={chrom}/data.parquet",
            geno_result.output_dir.display()
        );
        let stats = builder::build_chromosome(
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

    let manifest = CohortManifest {
        version: 4,
        key,
        n_samples,
        n_variants: total_variants,
        chromosomes: chrom_infos,
        created_at: now_string(),
        cohort_version: env!("CARGO_PKG_VERSION").to_string(),
    };

    let manifest_json = serde_json::to_string_pretty(&manifest)
        .map_err(|e| CohortError::Resource(format!("Manifest serialize: {e}")))?;
    write_atomic(&staging.join("manifest.json"), manifest_json.as_bytes())?;

    atomic_dir_swap(&staging, store_dir)?;

    out.success(&format!(
        "Genotype store: {} variants x {} samples, {} genes, {} carriers",
        total_variants, n_samples, total_genes, total_carriers,
    ));

    Ok(manifest)
}

/// `{store_dir}.staging` — sibling so the recovery code finds it without
/// an extra config knob.
pub fn staging_path(store_dir: &Path) -> Result<PathBuf, CohortError> {
    let parent = store_dir.parent().ok_or_else(|| {
        CohortError::Resource(format!(
            "store dir '{}' has no parent",
            store_dir.display()
        ))
    })?;
    let name = store_dir.file_name().ok_or_else(|| {
        CohortError::Resource(format!(
            "store dir '{}' has no file name",
            store_dir.display()
        ))
    })?;
    let mut staging_name = name.to_os_string();
    staging_name.push(".staging");
    Ok(parent.join(staging_name))
}

/// Replace `final_path` with `staging`. Not strictly atomic — there's a
/// window between `remove_dir_all` and `rename` where neither exists; a
/// crash there is recovered by `finish_interrupted_swap` on the next run.
pub fn atomic_dir_swap(staging: &Path, final_path: &Path) -> Result<(), CohortError> {
    if final_path.exists() {
        fs::remove_dir_all(final_path).map_err(|e| {
            CohortError::Resource(format!(
                "remove old store {}: {e}",
                final_path.display()
            ))
        })?;
    }
    fs::rename(staging, final_path).map_err(|e| {
        CohortError::Resource(format!(
            "rename {} -> {}: {e}",
            staging.display(),
            final_path.display()
        ))
    })?;
    fsync_parent(final_path);
    Ok(())
}

/// Called at the top of `build`. Three cases:
/// - staging has a `manifest.json` and final doesn't: a prior run crashed
///   inside `atomic_dir_swap`, between the `remove_dir_all` and the
///   `rename`. Finish the rename so we don't throw away a good build.
/// - staging exists with no manifest: prior build died mid-way, the
///   contents are garbage. Wipe it.
/// - staging absent: nothing to do.
pub fn finish_interrupted_swap(staging: &Path, final_path: &Path) -> Result<(), CohortError> {
    if !staging.exists() {
        return Ok(());
    }
    let staging_complete = staging.join("manifest.json").exists();
    if staging_complete && !final_path.exists() {
        fs::rename(staging, final_path).map_err(|e| {
            CohortError::Resource(format!(
                "finish swap {} -> {}: {e}",
                staging.display(),
                final_path.display()
            ))
        })?;
        fsync_parent(final_path);
        return Ok(());
    }
    fs::remove_dir_all(staging).map_err(|e| {
        CohortError::Resource(format!(
            "remove stale staging {}: {e}",
            staging.display()
        ))
    })
}

use crate::store::manifest::{fsync_parent, write_atomic};

/// Exclusive build lock, file-based. Held for the duration of a cohort
/// build so concurrent runs serialize on the staging swap. Drops the lock
/// (and removes the file) on `Drop` so a panic mid-build still releases.
struct BuildLock {
    #[cfg(unix)]
    file: std::fs::File,
    path: PathBuf,
}

impl BuildLock {
    fn acquire(path: &Path) -> Result<Self, CohortError> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent).map_err(|e| {
                CohortError::Resource(format!("create {}: {e}", parent.display()))
            })?;
        }
        let file = std::fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(false)
            .open(path)
            .map_err(|e| CohortError::Resource(format!("open lock {}: {e}", path.display())))?;
        #[cfg(unix)]
        {
            use std::os::unix::io::AsRawFd;
            let rc = unsafe { libc::flock(file.as_raw_fd(), libc::LOCK_EX | libc::LOCK_NB) };
            if rc != 0 {
                return Err(CohortError::Resource(format!(
                    "another cohort build holds the lock at {}; wait for it to finish or remove the file",
                    path.display()
                )));
            }
        }
        Ok(Self {
            #[cfg(unix)]
            file,
            path: path.to_path_buf(),
        })
    }
}

impl Drop for BuildLock {
    fn drop(&mut self) {
        #[cfg(unix)]
        {
            use std::os::unix::io::AsRawFd;
            let _ = unsafe { libc::flock(self.file.as_raw_fd(), libc::LOCK_UN) };
        }
        let _ = fs::remove_file(&self.path);
    }
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
    fn describe_miss_humanises_each_variant() {
        assert!(describe_miss(&ProbeReason::NoManifest).contains("manifest"));
        assert!(describe_miss(&ProbeReason::ContentKeyChanged).contains("fingerprint"));
        assert!(describe_miss(&ProbeReason::SchemaVersionMismatch {
            found: 3,
            expected: 4
        })
        .contains("v3"));
    }

    #[test]
    fn staging_path_is_sibling_with_staging_suffix() {
        let store = PathBuf::from("/tmp/foo/store");
        let staging = staging_path(&store).unwrap();
        assert_eq!(staging, PathBuf::from("/tmp/foo/store.staging"));
    }

    #[test]
    fn atomic_dir_swap_round_trip() {
        let dir = tempfile::tempdir().unwrap();
        let staging = dir.path().join("store.staging");
        let final_path = dir.path().join("store");

        std::fs::create_dir_all(&staging).unwrap();
        std::fs::write(staging.join("manifest.json"), b"{}").unwrap();

        atomic_dir_swap(&staging, &final_path).unwrap();

        assert!(final_path.exists());
        assert!(!staging.exists());
        assert!(final_path.join("manifest.json").exists());
    }

    #[test]
    fn atomic_dir_swap_replaces_existing_final() {
        let dir = tempfile::tempdir().unwrap();
        let staging = dir.path().join("store.staging");
        let final_path = dir.path().join("store");

        std::fs::create_dir_all(&final_path).unwrap();
        std::fs::write(final_path.join("old.txt"), b"old").unwrap();
        std::fs::create_dir_all(&staging).unwrap();
        std::fs::write(staging.join("new.txt"), b"new").unwrap();

        atomic_dir_swap(&staging, &final_path).unwrap();

        assert!(final_path.join("new.txt").exists());
        assert!(!final_path.join("old.txt").exists());
        assert!(!staging.exists());
    }

    #[test]
    fn finish_interrupted_swap_completes_when_manifest_present() {
        // Prior run crashed between remove_dir_all(final) and rename(staging,
        // final). Recovery: finish the swap.
        let dir = tempfile::tempdir().unwrap();
        let staging = dir.path().join("store.staging");
        let final_path = dir.path().join("store");

        std::fs::create_dir_all(&staging).unwrap();
        std::fs::write(staging.join("manifest.json"), b"{}").unwrap();

        finish_interrupted_swap(&staging, &final_path).unwrap();

        assert!(final_path.exists());
        assert!(final_path.join("manifest.json").exists());
        assert!(!staging.exists());
    }

    #[test]
    fn finish_interrupted_swap_cleans_unfinished_staging() {
        // Prior run crashed mid-build (no manifest in staging). Recovery:
        // wipe the staging dir, leave the existing final alone.
        let dir = tempfile::tempdir().unwrap();
        let staging = dir.path().join("store.staging");
        let final_path = dir.path().join("store");

        std::fs::create_dir_all(&final_path).unwrap();
        std::fs::write(final_path.join("manifest.json"), b"{}").unwrap();
        std::fs::create_dir_all(&staging).unwrap();
        std::fs::write(staging.join("garbage.txt"), b"x").unwrap();

        finish_interrupted_swap(&staging, &final_path).unwrap();

        assert!(!staging.exists());
        assert!(final_path.join("manifest.json").exists());
    }

    #[test]
    fn finish_interrupted_swap_no_op_when_staging_absent() {
        let dir = tempfile::tempdir().unwrap();
        let staging = dir.path().join("store.staging");
        let final_path = dir.path().join("store");
        finish_interrupted_swap(&staging, &final_path).unwrap();
        assert!(!staging.exists());
        assert!(!final_path.exists());
    }
}
