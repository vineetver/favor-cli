use std::collections::HashMap;
use std::io::{IsTerminal, Read, Write};
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Sha256, Digest};

use crate::cli::DataAction;
use crate::config::{Config, Tier};
use crate::error::FavorError;
use crate::output::Output;
use super::Pack;

pub const REMOTE_BASE_URL: &str =
    "https://minio-s3-favor-4ee4be.apps.shift.nerc.mghpcc.org/favor-hg38";

const CHROMOSOMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "X", "Y",
];

const MAX_RETRIES: u32 = 3;

/// The manifest schema version this CLI understands.
/// Manifests with a higher version are rejected.
pub const MANIFEST_SCHEMA_VERSION: u32 = 2;

/// Map Tier to MinIO path. Sum type — no invalid input possible.
fn remote_tier_path(tier: Tier) -> &'static str {
    match tier {
        Tier::Full => "stage16/v2",
        Tier::Base => "base",
    }
}

pub fn human_size(bytes: u64) -> String {
    const UNITS: &[&str] = &["B", "KiB", "MiB", "GiB", "TiB"];
    let mut size = bytes as f64;
    for unit in UNITS {
        if size < 1024.0 {
            return format!("{size:.1} {unit}");
        }
        size /= 1024.0;
    }
    format!("{size:.1} PiB")
}

/// A single file entry in a pack manifest.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FileEntry {
    pub size: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sha256: Option<String>,
}

/// Pack manifest (v2 canonical format).
/// The CLI also accepts v1 manifests (files as bare u64 sizes) via `parse_pack_manifest`.
#[derive(Debug, Serialize, Deserialize)]
pub struct PackManifest {
    pub schema_version: u32,
    pub pack_id: String,
    pub version: String,
    #[serde(default)]
    pub base_dir: String,
    pub total_size: u64,
    pub files: HashMap<String, FileEntry>,
}

#[derive(Deserialize)]
struct PackManifestV1 {
    pack_id: String,
    version: String,
    total_size: u64,
    files: HashMap<String, u64>,
}

/// Parse a manifest body, handling both v1 (files as u64) and v2 (files as FileEntry).
fn parse_pack_manifest(body: &str) -> Result<PackManifest, FavorError> {
    let raw: Value = serde_json::from_str(body).map_err(|e| {
        FavorError::Resource(format!("Invalid manifest JSON: {e}"))
    })?;

    let schema_version = raw.get("schema_version")
        .and_then(|v| v.as_u64())
        .unwrap_or(1) as u32;

    if schema_version > MANIFEST_SCHEMA_VERSION {
        return Err(FavorError::Resource(format!(
            "Manifest schema version {schema_version} requires a newer CLI (this CLI supports up to v{MANIFEST_SCHEMA_VERSION}). \
             Update favor or download the latest release."
        )));
    }

    if schema_version >= 2 {
        // v2: files are {size, sha256}
        let manifest: PackManifest = serde_json::from_str(body).map_err(|e| {
            FavorError::Resource(format!("Invalid v2 manifest: {e}"))
        })?;
        return Ok(manifest);
    }

    // v1: files are bare u64 sizes
    let v1: PackManifestV1 = serde_json::from_str(body).map_err(|e| {
        FavorError::Resource(format!("Invalid v1 manifest: {e}"))
    })?;
    Ok(PackManifest {
        schema_version: 1,
        pack_id: v1.pack_id,
        version: v1.version,
        base_dir: String::new(), // v1 didn't have this
        total_size: v1.total_size,
        files: v1.files.into_iter()
            .map(|(k, size)| (k, FileEntry { size, sha256: None }))
            .collect(),
    })
}

#[derive(Debug, Serialize, Deserialize)]
struct AnnotationManifest {
    tier: String,
    version: String,
    total_size: u64,
    chromosomes: HashMap<String, ChromEntry>,
}

#[derive(Debug, Serialize, Deserialize)]
struct ChromEntry {
    size: u64,
}

#[derive(Debug, Serialize)]
#[serde(tag = "state")]
enum ChromState {
    Valid { chrom: String, path: PathBuf, size: u64 },
    SizeMismatch { chrom: String, path: PathBuf, local_size: u64, expected_size: u64 },
    Partial { chrom: String, tmp_path: PathBuf, downloaded: u64, expected_size: u64 },
    Missing { chrom: String, expected_size: u64 },
}

impl ChromState {
    fn chrom(&self) -> &str {
        match self {
            Self::Valid { chrom, .. }
            | Self::SizeMismatch { chrom, .. }
            | Self::Partial { chrom, .. }
            | Self::Missing { chrom, .. } => chrom,
        }
    }

    fn needs_download(&self) -> bool {
        !matches!(self, Self::Valid { .. })
    }
}

struct DataState {
    tier: Tier,
    data_dir: PathBuf,
    manifest: AnnotationManifest,
    chromosomes: Vec<ChromState>,
}

impl DataState {
    fn valid_count(&self) -> usize {
        self.chromosomes.iter().filter(|c| matches!(c, ChromState::Valid { .. })).count()
    }

    fn is_complete(&self) -> bool {
        self.valid_count() == 24
    }

    fn needs_download(&self) -> Vec<&ChromState> {
        self.chromosomes.iter().filter(|c| c.needs_download()).collect()
    }

    fn total_to_download(&self) -> u64 {
        self.chromosomes.iter().filter_map(|c| match c {
            ChromState::Missing { expected_size, .. } => Some(*expected_size),
            ChromState::SizeMismatch { expected_size, .. } => Some(*expected_size),
            ChromState::Partial { expected_size, downloaded, .. } => Some(expected_size - downloaded),
            ChromState::Valid { .. } => None,
        }).sum()
    }
}

#[derive(Debug)]
enum PackFileState {
    Valid { rel_path: String, size: u64 },
    SizeMismatch { rel_path: String, expected_size: u64 },
    Missing { rel_path: String, expected_size: u64 },
}

impl PackFileState {
    fn needs_download(&self) -> bool {
        !matches!(self, Self::Valid { .. })
    }
    fn rel_path(&self) -> &str {
        match self {
            Self::Valid { rel_path, .. }
            | Self::SizeMismatch { rel_path, .. }
            | Self::Missing { rel_path, .. } => rel_path,
        }
    }
    fn expected_size(&self) -> u64 {
        match self {
            Self::Valid { size, .. } => *size,
            Self::SizeMismatch { expected_size, .. }
            | Self::Missing { expected_size, .. } => *expected_size,
        }
    }
}

fn fetch_annotation_manifest(tier: Tier) -> Result<AnnotationManifest, FavorError> {
    let remote_path = remote_tier_path(tier);
    let url = format!("{REMOTE_BASE_URL}/{remote_path}/manifest.json");

    let response = ureq::get(&url).call().map_err(|e| {
        FavorError::Resource(format!("Failed to fetch manifest: {e}"))
    })?;

    let body = response.into_body().read_to_string().map_err(|e| {
        FavorError::Resource(format!("Failed to read manifest: {e}"))
    })?;
    let manifest: AnnotationManifest = serde_json::from_str(&body).map_err(|e| {
        FavorError::Resource(format!("Invalid manifest JSON: {e}"))
    })?;

    Ok(manifest)
}

fn assess_local_annotations(data_dir: &Path, manifest: &AnnotationManifest) -> Vec<ChromState> {
    CHROMOSOMES.iter().map(|chrom| {
        let expected = manifest.chromosomes.get(*chrom);
        let expected_size = expected.map(|e| e.size).unwrap_or(0);
        let parquet = data_dir.join(format!("chromosome={chrom}")).join("sorted.parquet");
        let tmp = parquet.with_extension("parquet.tmp");

        if parquet.exists() {
            let local_size = parquet.metadata().map(|m| m.len()).unwrap_or(0);
            if local_size == expected_size {
                ChromState::Valid { chrom: chrom.to_string(), path: parquet, size: local_size }
            } else {
                ChromState::SizeMismatch {
                    chrom: chrom.to_string(), path: parquet,
                    local_size, expected_size,
                }
            }
        } else if tmp.exists() {
            let downloaded = tmp.metadata().map(|m| m.len()).unwrap_or(0);
            ChromState::Partial {
                chrom: chrom.to_string(), tmp_path: tmp,
                downloaded, expected_size,
            }
        } else {
            ChromState::Missing { chrom: chrom.to_string(), expected_size }
        }
    }).collect()
}

fn build_data_state(config: &Config) -> Result<DataState, FavorError> {
    let tier = config.data.tier;
    let data_dir = config.annotations_dir();

    let manifest = fetch_annotation_manifest(tier)?;
    let chromosomes = assess_local_annotations(&data_dir, &manifest);

    Ok(DataState { tier, data_dir, manifest, chromosomes })
}

pub fn fetch_pack_manifest(pack_id: &str) -> Result<PackManifest, FavorError> {
    let pack = super::Pack::find(pack_id);
    let base = pack.map(|p| p.base_url(REMOTE_BASE_URL)).unwrap_or(REMOTE_BASE_URL);
    let url = format!("{base}/packs/{pack_id}/manifest.json");

    let response = ureq::get(&url).call().map_err(|e| {
        FavorError::Resource(format!("Failed to fetch pack manifest for '{pack_id}': {e}"))
    })?;

    let body = response.into_body().read_to_string().map_err(|e| {
        FavorError::Resource(format!("Failed to read pack manifest: {e}"))
    })?;

    let manifest = parse_pack_manifest(&body)?;

    // Cross-check: manifest pack_id must match requested pack_id
    if manifest.pack_id != pack_id {
        return Err(FavorError::Resource(format!(
            "Manifest pack_id '{}' does not match requested '{pack_id}'",
            manifest.pack_id,
        )));
    }

    // Cross-check: if v2 manifest has base_dir, it must match the compiled-in registry
    if manifest.schema_version >= 2 {
        if let Some(pack) = Pack::find(pack_id) {
            if manifest.base_dir != pack.base_dir {
                return Err(FavorError::Resource(format!(
                    "Manifest base_dir '{}' does not match CLI registry '{}' for pack '{}'. \
                     This indicates a publish/registry mismatch — re-run `favor data publish`.",
                    manifest.base_dir, pack.base_dir, pack_id,
                )));
            }
        }
    }

    Ok(manifest)
}

fn assess_pack_state(local_base: &Path, manifest: &PackManifest) -> Vec<PackFileState> {
    let mut states: Vec<PackFileState> = manifest.files.iter().map(|(rel_path, entry)| {
        let local_path = local_base.join(rel_path);
        if local_path.exists() {
            let local_size = local_path.metadata().map(|m| m.len()).unwrap_or(0);
            if local_size == entry.size {
                PackFileState::Valid { rel_path: rel_path.clone(), size: local_size }
            } else {
                PackFileState::SizeMismatch { rel_path: rel_path.clone(), expected_size: entry.size }
            }
        } else {
            PackFileState::Missing { rel_path: rel_path.clone(), expected_size: entry.size }
        }
    }).collect();
    states.sort_by(|a, b| a.rel_path().cmp(b.rel_path()));
    states
}

pub fn run(action: DataAction, output: &dyn Output) -> Result<(), FavorError> {
    match action {
        DataAction::Status => status(output),
        DataAction::Pull { full, dry_run, yes, pack } => {
            if let Some(pack_id) = pack {
                if pack_id == "all" {
                    pull_all_packs(dry_run, yes, output)
                } else {
                    pull_pack(&pack_id, dry_run, yes, output)
                }
            } else {
                pull(full, dry_run, yes, output)
            }
        }
        DataAction::Verify { pack, checksums } => verify(pack, checksums, output),
        DataAction::Publish { source_root, pack, dry_run } => {
            super::publish::publish(source_root, pack, dry_run, output)
        }
    }
}

fn status(output: &dyn Output) -> Result<(), FavorError> {
    let config = Config::load_configured()?;

    let state = build_data_state(&config)?;

    let rows: Vec<Vec<String>> = state.chromosomes.iter().map(|c| match c {
        ChromState::Valid { chrom, size, .. } => vec![
            format!("chr{chrom}"), "valid".into(), human_size(*size), String::new(),
        ],
        ChromState::SizeMismatch { chrom, local_size, expected_size, .. } => vec![
            format!("chr{chrom}"), "MISMATCH".into(),
            format!("{} / {}", human_size(*local_size), human_size(*expected_size)),
            "will re-download".into(),
        ],
        ChromState::Partial { chrom, downloaded, expected_size, .. } => vec![
            format!("chr{chrom}"), "partial".into(),
            format!("{} / {}", human_size(*downloaded), human_size(*expected_size)),
            "will resume".into(),
        ],
        ChromState::Missing { chrom, expected_size } => vec![
            format!("chr{chrom}"), "missing".into(), human_size(*expected_size), "will download".into(),
        ],
    }).collect();

    output.table(&["Chrom", "State", "Size", "Action"], &rows);

    if !config.data.packs.is_empty() {
        output.status(&format!("Installed packs: {}", config.data.packs.join(", ")));
    }

    if state.is_complete() {
        output.success(&format!(
            "favor-{}: {}/24 chromosomes valid ({})",
            state.tier, state.valid_count(), human_size(state.manifest.total_size),
        ));
    } else {
        let remaining = state.total_to_download();
        output.status(&format!(
            "favor-{}: {}/24 valid, {} to download",
            state.tier, state.valid_count(), human_size(remaining),
        ));
    }

    output.result_json(&json!({
        "tier": state.tier.as_str(),
        "valid": state.valid_count(),
        "total": 24,
        "complete": state.is_complete(),
        "to_download_bytes": state.total_to_download(),
        "packs": config.data.packs,
    }));

    Ok(())
}

fn pull(full: bool, dry_run: bool, yes: bool, output: &dyn Output) -> Result<(), FavorError> {
    let mut config = Config::load_configured()?;

    if full && config.data.tier != Tier::Full {
        config.data.tier = Tier::Full;
        config.save()?;
    }

    let state = build_data_state(&config)?;

    if state.is_complete() {
        output.success(&format!(
            "favor-{}: all 24 chromosomes valid ({})",
            state.tier, human_size(state.manifest.total_size),
        ));
        return Ok(());
    }

    let to_download = state.needs_download();
    let total_bytes = state.total_to_download();

    output.status(&format!(
        "favor-{}: {}/24 valid, downloading {} ({})",
        state.tier,
        state.valid_count(),
        to_download.len(),
        human_size(total_bytes),
    ));

    if dry_run {
        output.result_json(&json!({
            "tier": state.tier.as_str(),
            "valid": state.valid_count(),
            "to_download": to_download.iter().map(|c| c.chrom()).collect::<Vec<_>>(),
            "to_download_bytes": total_bytes,
        }));
        return Ok(());
    }

    if !yes && total_bytes > 0 && std::io::stdin().is_terminal()
        && !confirm_download(total_bytes)
    {
        output.warn("Download cancelled.");
        return Ok(());
    }

    let remote_path = remote_tier_path(state.tier);
    let show_progress = std::io::stderr().is_terminal();

    for (i, chrom_state) in to_download.iter().enumerate() {
        let chrom = chrom_state.chrom();
        let target_dir = state.data_dir.join(format!("chromosome={chrom}"));
        let target = target_dir.join("sorted.parquet");
        let tmp = target.with_extension("parquet.tmp");
        let url = format!("{REMOTE_BASE_URL}/{remote_path}/chromosome={chrom}/sorted.parquet");

        std::fs::create_dir_all(&target_dir)?;

        if matches!(chrom_state, ChromState::SizeMismatch { .. }) {
            let _ = std::fs::remove_file(&target);
        }

        let resume_from = if tmp.exists() {
            tmp.metadata().map(|m| m.len()).unwrap_or(0)
        } else {
            0
        };

        let expected_size = match chrom_state {
            ChromState::Missing { expected_size, .. }
            | ChromState::SizeMismatch { expected_size, .. }
            | ChromState::Partial { expected_size, .. } => *expected_size,
            _ => 0,
        };

        let label = format!("{chrom} ({}/{})", i + 1, to_download.len());

        let result = download_with_retry(&url, &tmp, resume_from, expected_size, &label, show_progress, output);

        match result {
            Ok(()) => {
                let actual = tmp.metadata().map(|m| m.len()).unwrap_or(0);
                if actual != expected_size && expected_size > 0 {
                    output.warn(&format!(
                        "chr{chrom}: size mismatch ({} vs expected {}), keeping partial",
                        human_size(actual), human_size(expected_size),
                    ));
                } else {
                    std::fs::rename(&tmp, &target)?;
                    output.success(&format!("chr{chrom} ({}/{})", i + 1, to_download.len()));
                }
            }
            Err(e) => {
                output.warn(&format!("chr{chrom} failed: {e} (partial saved, will resume next run)"));
            }
        }
    }

    let final_state = build_data_state(&config)?;
    if final_state.is_complete() {
        output.success(&format!(
            "favor-{} complete: 24/24 chromosomes ({})",
            state.tier, human_size(state.manifest.total_size),
        ));
    } else {
        let remaining = final_state.needs_download().len();
        output.warn(&format!(
            "{remaining} chromosomes still incomplete. Run `favor data pull` to retry.",
        ));
    }

    Ok(())
}

fn ensure_pack_registered(config: &mut Config, pack: &Pack) -> Result<(), FavorError> {
    if !pack.always_installed && !config.data.packs.contains(&pack.id.to_string()) {
        config.data.packs.push(pack.id.to_string());
        config.save()?;
    }
    Ok(())
}

pub fn pull_pack(pack_id: &str, dry_run: bool, yes: bool, output: &dyn Output) -> Result<(), FavorError> {
    let pack = Pack::find(pack_id).ok_or_else(|| {
        let available: Vec<_> = Pack::all().iter().map(|p| p.id).collect();
        FavorError::Input(format!(
            "Unknown pack '{}'. Available: {}",
            pack_id,
            available.join(", "),
        ))
    })?;

    let mut config = Config::load_configured()?;

    let local_base = pack.local_dir(&config.root_dir());

    output.status(&format!("Pack '{}' — {} ({})", pack.name, pack.source, pack.size_human));

    let manifest = fetch_pack_manifest(pack_id)?;
    let states = assess_pack_state(&local_base, &manifest);

    let to_download: Vec<_> = states.iter().filter(|s| s.needs_download()).collect();
    let valid = states.len() - to_download.len();

    if to_download.is_empty() {
        output.success(&format!("Pack '{}': {}/{} files valid ({})",
            pack_id, valid, states.len(), human_size(manifest.total_size)));
        ensure_pack_registered(&mut config, pack)?;
        return Ok(());
    }

    let total_bytes: u64 = to_download.iter().map(|s| s.expected_size()).sum();

    output.status(&format!("{}/{} valid, {} files to download ({})",
        valid, states.len(), to_download.len(), human_size(total_bytes)));

    if dry_run {
        output.result_json(&json!({
            "pack": pack_id,
            "valid": valid,
            "total_files": states.len(),
            "to_download": to_download.len(),
            "to_download_bytes": total_bytes,
        }));
        return Ok(());
    }

    if !yes && total_bytes > 0 && std::io::stdin().is_terminal()
        && !confirm_download(total_bytes)
    {
        output.warn("Download cancelled.");
        return Ok(());
    }

    let show_progress = std::io::stderr().is_terminal();

    for (i, file_state) in to_download.iter().enumerate() {
        let rel_path = file_state.rel_path();
        let local_path = local_base.join(rel_path);
        let tmp_path = local_path.with_extension("parquet.tmp");
        let pack = super::Pack::find(pack_id);
        let base = pack.map(|p| p.base_url(REMOTE_BASE_URL)).unwrap_or(REMOTE_BASE_URL);
        let url = format!("{base}/packs/{pack_id}/{rel_path}");

        if let Some(parent) = local_path.parent() {
            std::fs::create_dir_all(parent)?;
        }

        if matches!(file_state, PackFileState::SizeMismatch { .. }) {
            let _ = std::fs::remove_file(&local_path);
        }

        let resume_from = if tmp_path.exists() {
            tmp_path.metadata().map(|m| m.len()).unwrap_or(0)
        } else {
            0
        };

        let expected_size = file_state.expected_size();
        let table_name = rel_path.split('/').next().unwrap_or(rel_path);
        let label = format!("{table_name} ({}/{})", i + 1, to_download.len());

        match download_with_retry(&url, &tmp_path, resume_from, expected_size, &label, show_progress, output) {
            Ok(()) => {
                let actual = tmp_path.metadata().map(|m| m.len()).unwrap_or(0);
                if actual != expected_size && expected_size > 0 {
                    output.warn(&format!(
                        "{rel_path}: size mismatch ({} vs expected {}), keeping partial",
                        human_size(actual), human_size(expected_size),
                    ));
                } else {
                    std::fs::rename(&tmp_path, &local_path)?;
                }
            }
            Err(e) => {
                output.warn(&format!("{rel_path} failed: {e} (partial saved)"));
            }
        }
    }

    ensure_pack_registered(&mut config, pack)?;

    let final_states = assess_pack_state(&local_base, &manifest);
    let final_valid = final_states.iter().filter(|s| matches!(s, PackFileState::Valid { .. })).count();

    if final_valid == final_states.len() {
        output.success(&format!("Pack '{}' complete: {}/{} files ({})",
            pack_id, final_valid, final_states.len(), human_size(manifest.total_size)));
    } else {
        let remaining = final_states.len() - final_valid;
        output.warn(&format!("{remaining} files incomplete. Run `favor data pull --pack {pack_id}` to retry."));
    }

    Ok(())
}

fn pull_all_packs(dry_run: bool, yes: bool, output: &dyn Output) -> Result<(), FavorError> {
    let config = Config::load_configured()?;

    if config.data.packs.is_empty() {
        output.warn("No packs configured. Use `favor data pull --pack <name>` or run `favor setup`.");
        return Ok(());
    }

    let pack_ids = config.data.packs.clone();
    output.status(&format!("Pulling {} configured packs: {}", pack_ids.len(), pack_ids.join(", ")));

    for pack_id in &pack_ids {
        if let Err(e) = pull_pack(pack_id, dry_run, yes, output) {
            output.warn(&format!("Pack '{pack_id}' failed: {e}"));
        }
    }

    Ok(())
}

fn verify(pack_filter: Option<String>, checksums: bool, output: &dyn Output) -> Result<(), FavorError> {
    let config = Config::load_configured()?;

    let packs_to_verify: Vec<&Pack> = if let Some(ref id) = pack_filter {
        let pack = Pack::find(id).ok_or_else(|| {
            FavorError::Input(format!("Unknown pack '{id}'"))
        })?;
        vec![pack]
    } else {
        // Verify all installed packs
        let mut packs = Vec::new();
        for pack in Pack::required() {
            packs.push(pack);
        }
        for pack_id in &config.data.packs {
            if let Some(pack) = Pack::find(pack_id) {
                if !pack.always_installed {
                    packs.push(pack);
                }
            }
        }
        packs
    };

    let mut all_ok = true;

    for pack in &packs_to_verify {
        output.status(&format!("Verifying pack '{}'...", pack.id));

        let manifest = match fetch_pack_manifest(pack.id) {
            Ok(m) => m,
            Err(e) => {
                output.warn(&format!("  {} — failed to fetch manifest: {e}", pack.id));
                all_ok = false;
                continue;
            }
        };

        output.status(&format!("  schema v{}, {} files", manifest.schema_version, manifest.files.len()));

        // base_dir cross-check (v2 manifests validated in fetch_pack_manifest already)
        if manifest.schema_version >= 2 {
            output.status(&format!("  base_dir: '{}' (matches registry)", manifest.base_dir));
        } else {
            output.status("  base_dir: (v1 manifest, no cross-check)");
        }

        let local_base = pack.local_dir(&config.root_dir());
        let states = assess_pack_state(&local_base, &manifest);

        let valid = states.iter().filter(|s| matches!(s, PackFileState::Valid { .. })).count();
        let mismatched: Vec<_> = states.iter()
            .filter(|s| matches!(s, PackFileState::SizeMismatch { .. }))
            .collect();
        let missing: Vec<_> = states.iter()
            .filter(|s| matches!(s, PackFileState::Missing { .. }))
            .collect();

        if !mismatched.is_empty() {
            all_ok = false;
            for s in &mismatched {
                output.warn(&format!("  SIZE MISMATCH: {}", s.rel_path()));
            }
        }
        if !missing.is_empty() {
            all_ok = false;
            for s in &missing {
                output.warn(&format!("  MISSING: {}", s.rel_path()));
            }
        }

        // Checksum verification (opt-in, reads every byte)
        let mut checksum_failures = 0;
        if checksums && manifest.schema_version >= 2 {
            let files_with_hashes: Vec<_> = manifest.files.iter()
                .filter(|(_, entry)| entry.sha256.is_some())
                .collect();

            if files_with_hashes.is_empty() {
                output.status("  checksums: none in manifest (skipped)");
            } else {
                output.status(&format!("  verifying {} checksums...", files_with_hashes.len()));
                let show_progress = std::io::stderr().is_terminal();
                let pb = if show_progress {
                    let pb = indicatif::ProgressBar::new(files_with_hashes.len() as u64);
                    pb.set_style(
                        indicatif::ProgressStyle::default_bar()
                            .template("  [{bar:40}] {pos}/{len} files")
                            .expect("hardcoded template")
                            .progress_chars("=> "),
                    );
                    Some(pb)
                } else {
                    None
                };

                for (rel_path, entry) in &files_with_hashes {
                    let local_path = local_base.join(rel_path);
                    if local_path.exists() {
                        if let Some(expected_hash) = &entry.sha256 {
                            match sha256_file(&local_path) {
                                Ok(actual_hash) => {
                                    if &actual_hash != expected_hash {
                                        output.warn(&format!("  CHECKSUM MISMATCH: {rel_path}"));
                                        checksum_failures += 1;
                                        all_ok = false;
                                    }
                                }
                                Err(e) => {
                                    output.warn(&format!("  HASH ERROR: {rel_path}: {e}"));
                                    checksum_failures += 1;
                                    all_ok = false;
                                }
                            }
                        }
                    }
                    if let Some(pb) = &pb {
                        pb.inc(1);
                    }
                }
                if let Some(pb) = &pb {
                    pb.finish_and_clear();
                }
            }
        } else if checksums && manifest.schema_version < 2 {
            output.status("  checksums: v1 manifest has no hashes (skipped)");
        }

        if mismatched.is_empty() && missing.is_empty() && checksum_failures == 0 {
            output.success(&format!("  {}: {}/{} files valid ({})",
                pack.id, valid, states.len(), human_size(manifest.total_size)));
        } else {
            output.warn(&format!("  {}: {}/{} valid, {} mismatched, {} missing, {} checksum failures",
                pack.id, valid, states.len(), mismatched.len(), missing.len(), checksum_failures));
        }
    }

    if all_ok {
        output.success("All packs verified.");
    } else {
        return Err(FavorError::DataMissing("Verification failed. Run `favor data pull` to fix.".to_string()));
    }

    Ok(())
}

pub fn sha256_file(path: &Path) -> Result<String, FavorError> {
    let mut hasher = Sha256::new();
    let mut file = std::fs::File::open(path)?;
    let mut buf = [0u8; 256 * 1024];
    loop {
        let n = file.read(&mut buf)?;
        if n == 0 { break; }
        hasher.update(&buf[..n]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn confirm_download(total_bytes: u64) -> bool {
    eprint!("  Download {}? [Y/n] ", human_size(total_bytes));
    let _ = std::io::stderr().flush();
    let mut input = String::new();
    if std::io::stdin().read_line(&mut input).is_ok() {
        let trimmed = input.trim().to_lowercase();
        trimmed.is_empty() || trimmed == "y" || trimmed == "yes"
    } else {
        false
    }
}

fn download_with_retry(
    url: &str,
    tmp_path: &Path,
    resume_from: u64,
    expected_size: u64,
    label: &str,
    show_progress: bool,
    output: &dyn Output,
) -> Result<(), FavorError> {
    let mut offset = resume_from;

    for attempt in 1..=MAX_RETRIES {
        match download_chunk(url, tmp_path, offset, expected_size, label, show_progress) {
            Ok(()) => return Ok(()),
            Err(e) => {
                if tmp_path.exists() {
                    offset = tmp_path.metadata().map(|m| m.len()).unwrap_or(offset);
                }
                if attempt < MAX_RETRIES {
                    output.warn(&format!(
                        "{label} attempt {attempt}/{MAX_RETRIES} failed: {e}, resuming from {}",
                        human_size(offset),
                    ));
                    std::thread::sleep(std::time::Duration::from_secs(2u64.pow(attempt)));
                } else {
                    return Err(e);
                }
            }
        }
    }
    unreachable!()
}

fn download_chunk(
    url: &str,
    tmp_path: &Path,
    offset: u64,
    expected_size: u64,
    label: &str,
    show_progress: bool,
) -> Result<(), FavorError> {
    let mut request = ureq::get(url);
    if offset > 0 {
        request = request.header("Range", &format!("bytes={offset}-"));
    }

    let response = request.call().map_err(|e| {
        FavorError::Resource(format!("Download failed: {e}"))
    })?;

    let actual_offset = if offset > 0 && response.status() != 206 {
        // Server doesn't support Range — must start from scratch.
        // Truncate the existing tmp file instead of appending.
        0
    } else {
        offset
    };

    let pb = if show_progress {
        let pb = indicatif::ProgressBar::new(expected_size);
        pb.set_style(
            indicatif::ProgressStyle::default_bar()
                .template("  [{bar:40}] {bytes}/{total_bytes}  {msg}  ({bytes_per_sec}, {eta})")
                .expect("hardcoded template")
                .progress_chars("=> "),
        );
        pb.set_message(label.to_string());
        pb.set_position(actual_offset);
        Some(pb)
    } else {
        None
    };

    let mut file = if actual_offset > 0 {
        std::fs::OpenOptions::new().append(true).open(tmp_path)?
    } else {
        std::fs::File::create(tmp_path)?
    };

    let mut reader = response.into_body().into_reader();
    let mut buf = [0u8; 256 * 1024];
    loop {
        let n = reader.read(&mut buf).map_err(|e| {
            FavorError::Resource(format!("Read error: {e}"))
        })?;
        if n == 0 { break; }
        file.write_all(&buf[..n])?;
        if let Some(pb) = &pb {
            pb.inc(n as u64);
        }
    }

    if let Some(pb) = &pb {
        pb.finish_and_clear();
    }
    Ok(())
}

