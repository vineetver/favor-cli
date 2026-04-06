//! Pack publishing: admin-only operations for uploading annotation packs to MinIO.

use std::collections::HashMap;

use serde_json::json;

use crate::error::FavorError;
use crate::output::Output;

use super::transfer::{sha256_file, FileEntry, PackManifest, MANIFEST_SCHEMA_VERSION, REMOTE_BASE_URL};
use super::Pack;

const MINIO_ALIAS: &str = "favor";
const MINIO_BUCKET: &str = "favor-hg38";

pub fn publish(
    source_root: std::path::PathBuf,
    pack_filter: Option<String>,
    dry_run: bool,
    output: &dyn Output,
) -> Result<(), FavorError> {
    if !source_root.is_dir() {
        return Err(FavorError::Input(format!(
            "Source root does not exist: {}", source_root.display()
        )));
    }

    let packs: Vec<&Pack> = if let Some(ref id) = pack_filter {
        let pack = Pack::find(id).ok_or_else(|| {
            FavorError::Input(format!("Unknown pack '{id}'"))
        })?;
        vec![pack]
    } else {
        Pack::all().iter().collect()
    };

    output.status(&format!(
        "Publishing {} pack(s) from {}{}",
        packs.len(), source_root.display(),
        if dry_run { " (dry run)" } else { "" },
    ));

    for pack in &packs {
        publish_pack(pack, &source_root, dry_run, output)?;
    }

    if dry_run {
        output.success("Dry run complete. No files uploaded.");
    } else {
        output.success("All packs published.");
    }

    Ok(())
}

fn publish_pack(
    pack: &Pack,
    source_root: &std::path::Path,
    dry_run: bool,
    output: &dyn Output,
) -> Result<(), FavorError> {
    output.status(&format!("Pack '{}' — scanning...", pack.id));

    let source_base = pack.source_dir(source_root);
    if !source_base.is_dir() {
        return Err(FavorError::Input(format!(
            "Source directory does not exist for pack '{}': {}",
            pack.id, source_base.display()
        )));
    }

    let mut files: HashMap<String, FileEntry> = HashMap::new();

    for table in pack.tables {
        let table_dir = source_base.join(table);
        if !table_dir.is_dir() {
            output.warn(&format!("  table '{}' not found at {}", table, table_dir.display()));
            continue;
        }

        let walker = walkdir::WalkDir::new(&table_dir)
            .follow_links(true)
            .into_iter()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_type().is_file())
            .filter(|e| {
                let name = e.file_name().to_string_lossy();
                !name.starts_with('.') && !name.ends_with(".tmp")
            });

        for entry in walker {
            let rel_path = entry.path()
                .strip_prefix(&source_base)
                .map_err(|_| FavorError::Resource(
                    format!("Failed to compute relative path for {}", entry.path().display())
                ))?
                .to_string_lossy()
                .to_string();

            let size = entry.metadata()
                .map_err(|e| FavorError::Resource(
                    format!("Failed to read metadata for {}: {e}", entry.path().display())
                ))?
                .len();

            output.status(&format!("  hashing {rel_path}..."));
            let sha256 = sha256_file(entry.path())?;

            files.insert(rel_path, FileEntry { size, sha256: Some(sha256) });
        }
    }

    if files.is_empty() {
        output.warn(&format!("  Pack '{}': no files found, skipping", pack.id));
        return Ok(());
    }

    let total_size: u64 = files.values().map(|f| f.size).sum();

    let manifest = PackManifest {
        schema_version: MANIFEST_SCHEMA_VERSION,
        pack_id: pack.id.to_string(),
        version: pack.version.to_string(),
        base_dir: pack.base_dir.to_string(),
        total_size,
        files: files.clone(),
    };

    let manifest_json = serde_json::to_string_pretty(&manifest).map_err(|e| {
        FavorError::Resource(format!("Failed to serialize manifest: {e}"))
    })?;

    let total_gb = total_size as f64 / (1024.0 * 1024.0 * 1024.0);
    output.status(&format!("  {}: {} files, {:.1} GB, schema v{}", pack.id, manifest.files.len(), total_gb, MANIFEST_SCHEMA_VERSION));

    if dry_run {
        let manifest_path = format!("/tmp/favor-manifest-{}.json", pack.id);
        std::fs::write(&manifest_path, &manifest_json)
            .map_err(|e| FavorError::Resource(format!("Cannot write '{}': {e}", manifest_path)))?;
        output.status(&format!("  manifest written to {manifest_path}"));
        output.result_json(&json!({
            "pack_id": pack.id, "schema_version": MANIFEST_SCHEMA_VERSION,
            "files": manifest.files.len(), "total_size": total_size, "manifest_path": manifest_path,
        }));
        return Ok(());
    }

    output.status(&format!("  uploading {} files...", files.len()));
    let minio_pack_prefix = format!("{MINIO_ALIAS}/{MINIO_BUCKET}/packs/{}", pack.id);

    for (i, (rel_path, _entry)) in files.iter().enumerate() {
        let local_path = source_base.join(rel_path);
        let remote_path = format!("{minio_pack_prefix}/{rel_path}");

        let status = std::process::Command::new("mc")
            .args(["cp", &local_path.to_string_lossy(), &remote_path])
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .map_err(|e| FavorError::Resource(format!("Failed to run mc: {e}")))?;

        if !status.success() {
            return Err(FavorError::Resource(format!(
                "mc cp failed for {rel_path} (exit {})", status.code().unwrap_or(-1)
            )));
        }

        if (i + 1) % 10 == 0 || i + 1 == files.len() {
            output.status(&format!("  [{}/{}] uploaded", i + 1, files.len()));
        }
    }

    let tmp_manifest = format!("/tmp/favor-manifest-{}.json", pack.id);
    std::fs::write(&tmp_manifest, &manifest_json)
        .map_err(|e| FavorError::Resource(format!("Cannot write '{}': {e}", tmp_manifest)))?;

    let manifest_remote = format!("{minio_pack_prefix}/manifest.json");
    let status = std::process::Command::new("mc")
        .args(["cp", &tmp_manifest, &manifest_remote])
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map_err(|e| FavorError::Resource(format!("Failed to upload manifest: {e}")))?;

    if !status.success() {
        return Err(FavorError::Resource("mc cp failed for manifest.json".to_string()));
    }

    let _ = std::fs::remove_file(&tmp_manifest);

    output.status("  validating uploaded manifest...");
    let base = pack.base_url(REMOTE_BASE_URL);
    let url = format!("{base}/packs/{}/manifest.json", pack.id);
    let body = ureq::get(&url).call()
        .map_err(|e| FavorError::Resource(format!("Failed to re-fetch manifest after upload: {e}")))?
        .into_body()
        .read_to_string()
        .map_err(|e| FavorError::Resource(format!("Failed to read uploaded manifest: {e}")))?;

    let reread: PackManifest = serde_json::from_str(&body).map_err(|e| {
        FavorError::Resource(format!("Uploaded manifest is not valid v2 JSON: {e}"))
    })?;

    if reread.schema_version != MANIFEST_SCHEMA_VERSION {
        return Err(FavorError::Resource(format!(
            "Uploaded manifest has schema_version {}, expected {}", reread.schema_version, MANIFEST_SCHEMA_VERSION
        )));
    }
    if reread.files.len() != files.len() {
        return Err(FavorError::Resource(format!(
            "Uploaded manifest has {} files, expected {}", reread.files.len(), files.len()
        )));
    }

    output.success(&format!("  {} published: {} files, {:.1} GB", pack.id, files.len(), total_gb));
    Ok(())
}
