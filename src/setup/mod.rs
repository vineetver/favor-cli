mod tui;

use std::path::PathBuf;

use serde_json::json;

use crate::cli::DataAction;
use crate::config::{Config, DirProbe, Environment, ResourceConfig, Tier};
use crate::data::Pack;
use crate::error::FavorError;
use crate::output::Output;
use crate::output::OutputMode;
use crate::resource::Resources;

pub fn init(
    path: Option<PathBuf>,
    force: bool,
    out: &dyn Output,
    mode: &OutputMode,
) -> Result<(), FavorError> {
    let config = Config::load_configured()?;

    let project_dir = match path {
        Some(p) => {
            std::fs::create_dir_all(&p).map_err(|e| {
                FavorError::Resource(format!("Cannot create directory '{}': {e}", p.display()))
            })?;
            p.canonicalize()
                .map_err(|e| FavorError::Input(format!("Cannot resolve '{}': {e}", p.display())))?
        }
        None => std::env::current_dir().map_err(|e| {
            FavorError::Resource(format!("Cannot determine current directory: {e}"))
        })?,
    };

    let claude_path = project_dir.join("CLAUDE.md");
    let codex_dir = project_dir.join(".codex");
    let codex_path = codex_dir.join("instructions.md");

    let is_refresh = claude_path.exists();

    if is_refresh && !force && !mode.is_machine() {
        out.warn(
            "CLAUDE.md already exists. Use --force to overwrite, or re-run to refresh after installing new packs.",
        );
        return Ok(());
    }

    let probe = DirProbe::scan(&config.root_dir());
    let content = render_init_template(&config, &probe);

    std::fs::write(&claude_path, &content).map_err(|e| {
        FavorError::Resource(format!("Cannot write '{}': {e}", claude_path.display()))
    })?;
    std::fs::create_dir_all(&codex_dir).map_err(|e| {
        FavorError::Resource(format!("Cannot create '{}': {e}", codex_dir.display()))
    })?;
    std::fs::write(&codex_path, &content).map_err(|e| {
        FavorError::Resource(format!("Cannot write '{}': {e}", codex_path.display()))
    })?;

    let action = if is_refresh {
        "Refreshed"
    } else {
        "Initialized"
    };
    out.success(&format!(
        "{action} agent context in {}",
        project_dir.display()
    ));
    out.status("  CLAUDE.md (Claude Code)");
    out.status("  .codex/instructions.md (Codex)");
    out.status(&format!(
        "Open your coding agent in {} and ask a research question.",
        project_dir.display()
    ));

    out.result_json(&json!({
        "project_dir": project_dir.to_string_lossy(),
        "files": [claude_path.to_string_lossy(), codex_path.to_string_lossy()],
        "tier": config.data.tier.as_str(),
        "packs": config.data.packs,
    }));

    Ok(())
}

fn render_init_template(config: &Config, probe: &DirProbe) -> String {
    let pack_list = build_pack_list(config, probe);

    let environment = config
        .resources
        .environment
        .map(|e| e.as_str().to_string())
        .unwrap_or_else(|| "auto-detect".to_string());

    let memory_budget = config
        .resources
        .memory_budget
        .clone()
        .unwrap_or_else(|| "auto-detect".to_string());

    let chrom_count = match config.data.tier {
        Tier::Base => probe.base_chroms,
        Tier::Full => probe.full_chroms,
    };

    INIT_TEMPLATE
        .replace("{root_dir}", &config.data.root_dir)
        .replace("{tier}", config.data.tier.as_str())
        .replace("{tier_size}", config.data.tier.size_human())
        .replace("{chrom_count}", &chrom_count.to_string())
        .replace("{environment}", &environment)
        .replace("{memory_budget}", &memory_budget)
        .replace("{pack_list}", &pack_list)
}

fn build_pack_list(config: &Config, probe: &DirProbe) -> String {
    let lines: Vec<String> = Pack::all()
        .iter()
        .filter(|p| !p.always_installed)
        .filter(|p| {
            config.data.packs.contains(&p.id.to_string())
                || p.tables.iter().any(|t| probe.tissue_tables.contains(&t.to_string()))
        })
        .map(|p| format!("- **{}** ({}): {}", p.id, p.size_human, p.description))
        .collect();
    if lines.is_empty() {
        return "- No optional packs. Run `favor data pull --pack <id>` then `favor init --force`.".into();
    }
    lines.join("\n")
}

const INIT_TEMPLATE: &str = r#"# FAVOR Genomic Analysis

Always pass `--format json` to every favor command. Use `--dry-run` before heavy computation.

## System

- Data root: {root_dir}
- Tier: {tier} ({tier_size}, {chrom_count}/24 chromosomes)
- Environment: {environment}
- Memory budget: {memory_budget}

## Installed packs

{pack_list}

## Pipeline

```
variant file -> favor ingest -> favor annotate -> favor enrich --tissue X
                                                -> favor staar (rare-variant association)
                                                -> favor interpret (variant-to-gene)
```

## Commands

| Command | What it does |
|---------|-------------|
| `favor ingest <file>` | Normalize VCF/TSV/CSV to parquet with variant ID |
| `favor annotate <file>` | Add CADD, ClinVar, gnomAD, REVEL, aPC scores, regulatory marks |
| `favor enrich <file> --tissue <name>` | Add tissue eQTL, ChromBPNet, enhancer-gene links |
| `favor staar --genotypes <vcf> --phenotype <tsv> --trait-name <col> --annotations <parquet>` | Rare-variant burden test (STAAR-O) |
| `favor schema [table]` | Show column names and types |
| `favor manifest` | List installed data and available commands |

## STAAR usage

```bash
favor staar --dry-run --format json \
  --genotypes cohort.vcf.gz \
  --phenotype pheno.tsv \
  --trait-name LDL \
  --covariates age,sex,PC1,PC2 \
  --annotations annotated.parquet \
  --masks coding
```

Use `--dry-run` first to check memory. On HPC: `srun --mem=64G -c 8 favor staar ...`

## Interpreting results

- **CADD phred > 20**: top 1% most deleterious variants genome-wide
- **REVEL > 0.5**: likely pathogenic missense
- **STAAR-O < 2.5e-6**: genome-wide significant (Bonferroni for ~20K genes)
- **pLoF mask significant**: protein-destroying variants cause the trait (strongest evidence)
- **Synonymous mask significant**: suspicious — likely confounding, not biology

STAAR runs 6 tests (Burden, SKAT, ACAT-V x two beta weights) across 11 annotation channels:
- Burden: all variants push trait same direction
- SKAT: mixed-effect variants
- ACAT-V: one or two variants drive everything

## Querying parquet output

Use DataFusion or any parquet-capable SQL engine:

```sql
-- Damaging variants
SELECT vid, gencode.genes[1] AS gene, main.cadd.phred
FROM 'annotated.parquet' WHERE main.cadd.phred > 20;

-- Significant genes
SELECT * FROM 'staar_results/coding_pLoF_missense.parquet'
WHERE "STAAR-O" < 2.5e-6 ORDER BY "STAAR-O";

-- Tissue enrichment
SELECT a.vid, e.* FROM 'enriched/annotated.parquet' a
JOIN 'enriched/eqtl.parquet' e ON a.vid = e.vid
WHERE e.tissue_name LIKE '%Liver%';
```

## Output conventions

- Parquet files: main results
- `.meta.json`: parameters and counts
- stdout (--format json): structured result summary
- stderr: progress/status messages
- Exit codes: 0=ok, 1=input, 2=data missing, 3=resource, 4=analysis
"#;

pub fn setup(
    output: &dyn Output,
    mode: &OutputMode,
    cli_environment: Option<String>,
    cli_memory_budget: Option<String>,
) -> Result<(), FavorError> {
    if mode.is_machine() {
        return Err(FavorError::Input(
            "setup requires interactive mode — run without --format json".to_string(),
        ));
    }

    // 1. Pick tier — returns Tier directly, no index mapping
    let tier = match tui::select_tier().map_err(|e| FavorError::Internal(e.into()))? {
        Some(t) => t,
        None => {
            output.warn("Setup cancelled");
            return Ok(());
        }
    };

    // 2. Pick root directory — browser shows live data probe
    let cwd = std::env::current_dir().unwrap_or_else(|_| Config::default_root_dir());
    let root = match tui::select_directory("Select FAVOR root directory", &cwd)
        .map_err(|e| FavorError::Internal(e.into()))?
    {
        Some(p) => p,
        None => {
            output.warn("Setup cancelled");
            return Ok(());
        }
    };

    // 3. Pick add-on packs — detect what's already installed on disk
    let optional_packs = Pack::optional();
    let installed_pack_ids: Vec<String> = optional_packs
        .iter()
        .filter(|p| p.tables.iter().any(|t| p.local_dir(&root).join(t).is_dir()))
        .map(|p| p.id.to_string())
        .collect();
    let selected_packs = tui::select_packs(&optional_packs, &installed_pack_ids)
        .map_err(|e| FavorError::Internal(e.into()))?
        .unwrap_or_default(); // Esc = skip, not cancel

    // 4. Environment selection — HPC or workstation?
    let environment = if let Some(env_str) = &cli_environment {
        Some(env_str.parse::<Environment>()?)
    } else {
        tui::select_environment().map_err(|e| FavorError::Internal(e.into()))?
    };

    if environment == Some(Environment::Hpc) {
        let has_srun = std::process::Command::new("which")
            .arg("srun")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false);
        let has_sbatch = std::process::Command::new("which")
            .arg("sbatch")
            .stdout(std::process::Stdio::null())
            .stderr(std::process::Stdio::null())
            .status()
            .map(|s| s.success())
            .unwrap_or(false);

        if !has_srun && !has_sbatch {
            output.warn("srun/sbatch not found in PATH — are you sure this is HPC?");
            output.warn("Continuing anyway. You can re-run `favor setup` to change.");
        }
    }

    // 5. Memory budget selection
    let res_detect = Resources::detect();
    let memory_budget = if let Some(budget) = &cli_memory_budget {
        if ResourceConfig::parse_memory_bytes(budget).is_none() {
            return Err(FavorError::Input(format!(
                "Cannot parse memory budget '{budget}'. Use format like '16GB', '64G', '8192MB'."
            )));
        }
        Some(budget.clone())
    } else {
        tui::select_memory_budget(&res_detect).map_err(|e| FavorError::Internal(e.into()))?
    };

    // 6. Probe the chosen root to show what's already there
    let probe = DirProbe::scan(&root);

    output.status("FAVOR Configuration");
    output.status(&format!("  Root:       {}", root.display()));
    output.status(&format!("  Tier:       {tier}"));

    // Annotation status from probe
    let chrom_count = match tier {
        Tier::Full => probe.full_chroms,
        Tier::Base => probe.base_chroms,
    };
    if chrom_count == 24 {
        output.success(&format!("  {tier}: 24/24 chromosomes found"));
    } else if chrom_count > 0 {
        output.status(&format!(
            "  {tier}: {chrom_count}/24 found, will download remaining"
        ));
    } else {
        output.status(&format!(
            "  {tier}: will download ({tier_size})",
            tier_size = tier.size_human()
        ));
    }

    if !probe.tissue_tables.is_empty() {
        output.success(&format!(
            "  Tissue: {} tables found",
            probe.tissue_tables.len()
        ));
    }
    if !selected_packs.is_empty() {
        output.status(&format!("  Packs:  {}", selected_packs.join(", ")));
    }
    if let Some(env) = environment {
        output.status(&format!("  Env:    {}", env));
    }
    if let Some(budget) = &memory_budget {
        output.status(&format!("  Budget: {}", budget));
    }
    output.status(&format!(
        "  System: {} memory, {} threads ({})",
        res_detect.memory_human(),
        res_detect.threads,
        res_detect.environment()
    ));

    // 7. Save config — preserve existing resource customizations, overlay new ones
    let mut existing_resources = Config::load().map(|c| c.resources).unwrap_or_default();

    if environment.is_some() {
        existing_resources.environment = environment;
    }
    if memory_budget.is_some() {
        existing_resources.memory_budget = memory_budget;
    }

    let config = Config {
        data: crate::config::DataConfig {
            tier,
            root_dir: root.to_string_lossy().to_string(),
            packs: selected_packs.clone(),
        },
        resources: existing_resources,
    };
    config.save()?;
    output.success(&format!("Saved to {}", Config::config_path().display()));

    // 8. Download annotations (skips if already complete)
    output.status("Checking annotations...");
    if let Err(e) = crate::data::transfer::run(
        DataAction::Pull {
            full: tier == Tier::Full,
            dry_run: false,
            yes: true,
            pack: None,
        },
        output,
    ) {
        output.warn(
            "Config saved but annotation download incomplete. Run `favor data pull` to retry.",
        );
        return Err(e);
    }

    // 9. Download always-installed packs
    for pack in Pack::required() {
        if let Err(e) = crate::data::transfer::pull_pack(pack.id, false, true, output) {
            output.warn(&format!(
                "{} failed: {e}. Run `favor data pull --pack {}` to retry.",
                pack.name, pack.id
            ));
        }
    }

    // 10. Download selected packs
    for pack_id in &selected_packs {
        if let Err(e) = crate::data::transfer::pull_pack(pack_id, false, true, output) {
            output.warn(&format!(
                "Pack '{pack_id}' failed: {e}. Run `favor data pull --pack {pack_id}` to retry."
            ));
        }
    }

    output.success("Setup complete. Run: favor annotate <input.vcf>");
    Ok(())
}

pub fn uninstall(out: &dyn Output) -> Result<(), FavorError> {
    let binary = std::env::current_exe().unwrap_or_default();
    let config_dir = Config::config_dir();

    out.status("This will remove:");
    out.status(&format!("  Binary: {}", binary.display()));
    out.status(&format!("  Config: {}", config_dir.display()));
    out.warn("Data packs at your configured root directory will NOT be deleted.");

    if config_dir.exists() {
        std::fs::remove_dir_all(&config_dir).map_err(|e| {
            FavorError::Resource(format!(
                "Cannot remove config directory '{}': {e}",
                config_dir.display()
            ))
        })?;
        out.status("  Removed config directory");
    }

    let home = dirs::home_dir().unwrap_or_default();
    let install_dir = binary.parent().unwrap_or(std::path::Path::new(""));
    let install_str = install_dir.to_string_lossy();

    for rc in &[".bashrc", ".zshrc"] {
        let rc_path = home.join(rc);
        if rc_path.exists() {
            if let Ok(content) = std::fs::read_to_string(&rc_path) {
                let filtered: Vec<&str> = content
                    .lines()
                    .filter(|line| !line.contains(&*install_str) || !line.contains("PATH"))
                    .collect();
                if filtered.len() < content.lines().count() {
                    let _ = std::fs::write(&rc_path, filtered.join("\n") + "\n");
                    out.status(&format!("  Cleaned PATH from {}", rc));
                }
            }
        }
    }

    // Remove binary last (we're running it)
    if binary.exists() {
        std::fs::remove_file(&binary).map_err(|e| {
            FavorError::Resource(format!("Cannot remove binary '{}': {e}", binary.display()))
        })?;
    }

    out.success("Uninstalled. Data packs remain at your configured root directory.");
    Ok(())
}
