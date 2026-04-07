use std::path::PathBuf;

use crate::config::{Config, Environment, ResourceConfig};

/// Detected system resources for DataFusion and parallel execution.
///
/// Memory resolution priority:
/// 1. SLURM allocation (SLURM_MEM_PER_NODE) — srun jobs get full allocation
/// 2. Config memory_budget — user's chosen default from setup
/// 3. Auto-detect — cgroup → /proc/meminfo → 4GB fallback
#[derive(Debug, Clone)]
pub struct Resources {
    pub memory_bytes: u64,
    pub threads: usize,
    pub temp_dir: PathBuf,
    /// Configured environment, if set during setup
    environment: Option<Environment>,
}

impl Resources {
    /// Load config and detect resources in one call.
    /// If config doesn't exist or can't be read, falls back to auto-detect.
    pub fn detect_configured() -> Self {
        let config_resources = Config::load().map(|c| c.resources).unwrap_or_default();
        Self::detect_with_config(&config_resources)
    }

    /// Auto-detect resources without config. Used during setup (before config exists).
    pub fn detect() -> Self {
        Self {
            memory_bytes: detect_memory(),
            threads: detect_threads(None),
            temp_dir: detect_temp(None),
            environment: None,
        }
    }

    /// Detect resources respecting user's configured budget.
    ///
    /// - If inside a SLURM job, uses the SLURM allocation (overrides budget).
    /// - Otherwise, uses the configured memory_budget as the ceiling.
    /// - Falls back to auto-detect if no budget is configured.
    pub fn detect_with_config(config: &ResourceConfig) -> Self {
        let memory_bytes = if let Some(bytes) = slurm_memory() {
            bytes // 95% of SLURM allocation
        } else if let Some(budget) = config.memory_budget_bytes() {
            budget * 90 / 100
        } else {
            detect_memory()
        };

        Self {
            memory_bytes,
            threads: detect_threads(config.threads),
            temp_dir: detect_temp(Some(config)),
            environment: config.environment,
        }
    }

    pub fn memory_human(&self) -> String {
        let gb = self.memory_bytes as f64 / (1024.0 * 1024.0 * 1024.0);
        format!("{gb:.1} GiB")
    }

    pub fn environment(&self) -> &'static str {
        if std::env::var("SLURM_JOB_ID").is_ok() {
            return "SLURM job";
        }
        if let Some(env) = self.environment {
            return match env {
                Environment::Hpc => "HPC",
                Environment::Workstation => "workstation",
            };
        }
        if let Ok(hostname) = std::env::var("HOSTNAME") {
            if hostname.contains("login") {
                return "login node";
            }
            return "compute";
        }
        "local"
    }
}

/// SLURM allocation in bytes. Uses 95% — leaves headroom for slurmstepd + kernel.
fn slurm_memory() -> Option<u64> {
    let val = std::env::var("SLURM_MEM_PER_NODE").ok()?;
    let mb: u64 = val.trim_end_matches('M').parse().ok()?;
    Some(mb * 1024 * 1024 * 95 / 100)
}

/// cgroup memory limit, if any. 90% of the limit.
fn cgroup_memory() -> Option<u64> {
    const PATHS: &[&str] = &[
        "/sys/fs/cgroup/memory/memory.limit_in_bytes",
        "/sys/fs/cgroup/memory.max",
    ];
    PATHS.iter().find_map(|p| {
        let content = std::fs::read_to_string(p).ok()?;
        let trimmed = content.trim();
        if trimmed == "max" || trimmed == "9223372036854771712" {
            return None;
        }
        let bytes: u64 = trimmed.parse().ok()?;
        // Reject implausibly large values (1 TiB ceiling on cgroup-reported memory).
        (bytes < 1024 * 1024 * 1024 * 1024).then_some(bytes * 90 / 100)
    })
}

/// /proc/meminfo MemAvailable. 90% to leave room for OS.
fn meminfo_memory() -> Option<u64> {
    let content = std::fs::read_to_string("/proc/meminfo").ok()?;
    let kb: u64 = content
        .lines()
        .find_map(|l| l.strip_prefix("MemAvailable:"))?
        .split_whitespace()
        .next()?
        .parse()
        .ok()?;
    Some(kb * 1024 * 90 / 100)
}

fn detect_memory() -> u64 {
    slurm_memory()
        .or_else(cgroup_memory)
        .or_else(meminfo_memory)
        .unwrap_or(4 * 1024 * 1024 * 1024)
}

fn env_usize(name: &str) -> Option<usize> {
    std::env::var(name).ok()?.parse().ok()
}

fn cgroup_threads() -> Option<usize> {
    let q: i64 = std::fs::read_to_string("/sys/fs/cgroup/cpu/cpu.cfs_quota_us")
        .ok()?
        .trim()
        .parse()
        .ok()?;
    let p: i64 = std::fs::read_to_string("/sys/fs/cgroup/cpu/cpu.cfs_period_us")
        .ok()?
        .trim()
        .parse()
        .ok()?;
    (q > 0 && p > 0).then(|| ((q / p).max(1)) as usize)
}

fn detect_threads(config_threads: Option<usize>) -> usize {
    env_usize("COHORT_THREADS")
        .or_else(|| env_usize("SLURM_CPUS_PER_TASK"))
        .or(config_threads)
        .or_else(cgroup_threads)
        .unwrap_or_else(|| std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4))
}

fn detect_temp(config: Option<&ResourceConfig>) -> PathBuf {
    let existing = |s: &str| (!s.is_empty()).then(|| PathBuf::from(s)).filter(|p| p.exists());

    config
        .and_then(|c| existing(&c.temp_dir))
        .or_else(|| ["TMPDIR", "SCRATCH", "LOCAL_SCRATCH"].iter()
            .find_map(|v| std::env::var(v).ok().and_then(|s| existing(&s))))
        .or_else(|| existing("/scratch"))
        .unwrap_or_else(std::env::temp_dir)
}
