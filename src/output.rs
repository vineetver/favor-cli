//! Output trait + implementations. All commands use this — never write stdout directly.

use std::io::IsTerminal;
use std::sync::Arc;
use std::time::Duration;

use clap::ValueEnum;
use serde_json::json;

use crate::error::CohortError;

#[derive(Clone, Debug, ValueEnum)]
pub enum Format {
    Auto,
    Json,
    Human,
}

#[derive(Clone, Debug, PartialEq)]
pub enum OutputMode {
    Human,
    Machine,
}

impl OutputMode {
    pub fn detect(format: &Format) -> Self {
        match format {
            Format::Json => Self::Machine,
            Format::Human => Self::Human,
            Format::Auto => {
                if std::io::stdout().is_terminal() {
                    Self::Human
                } else {
                    Self::Machine
                }
            }
        }
    }

    pub fn is_machine(&self) -> bool {
        *self == Self::Machine
    }
}

#[derive(Clone)]
pub struct Progress(Option<Arc<indicatif::ProgressBar>>);

impl Progress {
    /// A no-op progress handle. Used by machine-mode output and by test
    /// `Output` impls that have no reason to render progress bars.
    pub fn noop() -> Self {
        Self(None)
    }

    #[inline]
    pub fn inc(&self, n: u64) {
        if let Some(p) = &self.0 {
            p.inc(n);
        }
    }

    pub fn finish(&self, msg: &str) {
        if let Some(p) = &self.0 {
            p.finish_with_message(msg.to_string());
        }
    }
}

/// All commands use this trait for output — never write stdout directly.
pub trait Output: Send + Sync {
    fn status(&self, msg: &str);
    fn success(&self, msg: &str);
    fn warn(&self, msg: &str);
    fn error(&self, err: &CohortError);
    fn result_json(&self, data: &serde_json::Value);
    fn table(&self, headers: &[&str], rows: &[Vec<String>]);
    /// `total == 0` is a spinner (unknown length). Machine mode is a no-op.
    fn progress(&self, total: u64, label: &str) -> Progress;
    fn is_cancelled(&self) -> bool {
        false
    }
}

pub fn bail_if_cancelled(out: &dyn Output) -> Result<(), CohortError> {
    if out.is_cancelled() {
        return Err(CohortError::Cancelled);
    }
    Ok(())
}

pub fn create(mode: &OutputMode) -> Box<dyn Output> {
    match mode {
        OutputMode::Human => Box::new(HumanOutput),
        OutputMode::Machine => Box::new(MachineOutput),
    }
}

struct HumanOutput;

impl Output for HumanOutput {
    fn status(&self, msg: &str) {
        eprintln!("  > {msg}");
    }

    fn success(&self, msg: &str) {
        eprintln!("  + {msg}");
    }

    fn warn(&self, msg: &str) {
        eprintln!("  ! {msg}");
    }

    fn error(&self, err: &CohortError) {
        eprintln!("  x {err}");
    }

    fn result_json(&self, data: &serde_json::Value) {
        if let Ok(json) = serde_json::to_string_pretty(data) {
            println!("{json}");
        }
    }

    fn progress(&self, total: u64, label: &str) -> Progress {
        let pb = if total == 0 {
            indicatif::ProgressBar::new_spinner()
        } else {
            indicatif::ProgressBar::new(total)
        };
        let style = if total == 0 {
            indicatif::ProgressStyle::default_spinner()
                .template("  {spinner} {msg} {pos} ({per_sec})")
                .unwrap()
        } else {
            indicatif::ProgressStyle::default_bar()
                .template("  {msg:<28} [{bar:30}] {pos}/{len} ({eta})")
                .unwrap()
                .progress_chars("=> ")
        };
        pb.set_style(style);
        pb.set_message(label.to_string());
        pb.enable_steady_tick(Duration::from_millis(120));
        Progress(Some(Arc::new(pb)))
    }

    fn table(&self, headers: &[&str], rows: &[Vec<String>]) {
        if headers.is_empty() {
            return;
        }
        let mut widths: Vec<usize> = headers.iter().map(|h| h.len()).collect();
        for row in rows {
            for (i, cell) in row.iter().enumerate() {
                if i < widths.len() {
                    widths[i] = widths[i].max(cell.len());
                }
            }
        }
        let header: String = headers
            .iter()
            .zip(&widths)
            .map(|(h, w)| format!("  {h:<w$}"))
            .collect::<Vec<_>>()
            .join("");
        eprintln!("{header}");
        let sep: String = widths
            .iter()
            .map(|w| format!("  {}", "-".repeat(*w)))
            .collect::<Vec<_>>()
            .join("");
        eprintln!("{sep}");
        for row in rows {
            let line: String = row
                .iter()
                .zip(&widths)
                .map(|(cell, w)| format!("  {cell:<w$}"))
                .collect::<Vec<_>>()
                .join("");
            eprintln!("{line}");
        }
    }
}

struct MachineOutput;

impl Output for MachineOutput {
    fn status(&self, msg: &str) {
        eprintln!("{}", json!({"level": "info", "message": msg}));
    }

    fn success(&self, msg: &str) {
        eprintln!("{}", json!({"level": "info", "message": msg}));
    }

    fn warn(&self, msg: &str) {
        eprintln!("{}", json!({"level": "warn", "message": msg}));
    }

    fn error(&self, err: &CohortError) {
        eprintln!(
            "{}",
            json!({"error": err.code_name(), "message": err.to_string(), "exit_code": err.exit_code()})
        );
    }

    fn result_json(&self, data: &serde_json::Value) {
        if let Ok(json) = serde_json::to_string(data) {
            println!("{json}");
        }
    }

    fn progress(&self, _total: u64, _label: &str) -> Progress {
        Progress::noop()
    }

    fn table(&self, headers: &[&str], rows: &[Vec<String>]) {
        let result: Vec<serde_json::Value> = rows
            .iter()
            .map(|row| {
                let mut obj = serde_json::Map::new();
                for (i, cell) in row.iter().enumerate() {
                    let key = headers.get(i).unwrap_or(&"");
                    obj.insert(key.to_string(), serde_json::Value::String(cell.clone()));
                }
                serde_json::Value::Object(obj)
            })
            .collect();
        if let Ok(json) = serde_json::to_string(&result) {
            println!("{json}");
        }
    }
}
