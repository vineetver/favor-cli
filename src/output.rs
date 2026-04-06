//! Output trait + implementations. All commands use this — never write stdout directly.

use std::io::IsTerminal;

use clap::ValueEnum;
use serde_json::json;

use crate::error::FavorError;

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
        if std::env::var("FAVOR_MACHINE").is_ok_and(|v| !v.is_empty()) {
            return Self::Machine;
        }
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

/// All commands use this trait for output — never write stdout directly.
pub trait Output {
    fn status(&self, msg: &str);
    fn success(&self, msg: &str);
    fn warn(&self, msg: &str);
    fn error(&self, err: &FavorError);
    fn result_json(&self, data: &serde_json::Value);
    fn table(&self, headers: &[&str], rows: &[Vec<String>]);
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

    fn error(&self, err: &FavorError) {
        eprintln!("  x {err}");
    }

    fn result_json(&self, data: &serde_json::Value) {
        if let Ok(json) = serde_json::to_string_pretty(data) {
            println!("{json}");
        }
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

    fn error(&self, err: &FavorError) {
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
