use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::Sender;
use std::sync::{Arc, Mutex};

use indicatif::{ProgressBar, ProgressDrawTarget};

use crate::error::CohortError;
use crate::output::{Output, Progress};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LogLevel {
    Status,
    Success,
    Warn,
    Error,
}

#[derive(Debug, Clone)]
pub struct LogLine {
    pub level: LogLevel,
    pub message: String,
}

#[derive(Debug, Clone)]
pub struct ProgressSnapshot {
    pub id: u64,
    pub position: u64,
    pub length: u64,
    pub message: String,
    pub indeterminate: bool,
    pub finished: bool,
}

pub type BarRegistry = Arc<Mutex<Vec<Arc<ProgressBar>>>>;

pub struct TuiOutput {
    log_tx: Sender<LogLine>,
    bars: BarRegistry,
    cancel: Mutex<Arc<AtomicBool>>,
}

impl TuiOutput {
    pub fn new(log_tx: Sender<LogLine>, bars: BarRegistry) -> Self {
        Self {
            log_tx,
            bars,
            cancel: Mutex::new(Arc::new(AtomicBool::new(false))),
        }
    }

    pub fn arm_cancel(&self) -> Arc<AtomicBool> {
        let fresh = Arc::new(AtomicBool::new(false));
        *self.cancel.lock().unwrap() = Arc::clone(&fresh);
        fresh
    }

    fn cancel_handle(&self) -> Arc<AtomicBool> {
        Arc::clone(&self.cancel.lock().unwrap())
    }

    fn emit(&self, level: LogLevel, message: String) {
        let _ = self.log_tx.send(LogLine { level, message });
    }
}

impl Output for TuiOutput {
    fn status(&self, msg: &str) {
        self.emit(LogLevel::Status, msg.to_string());
    }

    fn success(&self, msg: &str) {
        self.emit(LogLevel::Success, msg.to_string());
    }

    fn warn(&self, msg: &str) {
        self.emit(LogLevel::Warn, msg.to_string());
    }

    fn error(&self, err: &CohortError) {
        self.emit(LogLevel::Error, err.to_string());
    }

    fn result_json(&self, data: &serde_json::Value) {
        let _ = data;
    }

    fn table(&self, headers: &[&str], rows: &[Vec<String>]) {
        let _ = (headers, rows);
    }

    fn progress(&self, total: u64, label: &str) -> Progress {
        // total == 0 means indeterminate; sentinel length u64::MAX is flagged as such
        // during sampling. new_spinner is avoided because its position() returns tick
        // count, not caller-driven progress.
        let pb = if total == 0 {
            ProgressBar::new(u64::MAX)
        } else {
            ProgressBar::new(total)
        };
        pb.set_draw_target(ProgressDrawTarget::hidden());
        pb.set_message(label.to_string());
        let arc = Arc::new(pb);
        self.bars.lock().unwrap().push(arc.clone());
        Progress::from_arc(arc)
    }

    fn is_cancelled(&self) -> bool {
        self.cancel_handle().load(Ordering::Relaxed)
    }
}
