pub mod action;
pub mod app;
pub mod event;
pub mod output;
pub mod screen;
pub mod screens;
mod smoke;
pub mod state;
pub mod term;
pub mod theme;
pub mod widgets;

use std::path::PathBuf;
use std::sync::{mpsc, Arc, Mutex};

use crate::config::Environment;
use crate::error::CohortError;
use crate::output::{Output, OutputMode};

use self::app::{make_terminal, App};
use self::output::{LogLine, TuiOutput};
use self::screens::setup::{SetupOutcome, SetupScreen};
use self::screens::workspace::WorkspaceScreen;
use self::term::TermGuard;

pub fn run(cwd: PathBuf, _out: &dyn Output, _mode: &OutputMode) -> Result<(), CohortError> {
    let _guard = TermGuard::enter().map_err(|e| CohortError::Internal(e.into()))?;
    let mut terminal = make_terminal().map_err(|e| CohortError::Internal(e.into()))?;
    let (log_tx, log_rx) = mpsc::channel::<LogLine>();
    let bars = Arc::new(Mutex::new(Vec::new()));
    let tui_out = Arc::new(TuiOutput::new(log_tx, bars.clone()));
    if std::env::var("COHORT_TUI_SMOKE").is_ok_and(|v| !v.is_empty()) {
        smoke::spawn_demo(tui_out.clone());
    }
    let initial = Box::new(WorkspaceScreen::new(cwd));
    App::new(initial, log_rx, bars, tui_out).run(&mut terminal)
}

pub fn run_setup_only(
    initial_env: Option<Environment>,
    initial_memory: Option<String>,
) -> Result<Option<SetupOutcome>, CohortError> {
    let _guard = TermGuard::enter().map_err(|e| CohortError::Internal(e.into()))?;
    let mut terminal = make_terminal().map_err(|e| CohortError::Internal(e.into()))?;
    let (log_tx, log_rx) = mpsc::channel::<LogLine>();
    let bars = Arc::new(Mutex::new(Vec::new()));
    let tui_out = Arc::new(TuiOutput::new(log_tx, bars.clone()));
    let sink: Arc<Mutex<Option<SetupOutcome>>> = Arc::new(Mutex::new(None));
    let screen = SetupScreen::with_sink(initial_env, initial_memory, Arc::clone(&sink));
    App::new(Box::new(screen), log_rx, bars, tui_out).run(&mut terminal)?;
    let outcome = sink.lock().unwrap().take();
    Ok(outcome)
}
