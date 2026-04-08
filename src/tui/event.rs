use crossterm::event::KeyEvent;

use crate::error::CohortError;

use super::output::{LogLine, ProgressSnapshot};

pub enum AppEvent {
    Key(KeyEvent),
    Tick,
    Log(LogLine),
    ProgressUpdate(ProgressSnapshot),
    ProgressSweep(u64),
    CommandDone(Result<(), CohortError>),
}
