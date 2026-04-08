use std::io;

use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use crossterm::ExecutableCommand;

/// RAII guard for raw mode + alternate screen. Restores terminal on Drop,
/// including panic unwind.
pub struct TermGuard;

impl TermGuard {
    pub fn enter() -> io::Result<Self> {
        enable_raw_mode()?;
        io::stdout().execute(EnterAlternateScreen)?;
        Ok(TermGuard)
    }
}

impl Drop for TermGuard {
    fn drop(&mut self) {
        let _ = io::stdout().execute(LeaveAlternateScreen);
        let _ = disable_raw_mode();
    }
}
