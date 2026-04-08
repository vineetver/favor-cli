use std::collections::VecDeque;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Modifier, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Clear, Paragraph, Widget};

use crate::tui::output::{LogLevel, LogLine, ProgressSnapshot};
use crate::tui::screen::RunRequest;
use crate::tui::theme;

const TAIL_LINES: usize = 5;

pub enum RunOutcome {
    Cancelled,
    Succeeded { artifact: PathBuf },
    Failed { error: String },
}

enum Phase {
    Running,
    Succeeded,
    Failed(String),
}

pub struct RunOverlay {
    title: String,
    artifact: PathBuf,
    cancel: Arc<AtomicBool>,
    phase: Phase,
    cancel_requested: bool,
    stage: String,
    progress: Option<ProgressSnapshot>,
    tail: VecDeque<LogLine>,
    finished: bool,
}

impl RunOverlay {
    pub fn launch(request: &RunRequest, cancel: Arc<AtomicBool>) -> Self {
        Self {
            title: request.description(),
            artifact: request.expected_artifact(),
            cancel,
            phase: Phase::Running,
            cancel_requested: false,
            stage: "starting".to_string(),
            progress: None,
            tail: VecDeque::with_capacity(TAIL_LINES),
            finished: false,
        }
    }

    pub fn is_finished(&self) -> bool {
        self.finished
    }

    pub fn push_log(&mut self, line: LogLine) {
        if matches!(line.level, LogLevel::Status | LogLevel::Success) {
            self.stage = line.message.clone();
        }
        if self.tail.len() == TAIL_LINES {
            self.tail.pop_front();
        }
        self.tail.push_back(line);
    }

    pub fn push_progress(&mut self, snap: ProgressSnapshot) {
        if !snap.message.is_empty() {
            self.stage = snap.message.clone();
        }
        self.progress = Some(snap);
    }

    pub fn on_command_done(&mut self, res: Result<(), String>) {
        self.phase = match res {
            Ok(()) => {
                if self.artifact.exists() {
                    Phase::Succeeded
                } else {
                    Phase::Failed(format!(
                        "command reported success but artifact missing: {}",
                        self.artifact.display()
                    ))
                }
            }
            Err(e) => Phase::Failed(e),
        };
    }

    pub fn handle(&mut self, key: KeyEvent) -> Option<RunOutcome> {
        let none = KeyModifiers::NONE;
        match (&self.phase, key.code, key.modifiers) {
            (Phase::Running, KeyCode::Char('c'), m) if m == none => {
                self.cancel_requested = true;
                self.cancel.store(true, Ordering::Relaxed);
                None
            }
            (Phase::Succeeded, KeyCode::Enter, _) => {
                self.finished = true;
                Some(RunOutcome::Succeeded {
                    artifact: self.artifact.clone(),
                })
            }
            (Phase::Succeeded, KeyCode::Esc, _) => {
                self.finished = true;
                Some(RunOutcome::Succeeded {
                    artifact: self.artifact.clone(),
                })
            }
            (Phase::Failed(msg), KeyCode::Enter, _) | (Phase::Failed(msg), KeyCode::Esc, _) => {
                self.finished = true;
                Some(RunOutcome::Failed {
                    error: msg.clone(),
                })
            }
            (Phase::Failed(_), KeyCode::Char('r'), m) if m == none => {
                self.finished = true;
                Some(RunOutcome::Cancelled)
            }
            _ => None,
        }
    }

    pub fn render(&self, area: Rect, buf: &mut Buffer) {
        let region = centered_rect(area, 80, 60);
        Clear.render(region, buf);
        let (status_label, status_color) = match &self.phase {
            Phase::Running => ("running", theme::ACCENT),
            Phase::Succeeded => ("done", theme::OK),
            Phase::Failed(_) => ("failed", theme::BAD),
        };
        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" {} [{status_label}] ", self.title))
            .border_style(Style::default().fg(status_color));
        let inner = block.inner(region);
        block.render(region, buf);

        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(TAIL_LINES as u16),
                Constraint::Min(0),
                Constraint::Length(1),
            ])
            .split(inner);

        Paragraph::new(Line::from(Span::styled(
            format!(" {}", self.stage),
            Style::default().fg(theme::FG).bold(),
        )))
        .render(v[0], buf);

        Paragraph::new(self.progress_line()).render(v[2], buf);

        let tail_lines: Vec<Line> = self
            .tail
            .iter()
            .map(|l| {
                Line::from(Span::styled(
                    format!("  {}", l.message),
                    Style::default().fg(theme::MUTED),
                ))
            })
            .collect();
        Paragraph::new(tail_lines)
            .block(
                Block::default()
                    .borders(Borders::TOP | Borders::BOTTOM)
                    .border_style(Style::default().fg(theme::MUTED)),
            )
            .render(v[3], buf);

        let action = match (&self.phase, self.cancel_requested) {
            (Phase::Running, false) => "[c] cancel",
            (Phase::Running, true) => "cancelling at next stage boundary...",
            (Phase::Succeeded, _) => "[enter] view results",
            (Phase::Failed(_), _) => "[r] retry   [enter] dismiss",
        };
        Paragraph::new(Line::from(Span::styled(
            format!(" {action}"),
            Style::default().fg(theme::ACCENT).bold(),
        )))
        .render(v[5], buf);
    }

    fn progress_line(&self) -> Line<'static> {
        let width = 40usize;
        let (filled, label) = match (&self.phase, &self.progress) {
            (Phase::Succeeded, _) => (width, "100%".to_string()),
            (_, Some(snap)) if !snap.indeterminate && snap.length > 0 => {
                let ratio = (snap.position as f64 / snap.length as f64).clamp(0.0, 1.0);
                let f = (ratio * width as f64).round() as usize;
                (f, format!("{}/{}", snap.position, snap.length))
            }
            (_, Some(snap)) => (0, format!("{}", snap.position)),
            (_, None) => (0, String::new()),
        };
        let mut bar = String::with_capacity(width + 2);
        bar.push('[');
        for i in 0..width {
            if i < filled {
                bar.push('=');
            } else {
                bar.push(' ');
            }
        }
        bar.push(']');
        Line::from(vec![
            Span::raw(" "),
            Span::styled(bar, Style::default().fg(theme::ACCENT)),
            Span::raw(" "),
            Span::styled(label, Style::default().fg(theme::FG)),
        ])
    }
}

pub fn dim_area(area: Rect, buf: &mut Buffer) {
    let dim = Style::default().add_modifier(Modifier::DIM);
    for y in area.top()..area.bottom() {
        for x in area.left()..area.right() {
            buf[(x, y)].set_style(dim);
        }
    }
}

fn centered_rect(area: Rect, percent_x: u16, percent_y: u16) -> Rect {
    let v = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Percentage((100 - percent_y) / 2),
            Constraint::Percentage(percent_y),
            Constraint::Percentage((100 - percent_y) / 2),
        ])
        .split(area);
    Layout::default()
        .direction(Direction::Horizontal)
        .constraints([
            Constraint::Percentage((100 - percent_x) / 2),
            Constraint::Percentage(percent_x),
            Constraint::Percentage((100 - percent_x) / 2),
        ])
        .split(v[1])[1]
}
