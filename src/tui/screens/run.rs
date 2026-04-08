use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

pub enum RunState {
    Running,
    Ok,
    Err(String),
}

pub struct RunScreen {
    title: String,
    state: RunState,
    cancel_requested: bool,
}

impl RunScreen {
    pub fn new(title: String) -> Self {
        Self {
            title,
            state: RunState::Running,
            cancel_requested: false,
        }
    }
}

impl Screen for RunScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Run
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        match self.state {
            RunState::Running => KeyMap::new().bind(KeyCode::Char('c'), none, Action::RunCancelRequest),
            RunState::Ok | RunState::Err(_) => KeyMap::new()
                .bind(KeyCode::Enter, none, Action::RunReturn)
                .bind(KeyCode::Esc, none, Action::RunReturn),
        }
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::RunCancelRequest => {
                self.cancel_requested = true;
                Transition::Stay
            }
            Action::RunReturn => Transition::Pop,
            _ => Transition::Stay,
        }
    }

    fn on_other_event(&mut self, event: &AppEvent) -> Transition {
        match event {
            AppEvent::CommandDone(res) => {
                self.state = match res {
                    Ok(()) => RunState::Ok,
                    Err(e) => RunState::Err(e.to_string()),
                };
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail) {
        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Min(4),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let (status_label, status_color) = match &self.state {
            RunState::Running => ("running", theme::ACCENT),
            RunState::Ok => ("done", theme::OK),
            RunState::Err(_) => ("failed", theme::BAD),
        };
        let header = Paragraph::new(Line::from(vec![
            Span::styled(format!("  {}  ", self.title), Style::default().fg(theme::FG).bold()),
            Span::styled(format!("[{status_label}]"), Style::default().fg(status_color).bold()),
        ]));
        frame.render_widget(header, v[0]);

        log.draw(frame, v[1], "Output");

        let err_line = match &self.state {
            RunState::Err(msg) => Line::from(Span::styled(
                format!("  {msg}"),
                Style::default().fg(theme::BAD),
            )),
            _ => Line::from(""),
        };
        frame.render_widget(Paragraph::new(err_line), v[2]);

        let keys = match &self.state {
            RunState::Running if self.cancel_requested => {
                "cancel detaches view on completion (commands not interruptible)"
            }
            RunState::Running => "c cancel (advisory)  ctrl-c quit",
            RunState::Ok => "enter return to workspace",
            RunState::Err(_) => "enter return to workspace",
        };
        StatusBar {
            title: &self.title,
            keys,
        }
        .render(frame, v[3]);
    }
}
