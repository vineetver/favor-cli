use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::Rect;
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Paragraph, Widget};
use ratatui::Frame;

use crate::tui::action::{format_binding, Action, ActionScope, KeyMap};
use crate::tui::screen::{Screen, Transition};
use crate::tui::shell::{Binding, ErrorMessage, ScreenChrome, Shell};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;

pub struct HelpScreen {
    scopes: Vec<ActionScope>,
    focused: usize,
    error: Option<ErrorMessage>,
}

impl HelpScreen {
    pub fn new(origin: ActionScope) -> Self {
        let origin = if origin == ActionScope::Help {
            ActionScope::Global
        } else {
            origin
        };
        let mut scopes: Vec<ActionScope> = Vec::with_capacity(ActionScope::ordered().len());
        scopes.push(origin);
        for s in ActionScope::ordered() {
            if *s != origin && Action::all().iter().any(|a| a.scope() == *s) {
                scopes.push(*s);
            }
        }
        Self {
            scopes,
            focused: 0,
            error: None,
        }
    }

    fn cycle(&mut self) {
        if !self.scopes.is_empty() {
            self.focused = (self.focused + 1) % self.scopes.len();
        }
    }
}

impl Screen for HelpScreen {
    fn title(&self) -> &str {
        "Help"
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Help
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        KeyMap::new()
            .bind(KeyCode::Esc, none, Action::HelpClose)
            .bind(KeyCode::Char('q'), none, Action::HelpClose)
            .bind(KeyCode::Char('?'), none, Action::HelpClose)
            .bind(KeyCode::Tab, none, Action::HelpCycleScope)
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::HelpClose => Transition::Pop,
            Action::HelpCycleScope => {
                self.cycle();
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, _log: &LogTail) {
        let focused_scope = self.scopes.get(self.focused).copied().unwrap_or(ActionScope::Global);
        let title = format!("Help — {}", focused_scope.title());
        let hint = [
            Binding::new(
                (KeyCode::Tab, KeyModifiers::NONE),
                "cycle scope",
            ),
            Binding::new((KeyCode::Esc, KeyModifiers::NONE), "close"),
        ];
        let scopes = self.scopes.clone();
        let focused = self.focused;
        let chrome = ScreenChrome {
            title: &title,
            status: None,
            error: self.error.as_ref(),
            hint: &hint,
            graph: None,
        };
        let body = |inner: Rect, buf: &mut ratatui::buffer::Buffer| {
            let mut lines: Vec<Line> = Vec::new();
            for (i, scope) in scopes.iter().enumerate() {
                let is_focused = i == focused;
                let glyph = if is_focused { "▸ " } else { "  " };
                let header_style = if is_focused {
                    Style::default().fg(theme::ACCENT).bold()
                } else {
                    Style::default().fg(theme::MUTED).bold()
                };
                lines.push(Line::from(vec![
                    Span::styled(glyph, header_style),
                    Span::styled(scope.title().to_string(), header_style),
                ]));
                if !is_focused {
                    continue;
                }
                for action in Action::all().iter().filter(|a| a.scope() == *scope) {
                    let key = action
                        .default_key()
                        .map(|(c, m)| format_binding(c, m))
                        .unwrap_or_else(|| "—".to_string());
                    lines.push(Line::from(vec![
                        Span::styled("    ", Style::default()),
                        Span::styled(format!("{key:<14}"), Style::default().fg(theme::WARN)),
                        Span::styled(
                            format!("{:<26}", action.title()),
                            Style::default().fg(theme::FG),
                        ),
                        Span::styled(action.description(), Style::default().fg(theme::MUTED)),
                    ]));
                }
            }
            Paragraph::new(lines).render(inner, buf);
        };
        Shell::new(chrome, body).render(area, frame.buffer_mut());
    }
}
