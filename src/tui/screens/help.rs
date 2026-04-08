use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph};
use ratatui::Frame;

use crate::tui::action::{format_binding, Action, ActionScope, KeyMap};
use crate::tui::screen::{Screen, Transition};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;

pub struct HelpScreen {
    scopes: Vec<ActionScope>,
    focused: usize,
    error: Option<String>,
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
        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" Help — {} ", focused_scope.title()))
            .border_style(Style::default().fg(theme::ACCENT));
        let inner = block.inner(area);
        frame.render_widget(block, area);

        if inner.height < 4 || inner.width < 20 {
            let msg = Paragraph::new(Line::from(Span::styled(
                "terminal too small",
                Style::default().fg(theme::WARN),
            )));
            frame.render_widget(msg, inner);
            return;
        }

        let chunks = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(1),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(inner);

        let mut lines: Vec<Line> = Vec::new();
        for (i, scope) in self.scopes.iter().enumerate() {
            let is_focused = i == self.focused;
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
        frame.render_widget(Paragraph::new(lines), chunks[0]);

        let error_line = match &self.error {
            Some(msg) => Line::from(Span::styled(msg.clone(), Style::default().fg(theme::WARN))),
            None => Line::from(Span::raw("")),
        };
        frame.render_widget(Paragraph::new(error_line), chunks[1]);

        let hint = Line::from(vec![
            Span::styled("Tab", Style::default().fg(theme::ACCENT).bold()),
            Span::styled(" cycle scope   ", Style::default().fg(theme::MUTED)),
            Span::styled("Esc", Style::default().fg(theme::ACCENT).bold()),
            Span::styled(" close", Style::default().fg(theme::MUTED)),
        ]);
        frame.render_widget(Paragraph::new(hint), chunks[2]);
    }
}
