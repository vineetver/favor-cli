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
    active_scope: ActionScope,
}

impl HelpScreen {
    pub fn new(active_scope: ActionScope) -> Self {
        Self { active_scope }
    }

    fn ordered_scopes(&self) -> Vec<ActionScope> {
        let mut order: Vec<ActionScope> = Vec::with_capacity(8);
        order.push(ActionScope::Global);
        if self.active_scope != ActionScope::Global {
            order.push(self.active_scope);
        }
        for s in ActionScope::ordered() {
            if !order.contains(s) {
                order.push(*s);
            }
        }
        order
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
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::HelpClose => Transition::Pop,
            _ => Transition::Stay,
        }
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, _log: &LogTail) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" Help — {} active ", self.active_scope.title()))
            .border_style(Style::default().fg(theme::ACCENT));
        let inner = block.inner(area);
        frame.render_widget(block, area);

        let scopes = self.ordered_scopes();
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints(vec![Constraint::Min(1); scopes.len()])
            .split(inner);

        for (i, scope) in scopes.iter().enumerate() {
            let actions: Vec<&Action> = Action::all()
                .iter()
                .filter(|a| a.scope() == *scope)
                .collect();
            let header_style = if *scope == self.active_scope {
                Style::default().fg(theme::ACCENT).bold()
            } else {
                Style::default().fg(theme::WARN).bold()
            };
            let mut lines: Vec<Line> = Vec::with_capacity(actions.len() + 1);
            lines.push(Line::from(Span::styled(
                format!("  {}", scope.title()),
                header_style,
            )));
            for action in actions {
                let key = action
                    .default_key()
                    .map(|(c, m)| format_binding(c, m))
                    .unwrap_or_else(|| "—".to_string());
                lines.push(Line::from(vec![
                    Span::styled(format!("    {key:<14}"), Style::default().fg(theme::WARN)),
                    Span::styled(
                        format!("{:<26}", action.title()),
                        Style::default().fg(theme::FG),
                    ),
                    Span::styled(action.description(), Style::default().fg(theme::MUTED)),
                ]));
            }
            frame.render_widget(Paragraph::new(lines), layout[i]);
        }
    }
}
