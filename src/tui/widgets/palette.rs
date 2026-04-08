use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use nucleo_matcher::pattern::{CaseMatching, Normalization, Pattern};
use nucleo_matcher::{Config, Matcher, Utf32Str};
use ratatui::layout::Rect;
use ratatui::style::{Modifier, Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Clear, Padding, Paragraph};
use ratatui::Frame;

use crate::tui::action::{format_binding, Action, ActionScope, KeyMap};
use crate::tui::theme;

pub enum PaletteOutcome {
    Stay,
    Close,
    Run(Action),
}

pub struct Palette {
    query: String,
    matches: Vec<(Action, u32)>,
    cursor: usize,
    scope: ActionScope,
    matcher: Matcher,
}

enum Row {
    Header(ActionScope),
    Action(Action, usize),
}

impl Palette {
    pub fn open(scope: ActionScope) -> Self {
        let mut palette = Self {
            query: String::new(),
            matches: Vec::new(),
            cursor: 0,
            scope,
            matcher: Matcher::new(Config::DEFAULT),
        };
        palette.refresh();
        palette
    }

    pub fn handle_key(&mut self, key: KeyEvent) -> PaletteOutcome {
        match (key.code, key.modifiers) {
            (KeyCode::Esc, _) => PaletteOutcome::Close,
            (KeyCode::Enter, _) => match self.matches.get(self.cursor) {
                Some((action, _)) => PaletteOutcome::Run(*action),
                None => PaletteOutcome::Stay,
            },
            (KeyCode::Up, _) => {
                if self.cursor > 0 {
                    self.cursor -= 1;
                }
                PaletteOutcome::Stay
            }
            (KeyCode::Down, _) => {
                if self.cursor + 1 < self.matches.len() {
                    self.cursor += 1;
                }
                PaletteOutcome::Stay
            }
            (KeyCode::Char('n'), m) if m.contains(KeyModifiers::CONTROL) => {
                if self.cursor + 1 < self.matches.len() {
                    self.cursor += 1;
                }
                PaletteOutcome::Stay
            }
            (KeyCode::Char('p'), m) if m.contains(KeyModifiers::CONTROL) => {
                if self.cursor > 0 {
                    self.cursor -= 1;
                }
                PaletteOutcome::Stay
            }
            (KeyCode::Backspace, _) => {
                self.query.pop();
                self.refresh();
                PaletteOutcome::Stay
            }
            (KeyCode::Char(c), m) if !m.contains(KeyModifiers::CONTROL) && !m.contains(KeyModifiers::ALT) => {
                self.query.push(c);
                self.refresh();
                PaletteOutcome::Stay
            }
            _ => PaletteOutcome::Stay,
        }
    }

    fn refresh(&mut self) {
        let candidates: Vec<(Action, String)> = Action::all()
            .iter()
            .copied()
            .filter(|a| {
                let s = a.scope();
                s == ActionScope::Global || s == self.scope
            })
            .map(|a| (a, format!("{} {}", a.title(), a.description())))
            .collect();

        if self.query.is_empty() {
            self.matches = candidates.into_iter().map(|(a, _)| (a, 0)).collect();
        } else {
            let pattern = Pattern::parse(&self.query, CaseMatching::Smart, Normalization::Smart);
            let mut buf: Vec<char> = Vec::new();
            let mut scored: Vec<(Action, u32)> = candidates
                .into_iter()
                .filter_map(|(a, hay)| {
                    buf.clear();
                    let haystack = Utf32Str::new(&hay, &mut buf);
                    pattern.score(haystack, &mut self.matcher).map(|s| (a, s))
                })
                .collect();
            scored.sort_by(|a, b| b.1.cmp(&a.1));
            self.matches = scored;
        }

        if self.cursor >= self.matches.len() {
            self.cursor = self.matches.len().saturating_sub(1);
        }
    }

    fn rows(&self) -> Vec<Row> {
        let mut by_scope: Vec<(ActionScope, Vec<(Action, usize)>)> = Vec::new();
        for (idx, (action, _)) in self.matches.iter().enumerate() {
            let scope = action.scope();
            match by_scope.iter_mut().find(|(s, _)| *s == scope) {
                Some((_, v)) => v.push((*action, idx)),
                None => by_scope.push((scope, vec![(*action, idx)])),
            }
        }
        if self.query.is_empty() {
            by_scope.sort_by_key(|(s, _)| {
                ActionScope::ordered().iter().position(|x| x == s).unwrap_or(usize::MAX)
            });
        }
        let mut rows = Vec::with_capacity(self.matches.len() + by_scope.len());
        for (scope, items) in by_scope {
            rows.push(Row::Header(scope));
            for (a, i) in items {
                rows.push(Row::Action(a, i));
            }
        }
        rows
    }

    pub fn draw(&self, frame: &mut Frame, area: Rect) {
        let overlay = centered(area, 60, 60);
        frame.render_widget(Clear, overlay);

        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" command palette · {} ", self.scope.title().to_lowercase()))
            .border_style(Style::default().fg(theme::ACCENT))
            .padding(Padding::new(1, 1, 0, 0));
        let inner = block.inner(overlay);
        frame.render_widget(block, overlay);

        if inner.height < 3 || inner.width < 10 {
            return;
        }

        let prompt_area = Rect { x: inner.x, y: inner.y, width: inner.width, height: 1 };
        let body_height = inner.height.saturating_sub(2);
        let body_area = Rect {
            x: inner.x,
            y: inner.y + 1,
            width: inner.width,
            height: body_height,
        };
        let footer_area = Rect {
            x: inner.x,
            y: inner.y + inner.height - 1,
            width: inner.width,
            height: 1,
        };

        let prompt = Paragraph::new(Line::from(vec![
            Span::styled("> ", Style::default().fg(theme::WARN).bold()),
            Span::styled(self.query.as_str(), Style::default().fg(theme::FG)),
            Span::styled("_", Style::default().fg(theme::ACCENT)),
        ]));
        frame.render_widget(prompt, prompt_area);

        let rows = self.rows();
        let focus_row_in_rows = rows.iter().position(|r| matches!(r, Row::Action(_, i) if *i == self.cursor));
        let visible = body_area.height as usize;
        let start = match focus_row_in_rows {
            Some(pos) if pos >= visible => pos + 1 - visible,
            _ => 0,
        };

        let gutter_w: usize = 2;
        let key_w: usize = 12;
        let title_w: usize = 24;

        let lines: Vec<Line> = rows
            .iter()
            .skip(start)
            .take(visible)
            .map(|row| match row {
                Row::Header(scope) => Line::from(Span::styled(
                    format!("{:>w$}{}", "", scope.title().to_lowercase(), w = gutter_w),
                    Style::default().fg(theme::MUTED).add_modifier(Modifier::DIM),
                )),
                Row::Action(action, idx) => {
                    let is_sel = *idx == self.cursor;
                    let gutter = if is_sel {
                        Span::styled(
                            format!("{} ", theme::FOCUS_GLYPH),
                            Style::default().fg(theme::ACCENT).bold(),
                        )
                    } else {
                        Span::raw(format!("{:w$}", "", w = gutter_w))
                    };
                    let key_label = action
                        .default_key()
                        .map(|(c, m)| format_binding(c, m))
                        .unwrap_or_default();
                    let title_style = if is_sel {
                        Style::default().fg(theme::FG).bold()
                    } else {
                        Style::default().fg(theme::MUTED)
                    };
                    let desc_style = if is_sel {
                        Style::default().fg(theme::MUTED)
                    } else {
                        Style::default().fg(theme::MUTED).add_modifier(Modifier::DIM)
                    };
                    let key_style = if is_sel {
                        Style::default().fg(theme::WARN)
                    } else {
                        Style::default().fg(theme::MUTED).add_modifier(Modifier::DIM)
                    };
                    Line::from(vec![
                        gutter,
                        Span::styled(format!("{:<kw$}", key_label, kw = key_w), key_style),
                        Span::styled(format!("{:<tw$}", action.title(), tw = title_w), title_style),
                        Span::styled(action.description(), desc_style),
                    ])
                }
            })
            .collect();
        frame.render_widget(Paragraph::new(lines), body_area);

        let action_count = self.matches.len();
        let mut spans: Vec<Span> = vec![Span::styled(
            format!(
                "{} match{}    ",
                action_count,
                if action_count == 1 { "" } else { "es" }
            ),
            theme::hint_bar_style(),
        )];
        for (i, (code, mods, action)) in KeyMap::for_scope(ActionScope::Palette)
            .entries()
            .iter()
            .enumerate()
        {
            if i > 0 {
                spans.push(Span::styled("  ", theme::hint_bar_style()));
            }
            spans.push(Span::styled(
                format_binding(*code, *mods),
                Style::default().fg(theme::ACCENT).add_modifier(Modifier::BOLD),
            ));
            spans.push(Span::styled(format!(" {}", action.title()), theme::hint_bar_style()));
        }
        frame.render_widget(Paragraph::new(Line::from(spans)), footer_area);
    }
}

fn centered(area: Rect, pct_w: u16, pct_h: u16) -> Rect {
    let w = area.width * pct_w / 100;
    let h = area.height * pct_h / 100;
    let x = area.x + (area.width.saturating_sub(w)) / 2;
    let y = area.y + (area.height.saturating_sub(h)) / 2;
    Rect { x, y, width: w, height: h }
}
