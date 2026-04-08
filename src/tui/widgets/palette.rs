use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use nucleo_matcher::pattern::{CaseMatching, Normalization, Pattern};
use nucleo_matcher::{Config, Matcher, Utf32Str};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Clear, List, ListItem, Padding, Paragraph};
use ratatui::Frame;

use crate::tui::action::{format_binding, Action, ActionScope};
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
            self.matches = candidates
                .into_iter()
                .map(|(a, _)| (a, 0))
                .collect();
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

    pub fn draw(&self, frame: &mut Frame, area: Rect) {
        let overlay = centered(area, 60, 60);
        frame.render_widget(Clear, overlay);

        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" Command palette — {} ", self.scope.title()))
            .border_style(Style::default().fg(theme::ACCENT))
            .padding(Padding::new(1, 1, 0, 0));
        frame.render_widget(block, overlay);

        let inner = Rect {
            x: overlay.x + 2,
            y: overlay.y + 1,
            width: overlay.width.saturating_sub(4),
            height: overlay.height.saturating_sub(2),
        };

        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Min(3),
                Constraint::Length(1),
            ])
            .split(inner);

        let prompt = Paragraph::new(Line::from(vec![
            Span::styled("> ", Style::default().fg(theme::WARN).bold()),
            Span::styled(self.query.as_str(), Style::default().fg(theme::FG)),
            Span::styled("_", Style::default().fg(theme::ACCENT)),
        ]));
        frame.render_widget(prompt, layout[0]);

        let divider = Paragraph::new(Span::styled(
            "-".repeat(layout[1].width as usize),
            Style::default().fg(theme::MUTED),
        ));
        frame.render_widget(divider, layout[1]);

        let visible_height = layout[2].height as usize;
        let start = self.cursor.saturating_sub(visible_height.saturating_sub(1));
        let items: Vec<ListItem> = self
            .matches
            .iter()
            .enumerate()
            .skip(start)
            .take(visible_height)
            .map(|(i, (action, _))| {
                let is_sel = i == self.cursor;
                let marker = if is_sel { " > " } else { "   " };
                let title_style = if is_sel {
                    Style::default().fg(theme::ACCENT).bold()
                } else {
                    Style::default().fg(theme::FG)
                };
                let key_label = action
                    .default_key()
                    .map(|(c, m)| format_binding(c, m))
                    .unwrap_or_default();
                ListItem::new(Line::from(vec![
                    Span::raw(marker),
                    Span::styled(format!("{:<10}", key_label), Style::default().fg(theme::WARN)),
                    Span::styled(format!("{:<28}", action.title()), title_style),
                    Span::styled(action.description(), Style::default().fg(theme::MUTED)),
                ]))
            })
            .collect();
        let list = List::new(items);
        frame.render_widget(list, layout[2]);

        let footer = Paragraph::new(Line::from(vec![
            Span::styled(
                format!("{} match", self.matches.len()),
                Style::default().fg(theme::MUTED),
            ),
            Span::styled(
                if self.matches.len() == 1 { "" } else { "es" },
                Style::default().fg(theme::MUTED),
            ),
            Span::styled("    enter run    esc close", Style::default().fg(theme::MUTED)),
        ]));
        frame.render_widget(footer, layout[3]);
    }
}

fn centered(area: Rect, pct_w: u16, pct_h: u16) -> Rect {
    let w = area.width * pct_w / 100;
    let h = area.height * pct_h / 100;
    let x = area.x + (area.width.saturating_sub(w)) / 2;
    let y = area.y + (area.height.saturating_sub(h)) / 2;
    Rect { x, y, width: w, height: h }
}
