use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Paragraph, Widget};

use crate::tui::action::{format_binding, KeyBinding};
use crate::tui::theme;

#[derive(Clone, Debug)]
pub struct ErrorMessage {
    pub text: String,
}

impl ErrorMessage {
    #[allow(dead_code)]
    pub fn new(text: impl Into<String>) -> Self {
        Self { text: text.into() }
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Binding {
    pub key: KeyBinding,
    pub label: &'static str,
}

impl Binding {
    pub const fn new(key: KeyBinding, label: &'static str) -> Self {
        Self { key, label }
    }
}

pub struct ScreenChrome<'a> {
    pub title: &'a str,
    pub status: Option<&'a str>,
    pub error: Option<&'a ErrorMessage>,
    pub hint: &'a [Binding],
}

pub struct Shell<'a, F: FnOnce(Rect, &mut Buffer)> {
    pub chrome: ScreenChrome<'a>,
    pub body: F,
}

impl<'a, F: FnOnce(Rect, &mut Buffer)> Shell<'a, F> {
    pub fn new(chrome: ScreenChrome<'a>, body: F) -> Self {
        Self { chrome, body }
    }
}

impl<'a, F: FnOnce(Rect, &mut Buffer)> Widget for Shell<'a, F> {
    fn render(self, area: Rect, buf: &mut Buffer) {
        if area.height < 4 || area.width < 8 {
            Paragraph::new(Line::from(Span::styled(
                "terminal too small",
                Style::default().fg(theme::WARN),
            )))
            .render(area, buf);
            return;
        }

        let rows = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Min(1),
                Constraint::Length(1),
            ])
            .split(area);

        render_chrome_row(rows[0], buf, self.chrome.title, self.chrome.status);
        render_error_slot(rows[1], buf, self.chrome.error);
        (self.body)(rows[2], buf);
        render_hint_bar(rows[3], buf, self.chrome.hint);
    }
}

fn render_chrome_row(area: Rect, buf: &mut Buffer, title: &str, status: Option<&str>) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Min(1), Constraint::Length(status_width(status))])
        .split(area);

    let title_line = Line::from(vec![
        Span::styled(
            format!(" {} ", title),
            Style::default().fg(theme::ACCENT).add_modifier(Modifier::BOLD),
        ),
    ]);
    Paragraph::new(title_line).render(cols[0], buf);

    if let Some(s) = status {
        let status_line = Line::from(Span::styled(
            format!("{s} "),
            Style::default().fg(theme::MUTED),
        ));
        Paragraph::new(status_line)
            .alignment(ratatui::layout::Alignment::Right)
            .render(cols[1], buf);
    }
}

fn status_width(status: Option<&str>) -> u16 {
    status.map(|s| s.chars().count() as u16 + 1).unwrap_or(0)
}

fn render_error_slot(area: Rect, buf: &mut Buffer, error: Option<&ErrorMessage>) {
    let line = match error {
        Some(e) => Line::from(Span::styled(format!(" {} ", e.text), theme::error_slot_style())),
        None => Line::from(Span::raw("")),
    };
    Paragraph::new(line).render(area, buf);
}

fn render_hint_bar(area: Rect, buf: &mut Buffer, bindings: &[Binding]) {
    let mut spans: Vec<Span> = Vec::with_capacity(bindings.len() * 3);
    for (i, b) in bindings.iter().enumerate() {
        if i > 0 {
            spans.push(Span::styled("  ", theme::hint_bar_style()));
        }
        spans.push(Span::styled(
            format_binding(b.key.0, b.key.1),
            Style::default().fg(theme::ACCENT).add_modifier(Modifier::BOLD),
        ));
        spans.push(Span::styled(format!(" {}", b.label), theme::hint_bar_style()));
    }
    Paragraph::new(Line::from(spans)).render(area, buf);
}
