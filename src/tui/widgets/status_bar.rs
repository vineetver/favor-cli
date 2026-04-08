use ratatui::layout::Rect;
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::tui::theme;

pub struct StatusBar<'a> {
    pub title: &'a str,
    pub keys: &'a str,
}

impl<'a> StatusBar<'a> {
    pub fn render(&self, frame: &mut Frame, area: Rect) {
        let line = Line::from(vec![
            Span::styled(format!(" {} ", self.title), Style::default().fg(theme::ACCENT)),
            Span::styled("·", Style::default().fg(theme::MUTED)),
            Span::styled(format!(" {} ", self.keys), Style::default().fg(theme::MUTED)),
        ]);
        frame.render_widget(Paragraph::new(line), area);
    }
}
