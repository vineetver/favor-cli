use std::collections::VecDeque;

use ratatui::layout::Rect;
use ratatui::style::{Color, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph};
use ratatui::Frame;

use crate::tui::output::{LogLevel, LogLine, ProgressSnapshot};
use crate::tui::theme;

const SPINNER_GLYPHS: &[char] = &['|', '/', '-', '\\'];

enum Entry {
    Log(LogLine),
    Bar(ProgressSnapshot),
}

pub struct LogTail {
    entries: VecDeque<Entry>,
    cap: usize,
    spinner_tick: usize,
}

impl LogTail {
    pub fn new(cap: usize) -> Self {
        Self {
            entries: VecDeque::with_capacity(cap),
            cap,
            spinner_tick: 0,
        }
    }

    pub fn push_log(&mut self, line: LogLine) {
        self.push_entry(Entry::Log(line));
    }

    pub fn push_snapshot(&mut self, snap: ProgressSnapshot) {
        self.spinner_tick = self.spinner_tick.wrapping_add(1);
        for entry in self.entries.iter_mut() {
            if let Entry::Bar(existing) = entry {
                if existing.id == snap.id {
                    *existing = snap;
                    return;
                }
            }
        }
        self.push_entry(Entry::Bar(snap));
    }

    pub fn sweep(&mut self, id: u64) {
        self.entries.retain(|e| match e {
            Entry::Bar(snap) => snap.id != id,
            Entry::Log(_) => true,
        });
    }

    fn push_entry(&mut self, entry: Entry) {
        if self.entries.len() == self.cap {
            self.entries.pop_front();
        }
        self.entries.push_back(entry);
    }

    pub fn draw(&self, frame: &mut Frame, area: Rect, title: &str) {
        let lines: Vec<Line> = self
            .entries
            .iter()
            .map(|e| match e {
                Entry::Log(l) => {
                    let color = match l.level {
                        LogLevel::Status => theme::FG,
                        LogLevel::Success => theme::OK,
                        LogLevel::Warn => theme::WARN,
                        LogLevel::Error => theme::BAD,
                    };
                    Line::from(Span::styled(
                        format!("  {}", l.message),
                        Style::default().fg(color),
                    ))
                }
                Entry::Bar(snap) => self.render_bar(snap),
            })
            .collect();
        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" {title} "))
            .border_style(Style::default().fg(Color::DarkGray));
        frame.render_widget(Paragraph::new(lines).block(block), area);
    }

    fn render_bar(&self, snap: &ProgressSnapshot) -> Line<'static> {
        let label_color = if snap.finished { theme::OK } else { theme::FG };
        let label = Span::styled(
            format!("  > {:<24} ", snap.message),
            Style::default().fg(label_color),
        );
        if snap.indeterminate {
            let glyph = SPINNER_GLYPHS[self.spinner_tick % SPINNER_GLYPHS.len()];
            let body = format!("[{glyph}]  {}", snap.position);
            Line::from(vec![label, Span::styled(body, Style::default().fg(theme::ACCENT))])
        } else {
            let width = 20usize;
            let ratio = if snap.length == 0 {
                0.0
            } else {
                (snap.position as f64 / snap.length as f64).clamp(0.0, 1.0)
            };
            let filled = (ratio * width as f64).round() as usize;
            let mut bar = String::with_capacity(width + 2);
            bar.push('[');
            for i in 0..width {
                if i + 1 < filled {
                    bar.push('=');
                } else if i < filled {
                    bar.push('>');
                } else {
                    bar.push(' ');
                }
            }
            bar.push(']');
            let body = format!("{bar} {}/{}", snap.position, snap.length);
            Line::from(vec![label, Span::styled(body, Style::default().fg(theme::ACCENT))])
        }
    }
}
