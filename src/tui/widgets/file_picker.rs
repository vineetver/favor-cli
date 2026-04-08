use std::path::{Path, PathBuf};

use ratatui::layout::{Constraint, Direction, Layout};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, ListState, Padding, Paragraph};
use ratatui::Frame;

use crate::config::{DirProbe, ProbeStatus};
use crate::tui::theme;

pub const SELECT_ENTRY: &str = "[ Use this directory ]";

pub struct DirBrowserState {
    pub prompt: String,
    pub current_dir: PathBuf,
    pub entries: Vec<String>,
    pub list_state: ListState,
    pub typing_path: bool,
    pub input_buf: String,
    pub probe: DirProbe,
    pub show_files: bool,
}

impl DirBrowserState {
    pub fn new(prompt: &str, start: &Path) -> Self {
        Self::with_files(prompt, start, false)
    }

    pub fn with_files(prompt: &str, start: &Path, show_files: bool) -> Self {
        let current_dir = if start.is_dir() {
            start.to_path_buf()
        } else {
            start
                .parent()
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| PathBuf::from("/"))
        };
        let entries = list_entries(&current_dir, show_files);
        let mut list_state = ListState::default();
        if !entries.is_empty() {
            list_state.select(Some(0));
        }
        let probe = DirProbe::scan(&current_dir);

        Self {
            prompt: prompt.to_string(),
            current_dir,
            entries,
            list_state,
            typing_path: false,
            input_buf: String::new(),
            probe,
            show_files,
        }
    }

    pub fn navigate_to(&mut self, dir: PathBuf) {
        self.current_dir = dir;
        self.entries = list_entries(&self.current_dir, self.show_files);
        self.list_state
            .select(if self.entries.is_empty() { None } else { Some(0) });
        self.input_buf = self.current_dir.to_string_lossy().to_string();
        self.probe = DirProbe::scan(&self.current_dir);
    }

    pub fn go_parent(&mut self) {
        if let Some(parent) = self.current_dir.parent() {
            self.navigate_to(parent.to_path_buf());
        }
    }

    pub fn enter_selected(&mut self) -> Option<PathBuf> {
        let i = self.list_state.selected()?;
        let name = &self.entries[i];

        if name == SELECT_ENTRY {
            return Some(self.current_dir.clone());
        }
        if name == ".." {
            self.go_parent();
            return None;
        }

        let target = self.current_dir.join(name);
        if target.is_dir() {
            self.navigate_to(target);
            return None;
        }
        if self.show_files && target.is_file() {
            return Some(target);
        }
        None
    }

    pub fn select_up(&mut self) {
        if let Some(i) = self.list_state.selected() {
            if i > 0 {
                self.list_state.select(Some(i - 1));
            }
        }
    }

    pub fn select_down(&mut self) {
        if let Some(i) = self.list_state.selected() {
            if i + 1 < self.entries.len() {
                self.list_state.select(Some(i + 1));
            }
        }
    }
}

pub fn list_entries(dir: &Path, show_files: bool) -> Vec<String> {
    let mut entries = vec![SELECT_ENTRY.to_string(), "..".to_string()];
    let Ok(read) = std::fs::read_dir(dir) else {
        return entries;
    };
    let mut dirs: Vec<String> = Vec::new();
    let mut files: Vec<String> = Vec::new();
    for e in read.flatten() {
        let Ok(name) = e.file_name().into_string() else {
            continue;
        };
        let p = e.path();
        if p.is_dir() {
            dirs.push(name);
        } else if show_files && p.is_file() {
            files.push(name);
        }
    }
    dirs.sort();
    files.sort();
    entries.extend(dirs);
    entries.extend(files);
    entries
}

pub fn tab_complete(partial: &str) -> Option<String> {
    let path = Path::new(partial);
    let (parent, prefix) = if partial.ends_with('/') {
        if path.is_dir() {
            return list_first_child(path);
        }
        return None;
    } else {
        let parent = path.parent()?;
        let prefix = path.file_name()?.to_string_lossy().to_string();
        (parent, prefix)
    };

    if !parent.is_dir() {
        return None;
    }

    let matches: Vec<String> = std::fs::read_dir(parent)
        .ok()?
        .flatten()
        .filter(|e| e.path().is_dir())
        .filter_map(|e| {
            let name = e.file_name().into_string().ok()?;
            if name.starts_with(&prefix) {
                Some(name)
            } else {
                None
            }
        })
        .collect();

    match matches.len() {
        0 => None,
        1 => {
            let completed = parent.join(&matches[0]);
            Some(format!("{}/", completed.to_string_lossy()))
        }
        _ => {
            let lcp = longest_common_prefix(&matches);
            let completed = parent.join(&lcp);
            Some(completed.to_string_lossy().to_string())
        }
    }
}

fn list_first_child(dir: &Path) -> Option<String> {
    let mut entries: Vec<String> = std::fs::read_dir(dir)
        .ok()?
        .flatten()
        .filter(|e| e.path().is_dir())
        .filter_map(|e| e.file_name().into_string().ok())
        .collect();
    entries.sort();
    if entries.len() == 1 {
        Some(format!("{}/{}/", dir.to_string_lossy(), entries[0]))
    } else {
        None
    }
}

fn longest_common_prefix(strings: &[String]) -> String {
    if strings.is_empty() {
        return String::new();
    }
    let first = &strings[0];
    let mut len = first.len();
    for s in &strings[1..] {
        len = len.min(s.len());
        for (i, (a, b)) in first.bytes().zip(s.bytes()).enumerate() {
            if a != b {
                len = len.min(i);
                break;
            }
        }
    }
    first[..len].to_string()
}

pub fn draw(frame: &mut Frame, area: ratatui::layout::Rect, state: &mut DirBrowserState) {
    let h_split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(area);

    let left_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Min(6),
            Constraint::Length(2),
        ])
        .split(h_split[0]);

    let title = Paragraph::new(format!("  {}", state.prompt))
        .style(Style::default().fg(theme::ACCENT).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, left_layout[0]);

    if state.typing_path {
        let path_valid = Path::new(&state.input_buf).is_dir();
        let path_color = if path_valid { theme::OK } else { theme::BAD };
        let path_line = Line::from(vec![
            Span::styled(" > ", Style::default().fg(theme::WARN)),
            Span::styled(&state.input_buf, Style::default().fg(path_color)),
            Span::styled("_", Style::default().fg(theme::ACCENT)),
        ]);
        let hint = if path_valid {
            " Valid directory — enter to go "
        } else {
            " Not a directory "
        };
        let path_block = Paragraph::new(path_line).block(
            Block::default()
                .borders(Borders::ALL)
                .title(hint)
                .border_style(Style::default().fg(if path_valid { theme::OK } else { theme::WARN })),
        );
        frame.render_widget(path_block, left_layout[1]);
    } else {
        let path_line = Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled(
                state.current_dir.to_string_lossy().to_string(),
                Style::default().fg(theme::FG).bold(),
            ),
        ]);
        let path_block = Paragraph::new(path_line).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Current directory ")
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(path_block, left_layout[1]);
    }

    let items: Vec<ListItem> = state
        .entries
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let is_selected = state.list_state.selected() == Some(i);
            let (prefix, style) = match name.as_str() {
                s if s == SELECT_ENTRY => (
                    " ",
                    if is_selected {
                        Style::default().fg(theme::OK).bold()
                    } else {
                        Style::default().fg(theme::OK)
                    },
                ),
                ".." => (
                    " ..",
                    if is_selected {
                        Style::default().fg(theme::ACCENT).bold()
                    } else {
                        Style::default().fg(theme::MUTED)
                    },
                ),
                s if s.starts_with('.') => (
                    &name[..],
                    if is_selected {
                        Style::default().fg(theme::MUTED).bold()
                    } else {
                        Style::default().fg(theme::MUTED)
                    },
                ),
                _ => (
                    &name[..],
                    if is_selected {
                        Style::default().fg(theme::ACCENT).bold()
                    } else {
                        Style::default().fg(theme::FG)
                    },
                ),
            };
            let display = if name == SELECT_ENTRY || name == ".." {
                format!(" {prefix}")
            } else {
                format!(" /{name}")
            };
            ListItem::new(display).style(style)
        })
        .collect();

    let list = List::new(items)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Directories ")
                .border_style(Style::default().fg(theme::MUTED)),
        )
        .highlight_style(Style::default().bg(theme::MUTED).fg(theme::FG))
        .highlight_symbol(" > ");
    frame.render_stateful_widget(list, left_layout[2], &mut state.list_state);

    let help_text = if state.typing_path {
        "  type path    enter go    esc cancel"
    } else {
        "  enter open    space select    / type path    esc cancel"
    };
    let help = Paragraph::new(help_text).style(Style::default().fg(theme::MUTED));
    frame.render_widget(help, left_layout[3]);

    let right_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(6),
            Constraint::Length(3),
        ])
        .split(h_split[1]);

    let probe_header = if state.probe.has_any_data() {
        Paragraph::new("  FAVOR data detected")
            .style(Style::default().fg(theme::OK).bold())
            .block(Block::default().padding(Padding::top(1)))
    } else {
        Paragraph::new("  No FAVOR data at this path")
            .style(Style::default().fg(theme::MUTED))
            .block(Block::default().padding(Padding::top(1)))
    };
    frame.render_widget(probe_header, right_layout[0]);

    let summary = state.probe.summary_lines();
    let probe_lines: Vec<Line> = summary
        .iter()
        .map(|(text, status)| {
            let (icon, color) = match status {
                ProbeStatus::Good => ("  +  ", theme::OK),
                ProbeStatus::Partial => ("  ~  ", theme::WARN),
                ProbeStatus::Missing => ("  -  ", theme::MUTED),
                ProbeStatus::Info => ("     ", theme::MUTED),
            };
            Line::from(vec![
                Span::styled(icon, Style::default().fg(color)),
                Span::styled(text, Style::default().fg(color)),
            ])
        })
        .collect();

    let probe_detail = Paragraph::new(probe_lines).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Data scan ")
            .border_style(if state.probe.has_any_data() {
                Style::default().fg(theme::OK)
            } else {
                Style::default().fg(theme::MUTED)
            }),
    );
    frame.render_widget(probe_detail, right_layout[1]);

    let verdict = if state.probe.has_any_data() {
        let tier = state
            .probe
            .detected_tier()
            .map(|t| t.as_str())
            .unwrap_or("?");
        Paragraph::new(Line::from(vec![
            Span::styled("  Ready — ", Style::default().fg(theme::OK)),
            Span::styled(
                format!("{tier} tier detected"),
                Style::default().fg(theme::OK).bold(),
            ),
        ]))
    } else {
        Paragraph::new(Line::from(vec![Span::styled(
            "  Will download data after setup",
            Style::default().fg(theme::WARN),
        )]))
    };
    frame.render_widget(verdict, right_layout[2]);
}
