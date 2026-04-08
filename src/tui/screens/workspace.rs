use std::path::PathBuf;

use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, ListState, Paragraph};
use ratatui::Frame;

use crate::config::Config;
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::screen::{Screen, Transition};
use crate::tui::screens::setup::SetupScreen;
use crate::tui::screens::transform::TransformScreen;
use crate::tui::screens::variant::VariantScreen;
use crate::tui::state::artifacts::ArtifactKind;
use crate::tui::state::workspace::WorkspaceState;
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

pub struct WorkspaceScreen {
    title: String,
    state: WorkspaceState,
    list_state: ListState,
    notice: Option<String>,
}

impl WorkspaceScreen {
    pub fn new(cwd: PathBuf) -> Self {
        let mut extra_roots: Vec<PathBuf> = Vec::new();
        if let Ok(cfg) = Config::load() {
            let r = cfg.root_dir();
            if !r.as_os_str().is_empty() && r.exists() && r != cwd {
                extra_roots.push(r);
            }
        }
        let state = WorkspaceState::new(cwd, extra_roots);
        let mut list_state = ListState::default();
        if !state.artifacts.is_empty() {
            list_state.select(Some(0));
        }
        Self {
            title: "Workspace".to_string(),
            state,
            list_state,
            notice: None,
        }
    }

    fn sync_focus(&mut self) {
        self.list_state.select(if self.state.artifacts.is_empty() {
            None
        } else {
            Some(self.state.focus)
        });
    }
}

fn actions_for(kind: &ArtifactKind) -> &'static [&'static str] {
    match kind {
        ArtifactKind::RawVcf => &["ingest", "browse"],
        ArtifactKind::PhenotypeTsv => &["use as phenotype"],
        ArtifactKind::KinshipTsv => &["use as kinship"],
        ArtifactKind::ParquetFile => &["browse"],
        ArtifactKind::IngestedSet => &["annotate", "browse"],
        ArtifactKind::AnnotatedSet { .. } => &["staar", "enrich", "browse"],
        ArtifactKind::GenotypeStore => &["inspect"],
        ArtifactKind::StaarResults => &["browse"],
        ArtifactKind::AnnotationRoot => &["status"],
    }
}

fn fmt_size(n: u64) -> String {
    const KB: f64 = 1024.0;
    const MB: f64 = KB * 1024.0;
    const GB: f64 = MB * 1024.0;
    let n = n as f64;
    if n >= GB {
        format!("{:.1} GB", n / GB)
    } else if n >= MB {
        format!("{:.1} MB", n / MB)
    } else if n >= KB {
        format!("{:.1} KB", n / KB)
    } else {
        format!("{n} B")
    }
}

impl Screen for WorkspaceScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn on_focus(&mut self) {
        self.state.rescan();
        self.sync_focus();
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail) {
        if self.state.drain_scan() {
            self.sync_focus();
        }

        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(4),
                Constraint::Length(6),
                Constraint::Length(1),
            ])
            .split(area);

        let h = Layout::default()
            .direction(Direction::Horizontal)
            .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
            .split(v[0]);

        let items: Vec<ListItem> = self
            .state
            .artifacts
            .iter()
            .map(|a| {
                let line = Line::from(vec![
                    Span::styled(
                        format!(" {} ", a.kind.glyph()),
                        Style::default().fg(theme::WARN),
                    ),
                    Span::styled(
                        format!("{:<28}", a.display_name),
                        Style::default().fg(theme::FG),
                    ),
                    Span::styled(
                        a.path.to_string_lossy().into_owned(),
                        Style::default().fg(theme::MUTED),
                    ),
                ]);
                ListItem::new(line)
            })
            .collect();

        let list_title = if self.state.scanning {
            format!(
                " Artifacts ({}, scanning…) ",
                self.state.artifacts.len()
            )
        } else {
            format!(" Artifacts ({}) ", self.state.artifacts.len())
        };
        let list = List::new(items)
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title(list_title)
                    .border_style(Style::default().fg(theme::ACCENT)),
            )
            .highlight_style(Style::default().bg(theme::MUTED).fg(theme::FG))
            .highlight_symbol(" > ");
        frame.render_stateful_widget(list, h[0], &mut self.list_state);

        let detail_lines: Vec<Line> = match self.state.focused() {
            Some(a) => {
                let mut lines = vec![
                    Line::from(Span::styled(
                        format!("  {}", a.kind.title()),
                        Style::default().fg(theme::ACCENT).bold(),
                    )),
                    Line::from(""),
                    Line::from(vec![
                        Span::styled("  path  ", Style::default().fg(theme::MUTED)),
                        Span::styled(
                            a.path.to_string_lossy().into_owned(),
                            Style::default().fg(theme::FG),
                        ),
                    ]),
                ];
                if a.size_bytes > 0 {
                    lines.push(Line::from(vec![
                        Span::styled("  size  ", Style::default().fg(theme::MUTED)),
                        Span::styled(fmt_size(a.size_bytes), Style::default().fg(theme::FG)),
                    ]));
                }
                if let ArtifactKind::AnnotatedSet { tier } = &a.kind {
                    lines.push(Line::from(vec![
                        Span::styled("  tier  ", Style::default().fg(theme::MUTED)),
                        Span::styled(tier.as_str(), Style::default().fg(theme::FG)),
                    ]));
                }
                lines.push(Line::from(""));
                lines.push(Line::from(Span::styled(
                    "  actions",
                    Style::default().fg(theme::ACCENT),
                )));
                for action in actions_for(&a.kind) {
                    lines.push(Line::from(Span::styled(
                        format!("    {action}"),
                        Style::default().fg(theme::FG),
                    )));
                }
                lines
            }
            None => {
                let msg = if self.state.scanning {
                    "  Scanning current directory…"
                } else {
                    "  No artifacts found. Press r to rescan."
                };
                vec![Line::from(Span::styled(
                    msg,
                    Style::default().fg(theme::MUTED),
                ))]
            }
        };

        let detail = Paragraph::new(detail_lines).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Detail ")
                .border_style(Style::default().fg(theme::MUTED)),
        );
        frame.render_widget(detail, h[1]);

        log.draw(frame, v[1], "Log");

        let status_keys = match &self.notice {
            Some(msg) => msg.as_str(),
            None => "q quit  j/k move  enter transform  r rescan  s setup  ? help",
        };
        StatusBar {
            title: &self.title,
            keys: status_keys,
        }
        .render(frame, v[2]);
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Workspace
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        KeyMap::new()
            .bind(KeyCode::Char('q'), none, Action::Quit)
            .bind(KeyCode::Esc, none, Action::Quit)
            .bind(KeyCode::Down, none, Action::WorkspaceDown)
            .bind(KeyCode::Char('j'), none, Action::WorkspaceDown)
            .bind(KeyCode::Up, none, Action::WorkspaceUp)
            .bind(KeyCode::Char('k'), none, Action::WorkspaceUp)
            .bind(KeyCode::Char('r'), none, Action::WorkspaceRescan)
            .bind(KeyCode::Enter, none, Action::WorkspaceOpenFocused)
            .bind(KeyCode::Char('s'), none, Action::WorkspaceOpenSetup)
    }

    fn on_action(&mut self, action: Action) -> Transition {
        match action {
            Action::Quit => Transition::Quit,
            Action::WorkspaceDown => {
                self.state.move_focus(1);
                self.sync_focus();
                Transition::Stay
            }
            Action::WorkspaceUp => {
                self.state.move_focus(-1);
                self.sync_focus();
                Transition::Stay
            }
            Action::WorkspaceRescan => {
                self.state.rescan();
                self.sync_focus();
                self.notice = None;
                Transition::Stay
            }
            Action::WorkspaceOpenSetup => Transition::Push(Box::new(SetupScreen::new())),
            Action::WorkspaceOpenFocused => match self.state.focused() {
                Some(a) => match &a.kind {
                    ArtifactKind::RawVcf => {
                        self.notice = None;
                        Transition::Push(Box::new(TransformScreen::new_ingest(Some(a))))
                    }
                    ArtifactKind::IngestedSet => {
                        self.notice = None;
                        Transition::Push(Box::new(TransformScreen::new_annotate(a)))
                    }
                    ArtifactKind::AnnotatedSet { .. } => {
                        match VariantScreen::new_for_annotated_set(a.path.clone()) {
                            Ok(screen) => {
                                self.notice = None;
                                Transition::Push(Box::new(screen))
                            }
                            Err(e) => {
                                self.notice = Some(format!("cannot open: {e}"));
                                Transition::Stay
                            }
                        }
                    }
                    ArtifactKind::ParquetFile => {
                        match VariantScreen::new_for_parquet(a.path.clone()) {
                            Ok(screen) => {
                                self.notice = None;
                                Transition::Push(Box::new(screen))
                            }
                            Err(e) => {
                                self.notice = Some(format!("cannot open: {e}"));
                                Transition::Stay
                            }
                        }
                    }
                    other => {
                        self.notice = Some(format!("no transform for {}", other.title()));
                        Transition::Stay
                    }
                },
                None => Transition::Stay,
            },
            _ => Transition::Stay,
        }
    }
}
