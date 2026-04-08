use std::path::PathBuf;

use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, ListState, Paragraph, Widget};
use ratatui::Frame;

use crate::config::{Config, Tier};
use crate::tui::action::{Action, ActionScope};
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::screens::stage_view::StageView;
use crate::tui::screens::variant::{VariantScreen, VariantSpec};
use crate::tui::shell::graph_strip::compute_graph;
use crate::tui::shell::{ErrorMessage, ScreenChrome, Shell};
use crate::tui::stages::types::{ArtifactKind as StageArtifact, SessionCtx};
use crate::tui::stages::{self, Stage, SETUP_STAGE};
use crate::tui::state::artifacts::{Artifact, ArtifactKind};
use crate::tui::state::workspace::WorkspaceState;
use crate::tui::state::SessionState;
use crate::tui::theme;

pub struct WorkspaceScreen {
    title: String,
    state: WorkspaceState,
    list_state: ListState,
    error: Option<ErrorMessage>,
    pending_focus: Option<PathBuf>,
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
            error: None,
            pending_focus: None,
        }
    }

    fn try_apply_pending_focus(&mut self) {
        let Some(target) = self.pending_focus.as_ref() else {
            return;
        };
        if let Some(idx) = self.state.artifacts.iter().position(|a| &a.path == target) {
            self.state.focus = idx;
            self.pending_focus = None;
        }
    }

    fn browse_focused(&mut self) -> Transition {
        let Some(a) = self.state.focused() else {
            return Transition::Stay;
        };
        let path = a.path.clone();
        let spec = match &a.kind {
            ArtifactKind::AnnotatedSet { .. } => VariantSpec::for_annotated_set(path),
            ArtifactKind::ParquetFile => Ok(VariantSpec::for_file(path)),
            ArtifactKind::StaarResults => {
                let companion = self
                    .state
                    .artifacts
                    .iter()
                    .find(|x| matches!(x.kind, ArtifactKind::AnnotatedSet { .. }))
                    .map(|x| x.path.clone());
                VariantSpec::for_staar_results(path, companion)
            }
            _ => {
                self.error = Some(ErrorMessage {
                    text: format!("cannot browse {}", a.kind.title()),
                });
                return Transition::Stay;
            }
        };
        match spec.and_then(VariantScreen::open) {
            Ok(screen) => {
                self.error = None;
                Transition::Push(Box::new(screen))
            }
            Err(e) => {
                self.error = Some(ErrorMessage {
                    text: format!("cannot open: {e}"),
                });
                Transition::Stay
            }
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

fn classify(kind: &ArtifactKind) -> Option<StageArtifact> {
    match kind {
        ArtifactKind::RawVcf => Some(StageArtifact::VcfDir),
        ArtifactKind::IngestedSet => Some(StageArtifact::IngestedSet),
        ArtifactKind::AnnotatedSet { .. } => Some(StageArtifact::AnnotatedSet),
        ArtifactKind::GenotypeStore => Some(StageArtifact::GenotypeStore),
        ArtifactKind::StaarResults => Some(StageArtifact::StaarResults),
        _ => None,
    }
}

fn present_kinds(state: &WorkspaceState) -> Vec<StageArtifact> {
    let mut kinds: Vec<StageArtifact> = Vec::new();
    for a in &state.artifacts {
        if let Some(k) = classify(&a.kind) {
            if !kinds.contains(&k) {
                kinds.push(k);
            }
        }
    }
    kinds
}

fn next_stage_for(kind: &ArtifactKind) -> Option<&'static dyn Stage> {
    stages::next_for(classify(kind)?)
}

fn can_browse(kind: &ArtifactKind) -> bool {
    matches!(
        kind,
        ArtifactKind::AnnotatedSet { .. }
            | ArtifactKind::ParquetFile
            | ArtifactKind::StaarResults
    )
}

fn open_next(stage: &'static dyn Stage, art: &Artifact) -> Box<dyn Screen> {
    let cfg = Config::load().ok();
    let tier = cfg.as_ref().map(|c| c.data.tier).unwrap_or(Tier::Base);
    let data_root = cfg.map(|c| c.root_dir()).unwrap_or_default();
    let ctx = SessionCtx {
        data_root: &data_root,
        tier,
        focused: Some(art.path.as_path()),
    };
    Box::new(StageView::new(stage, ctx))
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

    fn draw(&mut self, frame: &mut Frame, area: Rect) {
        let present = present_kinds(&self.state);
        let focus_stage = self
            .state
            .focused()
            .and_then(|a| next_stage_for(&a.kind))
            .map(|s| s.id());
        let graph = compute_graph(&present, focus_stage);

        let list_title = if self.state.scanning {
            format!(" Artifacts ({}, scanning…) ", self.state.artifacts.len())
        } else {
            format!(" Artifacts ({}) ", self.state.artifacts.len())
        };

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
                    "  next",
                    Style::default().fg(theme::ACCENT),
                )));
                if let Some(stage) = next_stage_for(&a.kind) {
                    lines.push(Line::from(Span::styled(
                        format!("    {}", stage.label()),
                        Style::default().fg(theme::FG),
                    )));
                }
                if can_browse(&a.kind) {
                    lines.push(Line::from(Span::styled(
                        "    browse",
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

        let chrome = ScreenChrome {
            title: &self.title,
            status: None,
            error: self.error.as_ref(),
            scope: ActionScope::Workspace,
            graph: Some(graph),
        };

        let list_state = &mut self.list_state;
        let body = |inner: Rect, buf: &mut Buffer| {
            let cols = Layout::default()
                .direction(Direction::Horizontal)
                .constraints([Constraint::Percentage(55), Constraint::Percentage(45)])
                .split(inner);
            let list = List::new(items)
                .block(
                    Block::default()
                        .borders(Borders::ALL)
                        .title(list_title)
                        .border_style(Style::default().fg(theme::ACCENT)),
                )
                .highlight_style(Style::default().bg(theme::MUTED).fg(theme::FG))
                .highlight_symbol(" > ");
            ratatui::widgets::StatefulWidget::render(list, cols[0], buf, list_state);
            let detail = Paragraph::new(detail_lines).block(
                Block::default()
                    .borders(Borders::ALL)
                    .title(" Detail ")
                    .border_style(Style::default().fg(theme::MUTED)),
            );
            detail.render(cols[1], buf);
        };
        Shell::new(chrome, body).render(area, frame.buffer_mut());
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Workspace
    }

    fn on_other_event(&mut self, event: &AppEvent) -> Transition {
        if matches!(event, AppEvent::Tick) && self.state.drain_scan() {
            self.try_apply_pending_focus();
            self.sync_focus();
        }
        Transition::Stay
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
                self.error = None;
                Transition::Stay
            }
            Action::WorkspaceOpenSetup => {
                let cfg = Config::load().ok();
                let tier = cfg.as_ref().map(|c| c.data.tier).unwrap_or(Tier::Base);
                let data_root = cfg.map(|c| c.root_dir()).unwrap_or_default();
                let ctx = SessionCtx {
                    data_root: &data_root,
                    tier,
                    focused: None,
                };
                Transition::Push(Box::new(StageView::new(&SETUP_STAGE, ctx)))
            }
            Action::WorkspaceOpenFocused => {
                let Some(a) = self.state.focused() else {
                    return Transition::Stay;
                };
                if let Some(stage) = next_stage_for(&a.kind) {
                    self.error = None;
                    return Transition::Push(open_next(stage, a));
                }
                if can_browse(&a.kind) {
                    return self.browse_focused();
                }
                self.error = Some(ErrorMessage {
                    text: format!("no next stage for {}", a.kind.title()),
                });
                Transition::Stay
            }
            Action::WorkspaceBrowseFocused => self.browse_focused(),
            _ => Transition::Stay,
        }
    }

    fn contribute_session(&self, state: &mut SessionState) {
        state.cwd = self.state.cwd.clone();
        if let Some(a) = self.state.focused() {
            state.last_artifact = Some(a.path.clone());
        }
    }

    fn set_session_error(&mut self, msg: String) {
        self.error = Some(ErrorMessage { text: msg });
    }

    fn restore_session(&mut self, state: &SessionState) {
        if let Some(path) = &state.last_artifact {
            self.pending_focus = Some(path.clone());
            self.try_apply_pending_focus();
            self.sync_focus();
        }
    }
}
