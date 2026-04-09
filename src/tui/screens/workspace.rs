use std::path::PathBuf;

use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, Paragraph, Widget};
use ratatui::Frame;

use crate::config::{Config, Tier};
use crate::tui::action::{Action, ActionScope};
use crate::tui::screens::stage_view::StageViewState;
use crate::tui::screens::variant::{open_variant, VariantSpec};
use crate::tui::shell::graph_strip::compute_graph;
use crate::tui::shell::{ErrorMessage, ScreenChrome, Shell};
use crate::tui::stages::types::{ArtifactKind as StageArtifact, SessionCtx};
use crate::tui::stages::{self, Stage, SETUP_STAGE};
use crate::tui::state::app::{AppState, Outcome, View};
use crate::tui::state::artifacts::{Artifact, ArtifactKind};
use crate::tui::state::workspace::WorkspaceState;
use crate::tui::theme;

const TITLE: &str = "Workspace";

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

fn build_stage_view(stage: &'static dyn Stage, focused: Option<&Artifact>) -> StageViewState {
    let cfg = Config::load().ok();
    let tier = cfg.as_ref().map(|c| c.data.tier).unwrap_or(Tier::Base);
    let ctx = SessionCtx {
        tier,
        focused: focused.map(|a| a.path.as_path()),
    };
    StageViewState::new(stage, ctx)
}

pub fn render(state: &mut AppState, frame: &mut Frame, area: Rect) {
    let present = present_kinds(&state.workspace);
    let focus_stage = state
        .workspace
        .focused()
        .and_then(|a| next_stage_for(&a.kind))
        .map(|s| s.id());
    let graph = compute_graph(&present, focus_stage);

    let list_title = if state.workspace.scanning {
        format!(" Artifacts ({}, scanning…) ", state.workspace.artifacts.len())
    } else {
        format!(" Artifacts ({}) ", state.workspace.artifacts.len())
    };

    let items: Vec<ListItem> = state
        .workspace
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

    let detail_lines: Vec<Line> = match state.workspace.focused() {
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
            let msg = if state.workspace.scanning {
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
        title: TITLE,
        status: None,
        error: state.error.as_ref(),
        scope: ActionScope::Workspace,
        graph: Some(graph),
    };

    let list_state = &mut state.list_state;
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

pub fn on_tick(state: &mut AppState) {
    if state.workspace.drain_scan() {
        state.try_apply_pending_focus();
        state.sync_focus();
    }
}

pub fn handle_action(state: &mut AppState, action: Action) -> Outcome {
    match action {
        Action::Quit => Outcome::Quit,
        Action::WorkspaceDown => {
            state.workspace.move_focus(1);
            state.sync_focus();
            Outcome::Stay
        }
        Action::WorkspaceUp => {
            state.workspace.move_focus(-1);
            state.sync_focus();
            Outcome::Stay
        }
        Action::WorkspaceRescan => {
            state.workspace.rescan();
            state.sync_focus();
            state.error = None;
            Outcome::Stay
        }
        Action::WorkspaceOpenSetup => {
            let stage_view = build_stage_view(&SETUP_STAGE, None);
            state.view = View::Stage(stage_view);
            Outcome::Stay
        }
        Action::WorkspaceOpenFocused => {
            let Some(a) = state.workspace.focused() else {
                return Outcome::Stay;
            };
            if let Some(stage) = next_stage_for(&a.kind) {
                state.error = None;
                let view = build_stage_view(stage, Some(a));
                state.view = View::Stage(view);
                return Outcome::Stay;
            }
            if can_browse(&a.kind) {
                return browse_focused(state);
            }
            state.error = Some(ErrorMessage {
                text: format!("no next stage for {}", a.kind.title()),
            });
            Outcome::Stay
        }
        Action::WorkspaceBrowseFocused => browse_focused(state),
        _ => Outcome::Stay,
    }
}

fn browse_focused(state: &mut AppState) -> Outcome {
    let Some(a) = state.workspace.focused() else {
        return Outcome::Stay;
    };
    let path = a.path.clone();
    let companion: Option<PathBuf> = state
        .workspace
        .artifacts
        .iter()
        .find(|x| matches!(x.kind, ArtifactKind::AnnotatedSet { .. }))
        .map(|x| x.path.clone());
    let spec = match &a.kind {
        ArtifactKind::AnnotatedSet { .. } => VariantSpec::for_annotated_set(path),
        ArtifactKind::ParquetFile => Ok(VariantSpec::for_file(path)),
        ArtifactKind::StaarResults => VariantSpec::for_staar_results(path, companion),
        _ => {
            state.error = Some(ErrorMessage {
                text: format!("cannot browse {}", a.kind.title()),
            });
            return Outcome::Stay;
        }
    };
    match spec.and_then(open_variant) {
        Ok(view) => {
            state.error = None;
            state.view = View::Variant(view);
            Outcome::Stay
        }
        Err(e) => {
            state.error = Some(ErrorMessage {
                text: format!("cannot open: {e}"),
            });
            Outcome::Stay
        }
    }
}
