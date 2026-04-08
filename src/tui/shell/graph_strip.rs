#![allow(dead_code)]

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Paragraph, Widget};

use crate::tui::stages::{ArtifactKind, StageId, STAGES};
use crate::tui::theme;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NodeStatus {
    Done,
    Here,
    Ready,
    Queued,
}

impl NodeStatus {
    fn word(self) -> &'static str {
        match self {
            Self::Done => "done",
            Self::Here => "◉ here",
            Self::Ready => "ready",
            Self::Queued => "queued",
        }
    }

    fn bright(self) -> bool {
        matches!(self, Self::Done | Self::Here)
    }
}

#[derive(Debug, Clone)]
pub struct NodeView {
    pub id: StageId,
    pub label: &'static str,
    pub status: NodeStatus,
}

#[derive(Debug, Clone)]
pub struct GraphStripState {
    pub nodes: Vec<NodeView>,
    pub focus: Option<StageId>,
}

pub fn compute_graph(present: &[ArtifactKind], focus: Option<StageId>) -> GraphStripState {
    let mut nodes = Vec::with_capacity(STAGES.len());
    for stage in STAGES {
        let outs = stage.outputs();
        let ins = stage.inputs();
        let produced = !outs.is_empty() && outs.iter().all(|o| present.contains(o));
        let inputs_ok = ins.iter().all(|i| present.contains(i));
        let is_focus = focus.is_some_and(|f| f == stage.id());
        let status = if is_focus {
            NodeStatus::Here
        } else if produced {
            NodeStatus::Done
        } else if inputs_ok {
            NodeStatus::Ready
        } else {
            NodeStatus::Queued
        };
        nodes.push(NodeView {
            id: stage.id(),
            label: stage.label(),
            status,
        });
    }
    GraphStripState { nodes, focus }
}

const ARROW: &str = " → ";

fn node_cell_width(n: &NodeView) -> usize {
    n.label.chars().count().max(n.status.word().chars().count())
}

fn full_width(nodes: &[NodeView]) -> usize {
    if nodes.is_empty() {
        return 0;
    }
    let sum: usize = nodes.iter().map(node_cell_width).sum();
    sum + ARROW.chars().count() * (nodes.len() - 1)
}

fn node_style(status: NodeStatus) -> Style {
    let base = if status.bright() {
        Style::default().fg(theme::ACCENT).add_modifier(Modifier::BOLD)
    } else {
        Style::default().fg(theme::MUTED)
    };
    match status {
        NodeStatus::Done => base.fg(theme::OK).add_modifier(Modifier::BOLD),
        NodeStatus::Here => base,
        NodeStatus::Ready => Style::default().fg(theme::FG),
        NodeStatus::Queued => Style::default().fg(theme::MUTED),
    }
}

fn pad(s: &str, width: usize) -> String {
    let count = s.chars().count();
    if count >= width {
        s.to_string()
    } else {
        let mut out = String::from(s);
        for _ in 0..(width - count) {
            out.push(' ');
        }
        out
    }
}

fn windowed<'a>(
    nodes: &'a [NodeView],
    focus: Option<StageId>,
    max_width: usize,
) -> (Vec<&'a NodeView>, bool, bool) {
    let focus_idx = focus
        .and_then(|f| nodes.iter().position(|n| n.id == f))
        .unwrap_or(0);
    let marker_cost = 4;
    let mut lo = focus_idx;
    let mut hi = focus_idx;
    let mut width = node_cell_width(&nodes[focus_idx]);
    let budget = max_width.saturating_sub(marker_cost * 2);
    loop {
        let grew_left = if lo > 0 {
            let extra = node_cell_width(&nodes[lo - 1]) + ARROW.chars().count();
            if width + extra <= budget {
                lo -= 1;
                width += extra;
                true
            } else {
                false
            }
        } else {
            false
        };
        let grew_right = if hi + 1 < nodes.len() {
            let extra = node_cell_width(&nodes[hi + 1]) + ARROW.chars().count();
            if width + extra <= budget {
                hi += 1;
                width += extra;
                true
            } else {
                false
            }
        } else {
            false
        };
        if !grew_left && !grew_right {
            break;
        }
    }
    let left_overflow = lo > 0;
    let right_overflow = hi + 1 < nodes.len();
    (nodes[lo..=hi].iter().collect(), left_overflow, right_overflow)
}

pub fn render(state: &GraphStripState, label_area: Rect, status_area: Rect, buf: &mut Buffer) {
    if state.nodes.is_empty() {
        return;
    }
    let avail = label_area.width as usize;
    let (slice, left_over, right_over) = if full_width(&state.nodes) <= avail {
        (state.nodes.iter().collect::<Vec<_>>(), false, false)
    } else {
        windowed(&state.nodes, state.focus, avail)
    };

    let mut label_spans: Vec<Span> = Vec::new();
    let mut status_spans: Vec<Span> = Vec::new();
    let dim = Style::default().fg(theme::MUTED);

    if left_over {
        label_spans.push(Span::styled("< ", dim));
        status_spans.push(Span::styled("  ", dim));
    } else {
        label_spans.push(Span::raw("  "));
        status_spans.push(Span::raw("  "));
    }

    for (i, n) in slice.iter().enumerate() {
        if i > 0 {
            label_spans.push(Span::styled(ARROW, dim));
            status_spans.push(Span::raw("   "));
        }
        let w = node_cell_width(n);
        let style = node_style(n.status);
        label_spans.push(Span::styled(pad(n.label, w), style));
        status_spans.push(Span::styled(pad(n.status.word(), w), style));
    }

    if right_over {
        label_spans.push(Span::styled(" >", dim));
    }

    Paragraph::new(Line::from(label_spans)).render(label_area, buf);
    Paragraph::new(Line::from(status_spans)).render(status_area, buf);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compute_marks_done_ready_queued_and_here() {
        let present = [ArtifactKind::VcfDir, ArtifactKind::IngestedSet];
        let ingest_id = STAGES[0].id();
        let state = compute_graph(&present, Some(ingest_id));
        assert_eq!(state.nodes.len(), STAGES.len());
        let ingest = &state.nodes[0];
        assert_eq!(ingest.status, NodeStatus::Here);

        let state2 = compute_graph(&present, None);
        let ingest2 = &state2.nodes[0];
        assert_eq!(ingest2.status, NodeStatus::Done);
        let annotate = state2
            .nodes
            .iter()
            .find(|n| n.id == STAGES[1].id())
            .unwrap();
        assert_eq!(annotate.status, NodeStatus::Ready);
        let staar = state2
            .nodes
            .iter()
            .find(|n| n.id == STAGES[2].id())
            .unwrap();
        assert_eq!(staar.status, NodeStatus::Queued);
    }

    #[test]
    fn windowed_shrinks_to_focus_with_overflow_markers() {
        let state = compute_graph(&[], Some(STAGES[2].id()));
        let (slice, left, right) = windowed(&state.nodes, state.focus, 20);
        assert!(!slice.is_empty());
        assert!(slice.iter().any(|n| n.id == STAGES[2].id()));
        let _ = (left, right);
    }
}
