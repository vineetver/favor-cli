use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{Array, Int32Array, Int64Array, RecordBatch, StringArray};
use arrow::datatypes::SchemaRef;
use arrow::util::display::{ArrayFormatter, FormatOptions};
use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use parquet::arrow::ProjectionMask;
use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{
    Block, Borders, Clear, List, ListItem, ListState, Paragraph, StatefulWidget, Widget,
};
use ratatui::Frame;

use crate::column::Col;
use crate::error::CohortError;
use crate::store::cohort::sparse_g::SparseG;
use crate::store::cohort::variants::{CarrierList, VariantIndex};
use crate::tui::action::{Action, ActionScope};
use crate::tui::shell::{ErrorMessage, ScreenChrome, Shell};
use crate::tui::state::app::{AppState, Outcome, View};
use crate::tui::state::arrow_predicate::CompiledFilter;
use crate::tui::state::parquet_scroller::{ParquetScroller, RowFilterFactory};
use crate::tui::theme;

pub struct VariantSpec {
    pub parquet: PathBuf,
    pub title_prefix: Option<String>,
    pub summary: Option<String>,
    pub companion: Option<PathBuf>,
}

impl VariantSpec {
    pub fn for_file(parquet: PathBuf) -> Self {
        Self { parquet, title_prefix: None, summary: None, companion: None }
    }

    pub fn for_annotated_set(set_root: PathBuf) -> Result<Self, CohortError> {
        let parquet = pick_first_chrom_parquet(&set_root)?;
        let prefix = set_root
            .file_name()
            .map(|n| n.to_string_lossy().into_owned());
        Ok(Self { parquet, title_prefix: prefix, summary: None, companion: None })
    }

    pub fn for_staar_results(
        results_dir: PathBuf,
        companion: Option<PathBuf>,
    ) -> Result<Self, CohortError> {
        let parquet = first_results_parquet(&results_dir)?;
        let prefix = results_dir
            .file_name()
            .map(|n| n.to_string_lossy().into_owned());
        let summary = format!(
            "STAAR results · enter on a row opens {} filtered to that gene",
            companion
                .as_ref()
                .and_then(|p| p.file_name())
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "companion set".into())
        );
        Ok(Self { parquet, title_prefix: prefix, summary: Some(summary), companion })
    }
}

pub struct VariantViewState {
    pub title: String,
    pub summary: Option<String>,
    pub companion: Option<PathBuf>,
    pub scroller: ParquetScroller,
    pub display_columns: Vec<usize>,
    pub active_filter: Option<Arc<CompiledFilter>>,
    pub modal: VariantModal,
    pub chrom_cache: HashMap<String, ChromCacheEntry>,
    pub sample_names: Option<Vec<String>>,
    pub parent: Option<Box<VariantViewState>>,
}

pub enum VariantModal {
    None,
    Filter(FilterModalState),
    Columns(ColumnsModalState),
    Carrier(CarrierModalState),
}

pub struct FilterModalState {
    pub buf: String,
    pub parse_error: Option<String>,
}

pub struct ColumnsModalState {
    pub cursor: usize,
    pub selected: Vec<bool>,
}

pub struct CarrierModalState {
    pub vid: String,
    pub carriers: CarrierList,
    pub scroll: usize,
}

pub struct ChromCacheEntry {
    pub sparse_g: SparseG,
    pub index: VariantIndex,
}

impl VariantViewState {
    pub fn active_scope(&self) -> ActionScope {
        ActionScope::Variant
    }
}

pub fn open_variant(spec: VariantSpec) -> Result<VariantViewState, CohortError> {
    let scroller = ParquetScroller::open(&spec.parquet)?;
    let display_columns: Vec<usize> = (0..scroller.schema().fields().len()).collect();
    let leaf = spec
        .parquet
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_else(|| spec.parquet.display().to_string());
    let title = match spec.title_prefix {
        Some(p) => format!("{p} :: {leaf}"),
        None => leaf,
    };
    Ok(VariantViewState {
        title,
        summary: spec.summary,
        companion: spec.companion,
        scroller,
        display_columns,
        active_filter: None,
        modal: VariantModal::None,
        chrom_cache: HashMap::new(),
        sample_names: None,
        parent: None,
    })
}

fn with_initial_filter(
    mut view: VariantViewState,
    text: &str,
) -> Result<VariantViewState, CohortError> {
    let parsed = CompiledFilter::parse(text).map_err(CohortError::Input)?;
    let arc = Arc::new(parsed);
    let factory: Arc<dyn RowFilterFactory> = arc.clone();
    view.scroller.set_filter(Some(factory))?;
    view.active_filter = Some(arc);
    rebuild_projection(&mut view)?;
    Ok(view)
}

fn open_companion_for_focused_gene(
    view: &VariantViewState,
) -> Result<VariantViewState, CohortError> {
    let (batch, row) = view
        .scroller
        .focused_record()
        .ok_or_else(|| CohortError::Input("no focused row".into()))?;
    let schema = batch.schema();
    let idx = schema
        .index_of(Col::GeneName.as_str())
        .map_err(|_| CohortError::Input("row has no gene_name column".into()))?;
    let arr = batch
        .column(idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| CohortError::Input("gene_name column has unexpected type".into()))?;
    if arr.is_null(row) {
        return Err(CohortError::Input("gene_name is null".into()));
    }
    let gene = arr.value(row).to_string();
    let target = view
        .companion
        .clone()
        .ok_or_else(|| CohortError::Input("no companion artifact attached".into()))?;
    let prefix = target
        .file_name()
        .map(|n| n.to_string_lossy().into_owned());
    let spec = VariantSpec {
        parquet: target,
        title_prefix: prefix,
        summary: None,
        companion: None,
    };
    let text = format!("{} = \"{}\"", Col::GeneName.as_str(), gene);
    let next = open_variant(spec)?;
    with_initial_filter(next, &text)
}

fn rebuild_projection(view: &mut VariantViewState) -> Result<(), CohortError> {
    let n_fields = view.scroller.schema().fields().len();
    let mut wanted: Vec<bool> = vec![false; n_fields];
    for &i in &view.display_columns {
        if i < n_fields {
            wanted[i] = true;
        }
    }
    if let Some(filter) = &view.active_filter {
        let schema = view.scroller.schema().clone();
        for clause in filter.clauses() {
            if let Ok(idx) = schema.index_of(&clause.column) {
                wanted[idx] = true;
            }
        }
    }
    let leaves: Vec<usize> = wanted
        .iter()
        .enumerate()
        .filter_map(|(i, &k)| if k { Some(i) } else { None })
        .collect();
    let mask = ProjectionMask::leaves(view.scroller.metadata().parquet_schema(), leaves);
    view.scroller.set_projection(mask)
}

fn apply_filter_text(view: &mut VariantViewState, error: &mut Option<ErrorMessage>) {
    let VariantModal::Filter(modal) = &mut view.modal else {
        return;
    };
    let text = modal.buf.trim().to_string();
    if text.is_empty() {
        view.active_filter = None;
        if let Err(e) = view
            .scroller
            .set_filter(None)
            .and_then(|_| rebuild_projection(view))
        {
            *error = Some(ErrorMessage { text: format!("{e}") });
        }
        view.modal = VariantModal::None;
        return;
    }
    let parsed = match CompiledFilter::parse(&text) {
        Ok(c) => c,
        Err(e) => {
            modal.parse_error = Some(e);
            return;
        }
    };
    let arc = Arc::new(parsed);
    let factory: Arc<dyn RowFilterFactory> = arc.clone();
    if let Err(e) = view.scroller.set_filter(Some(factory)) {
        if let VariantModal::Filter(m) = &mut view.modal {
            m.parse_error = Some(format!("{e}"));
        }
        let _ = view.scroller.set_filter(None);
        return;
    }
    view.active_filter = Some(arc);
    if let Err(e) = rebuild_projection(view) {
        if let VariantModal::Filter(m) = &mut view.modal {
            m.parse_error = Some(format!("{e}"));
        }
        return;
    }
    view.modal = VariantModal::None;
}

fn clear_filter(view: &mut VariantViewState) {
    view.active_filter = None;
    let _ = view.scroller.set_filter(None);
    let _ = rebuild_projection(view);
}

fn open_filter(view: &mut VariantViewState) {
    let buf = view
        .scroller
        .filter_text()
        .map(|s| s.to_string())
        .unwrap_or_default();
    view.modal = VariantModal::Filter(FilterModalState { buf, parse_error: None });
}

fn open_column_picker(view: &mut VariantViewState) {
    let n = view.scroller.schema().fields().len();
    let mut selected = vec![false; n];
    for &i in &view.display_columns {
        if i < n {
            selected[i] = true;
        }
    }
    view.modal = VariantModal::Columns(ColumnsModalState { cursor: 0, selected });
}

fn close_column_picker(view: &mut VariantViewState, error: &mut Option<ErrorMessage>) {
    if let VariantModal::Columns(picker) =
        std::mem::replace(&mut view.modal, VariantModal::None)
    {
        view.display_columns = picker
            .selected
            .iter()
            .enumerate()
            .filter_map(|(i, &k)| if k { Some(i) } else { None })
            .collect();
        if view.display_columns.is_empty() {
            let n = view.scroller.schema().fields().len();
            view.display_columns = (0..n).collect();
        }
        if let Err(e) = rebuild_projection(view) {
            *error = Some(ErrorMessage { text: format!("{e}") });
        }
    }
}

fn open_carrier_view(view: &mut VariantViewState, error: &mut Option<ErrorMessage>) {
    let Some((batch, row)) = view.scroller.focused_record() else {
        return;
    };
    let schema = batch.schema();
    let vid = match find_vid(batch, &schema, row) {
        Some(v) => v,
        None => {
            *error = Some(ErrorMessage {
                text: "row has no vid and no (chromosome, position, ref_allele, alt_allele) tuple"
                    .into(),
            });
            return;
        }
    };
    let chrom = chrom_from_path_or_batch(view.scroller.path(), batch, &schema, row)
        .unwrap_or_default();
    let path = view.scroller.path().to_path_buf();
    match load_carriers(view, &path, &chrom, &vid) {
        Ok(carriers) => {
            view.modal = VariantModal::Carrier(CarrierModalState { vid, carriers, scroll: 0 });
        }
        Err(e) => {
            *error = Some(ErrorMessage { text: format!("{e}") });
        }
    }
}

fn close_modal(view: &mut VariantViewState, error: &mut Option<ErrorMessage>) {
    if matches!(view.modal, VariantModal::Columns(_)) {
        close_column_picker(view, error);
        return;
    }
    view.modal = VariantModal::None;
}

fn load_carriers(
    view: &mut VariantViewState,
    parquet_path: &Path,
    chrom: &str,
    vid: &str,
) -> Result<CarrierList, CohortError> {
    if chrom.is_empty() {
        return Err(CohortError::Input("row has no chromosome".into()));
    }
    let store_dir = locate_store(parquet_path).ok_or_else(|| {
        CohortError::DataMissing(format!(
            "no genotype store found near {}",
            parquet_path.display()
        ))
    })?;
    if view.sample_names.is_none() {
        view.sample_names = Some(load_sample_names(&store_dir)?);
    }
    let chrom_dir = store_dir.join(format!("chromosome={chrom}"));
    if !chrom_dir.exists() {
        return Err(CohortError::DataMissing(format!(
            "store has no chromosome {chrom}"
        )));
    }
    if !view.chrom_cache.contains_key(chrom) {
        let sparse_g = SparseG::open(&chrom_dir)?;
        let index = VariantIndex::load(&chrom_dir)?;
        view.chrom_cache
            .insert(chrom.to_string(), ChromCacheEntry { sparse_g, index });
    }
    let entry = view.chrom_cache.get(chrom).unwrap();
    let vvcf = entry
        .index
        .resolve_vid(vid)
        .ok_or_else(|| CohortError::DataMissing(format!("vid not in store: {vid}")))?;
    Ok(entry.sparse_g.load_variant(vvcf))
}

pub fn handle_modal_key(state: &mut AppState, key: KeyEvent) -> Outcome {
    let View::Variant(view) = &mut state.view else {
        return Outcome::Stay;
    };
    if matches!(key.code, KeyCode::Esc) {
        close_modal(view, &mut state.error);
        return Outcome::Stay;
    }
    match &mut view.modal {
        VariantModal::None => Outcome::Stay,
        VariantModal::Filter(modal) => {
            match (key.code, key.modifiers) {
                (KeyCode::Enter, _) => {
                    apply_filter_text(view, &mut state.error);
                }
                (KeyCode::Backspace, _) => {
                    modal.buf.pop();
                    modal.parse_error = None;
                }
                (KeyCode::Char('u'), m) if m.contains(KeyModifiers::CONTROL) => {
                    modal.buf.clear();
                    modal.parse_error = None;
                }
                (KeyCode::Char(c), m)
                    if !m.contains(KeyModifiers::CONTROL) && !m.contains(KeyModifiers::ALT) =>
                {
                    modal.buf.push(c);
                    modal.parse_error = None;
                }
                _ => {}
            }
            Outcome::Stay
        }
        VariantModal::Columns(picker) => {
            let n = picker.selected.len();
            match key.code {
                KeyCode::Enter => {
                    close_column_picker(view, &mut state.error);
                }
                KeyCode::Char('j') | KeyCode::Down => {
                    if n > 0 && picker.cursor + 1 < n {
                        picker.cursor += 1;
                    }
                }
                KeyCode::Char('k') | KeyCode::Up => {
                    if picker.cursor > 0 {
                        picker.cursor -= 1;
                    }
                }
                KeyCode::Char(' ') => {
                    if picker.cursor < n {
                        picker.selected[picker.cursor] = !picker.selected[picker.cursor];
                    }
                }
                _ => {}
            }
            Outcome::Stay
        }
        VariantModal::Carrier(panel) => {
            match key.code {
                KeyCode::Char('j') | KeyCode::Down => {
                    if panel.scroll + 1 < panel.carriers.len() {
                        panel.scroll += 1;
                    }
                }
                KeyCode::Char('k') | KeyCode::Up => {
                    panel.scroll = panel.scroll.saturating_sub(1);
                }
                _ => {}
            }
            Outcome::Stay
        }
    }
}

pub fn render(state: &mut AppState, frame: &mut Frame, area: Rect) {
    let View::Variant(view) = &state.view else {
        return;
    };
    let error_text = state.error.as_ref().map(|e| e.text.clone()).or_else(|| {
        match &view.modal {
            VariantModal::Filter(m) => m.parse_error.clone().map(|e| format!("error: {e}")),
            _ => None,
        }
    });
    let error = error_text.map(|text| ErrorMessage { text });
    let chrome = ScreenChrome {
        title: &view.title,
        status: None,
        error: error.as_ref(),
        scope: ActionScope::Variant,
        graph: None,
    };
    let modal = &view.modal;
    let body = |inner: Rect, buf: &mut Buffer| {
        let modal_height: u16 = match modal {
            VariantModal::Filter(_) => 2,
            VariantModal::Carrier(_) => 10,
            _ => 0,
        };
        let detail_height: u16 = if inner.height >= 20 { 8 } else { 0 };
        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(6),
                Constraint::Length(detail_height),
                Constraint::Length(modal_height),
            ])
            .split(inner);
        draw_table(view, v[0], buf);
        if detail_height > 0 && matches!(modal, VariantModal::None) {
            draw_detail(view, v[1], buf);
        }
        match modal {
            VariantModal::Filter(m) => draw_filter_bar(v[2], buf, m),
            VariantModal::Carrier(p) => draw_carrier_panel(view, v[2], buf, p),
            _ => {}
        }
        if let VariantModal::Columns(p) = modal {
            draw_column_picker(view, inner, buf, p);
        }
    };
    Shell::new(chrome, body).render(area, frame.buffer_mut());
}

pub fn handle_action(state: &mut AppState, action: Action) -> Outcome {
    state.error = None;
    if matches!(action, Action::VariantClose) {
        return pop_variant(state);
    }
    if matches!(action, Action::VariantOpenLinked) {
        return open_linked(state);
    }
    let View::Variant(view) = &mut state.view else {
        return Outcome::Stay;
    };
    let nav: Result<(), CohortError> = match action {
        Action::VariantScrollRowDown => {
            view.scroller.scroll_down();
            Ok(())
        }
        Action::VariantScrollRowUp => {
            view.scroller.scroll_up();
            Ok(())
        }
        Action::VariantPrevRowGroup => view.scroller.prev_row_group(),
        Action::VariantNextRowGroup => view.scroller.next_row_group(),
        Action::VariantJumpStart => view.scroller.goto_row_group(0),
        Action::VariantJumpEnd => {
            let last = view.scroller.row_group_count().saturating_sub(1);
            view.scroller.goto_row_group(last)
        }
        Action::VariantOpenFilter => {
            open_filter(view);
            Ok(())
        }
        Action::VariantFilterClear => {
            clear_filter(view);
            Ok(())
        }
        Action::VariantOpenColumnPicker => {
            open_column_picker(view);
            Ok(())
        }
        Action::VariantOpenCarrierView => {
            open_carrier_view(view, &mut state.error);
            Ok(())
        }
        Action::VariantCloseCarrierView | Action::VariantColumnPickerClose => {
            close_modal(view, &mut state.error);
            Ok(())
        }
        _ => Ok(()),
    };
    if let Err(e) = nav {
        state.error = Some(ErrorMessage { text: format!("{e}") });
    }
    Outcome::Stay
}

fn open_linked(state: &mut AppState) -> Outcome {
    let View::Variant(view) = &state.view else {
        return Outcome::Stay;
    };
    let next = match open_companion_for_focused_gene(view) {
        Ok(v) => v,
        Err(e) => {
            state.error = Some(ErrorMessage { text: format!("{e}") });
            return Outcome::Stay;
        }
    };
    let prev = std::mem::replace(&mut state.view, View::Workspace);
    let mut next = next;
    if let View::Variant(prev_view) = prev {
        next.parent = Some(Box::new(prev_view));
    }
    state.view = View::Variant(next);
    Outcome::Stay
}

fn pop_variant(state: &mut AppState) -> Outcome {
    let View::Variant(view) = &mut state.view else {
        return Outcome::Stay;
    };
    if let Some(parent) = view.parent.take() {
        *view = *parent;
        return Outcome::Stay;
    }
    state.view = View::Workspace;
    Outcome::Stay
}

fn draw_table(view: &VariantViewState, area: Rect, buf: &mut Buffer) {
    let block = Block::default()
        .borders(Borders::ALL)
        .title(format!(" {} ", header_text(view)))
        .border_style(Style::default().fg(theme::MUTED));
    let inner = block.inner(area);
    block.render(area, buf);

    let Some((batch, focus)) = view.scroller.focused_record() else {
        let msg = if view.scroller.row_group_count() == 0 {
            "  empty parquet"
        } else {
            "  no rows in this row group (filter excluded all)"
        };
        Paragraph::new(Line::from(Span::styled(
            msg,
            Style::default().fg(theme::MUTED),
        )))
        .render(inner, buf);
        return;
    };

    let schema = batch.schema();
    let cols: Vec<usize> = view
        .display_columns
        .iter()
        .filter_map(|i| schema.index_of(view.scroller.schema().field(*i).name()).ok())
        .collect();
    let gutter_w = theme::FOCUS_GLYPH.chars().count() + 1;
    let body_w = (inner.width as usize).saturating_sub(gutter_w);
    let widths = compute_widths(batch, &cols, body_w);
    let blank_gutter: String = " ".repeat(gutter_w);
    let mut header_spans = vec![Span::raw(blank_gutter.clone())];
    header_spans.extend(cols.iter().zip(widths.iter()).map(|(c, w)| {
        Span::styled(
            pad(schema.field(*c).name(), *w),
            Style::default().fg(theme::ACCENT).bold(),
        )
    }));
    let header_line = Line::from(header_spans);

    let visible_rows = inner.height.saturating_sub(1) as usize;
    let total = batch.num_rows();
    let start = focus
        .saturating_sub(visible_rows / 2)
        .min(total.saturating_sub(visible_rows.max(1)));
    let end = (start + visible_rows).min(total);

    let mut lines: Vec<Line> = Vec::with_capacity(end - start + 1);
    lines.push(header_line);
    for row in start..end {
        let is_focus = row == focus;
        let mut spans: Vec<Span> = Vec::with_capacity(cols.len() + 1);
        if is_focus {
            spans.push(Span::styled(
                format!("{} ", theme::FOCUS_GLYPH),
                Style::default().fg(theme::ACCENT).bold(),
            ));
        } else {
            spans.push(Span::raw(blank_gutter.clone()));
        }
        let style = if is_focus {
            Style::default().fg(theme::FG).bold()
        } else {
            Style::default().fg(theme::FG)
        };
        spans.extend(
            cols.iter()
                .zip(widths.iter())
                .map(|(c, w)| Span::styled(pad(&format_cell(batch, *c, row), *w), style)),
        );
        lines.push(Line::from(spans));
    }
    Paragraph::new(lines).render(inner, buf);
}

fn header_text(view: &VariantViewState) -> String {
    let total = view.scroller.total_rows_estimate();
    let rg = view.scroller.row_group_count();
    let cur = view.scroller.current_row_group() + 1;
    let visible = view.scroller.current_batch_len();
    let filter = view.scroller.filter_text().unwrap_or("none");
    let base = format!(
        "{} · row group {}/{} · {} visible · {} total est · filter: {}",
        view.title, cur, rg.max(1), visible, total, filter
    );
    match view.summary.as_deref() {
        Some(s) => format!("{base} · {s}"),
        None => base,
    }
}

fn draw_detail(view: &VariantViewState, area: Rect, buf: &mut Buffer) {
    let inner = area;
    let Some((batch, row)) = view.scroller.focused_record() else {
        return;
    };
    let schema = batch.schema();
    let pairs: Vec<(String, String)> = (0..schema.fields().len())
        .map(|i| {
            let name = schema.field(i).name().to_string();
            let val = format_cell(batch, i, row);
            (name, val)
        })
        .collect();

    let half = pairs.len().div_ceil(2);
    let max_rows = inner.height as usize;
    let mut lines: Vec<Line> = Vec::with_capacity(max_rows);
    let shown = half.min(max_rows);
    for i in 0..shown {
        let mut spans = Vec::new();
        spans.push(Span::styled(
            format!("  {:<22} ", pairs[i].0),
            Style::default().fg(theme::MUTED),
        ));
        spans.push(Span::styled(
            pad(&pairs[i].1, 30),
            Style::default().fg(theme::FG),
        ));
        if let Some(rhs) = pairs.get(i + half) {
            spans.push(Span::styled(
                format!("  {:<22} ", rhs.0),
                Style::default().fg(theme::MUTED),
            ));
            spans.push(Span::styled(
                pad(&rhs.1, 30),
                Style::default().fg(theme::FG),
            ));
        }
        lines.push(Line::from(spans));
    }
    if shown < half && shown > 0 {
        let hidden = half - shown;
        lines.pop();
        lines.push(Line::from(Span::styled(
            format!("  ↓ {hidden} more rows"),
            Style::default().fg(theme::MUTED),
        )));
    }
    Paragraph::new(lines).render(inner, buf);
}

fn draw_filter_bar(area: Rect, buf: &mut Buffer, modal: &FilterModalState) {
    let input = Line::from(vec![
        Span::styled(" filter ", Style::default().fg(theme::WARN).bold()),
        Span::styled("/ ", Style::default().fg(theme::MUTED)),
        Span::styled(modal.buf.as_str(), Style::default().fg(theme::FG)),
        Span::styled("_", Style::default().fg(theme::ACCENT)),
    ]);
    let hint = Line::from(Span::styled(
        " e.g. af < 0.01 AND consequence contains missense  (enter apply, esc cancel)",
        theme::hint_bar_style(),
    ));
    Paragraph::new(vec![input, hint]).render(area, buf);
}

fn draw_carrier_panel(
    view: &VariantViewState,
    area: Rect,
    buf: &mut Buffer,
    panel: &CarrierModalState,
) {
    let list = &panel.carriers;
    let mut lines: Vec<Line> = vec![Line::from(Span::styled(
        format!("  {} carriers of {}", list.len(), panel.vid),
        Style::default().fg(theme::ACCENT).bold(),
    ))];
    let names = view.sample_names.as_deref().unwrap_or(&[]);
    let visible = (area.height as usize).saturating_sub(1);
    if panel.scroll > 0 {
        lines.push(Line::from(Span::styled(
            format!("    ↑ {} above", panel.scroll),
            Style::default().fg(theme::MUTED),
        )));
    }
    let body_cap = visible.saturating_sub(usize::from(panel.scroll > 0));
    let total = list.entries.len();
    let below = total.saturating_sub((panel.scroll + body_cap).min(total));
    let body_take = body_cap.saturating_sub(usize::from(below > 0));
    for entry in list.entries.iter().skip(panel.scroll).take(body_take) {
        let name = names
            .get(entry.sample_idx as usize)
            .map(|s| s.as_str())
            .unwrap_or("?");
        lines.push(Line::from(vec![
            Span::styled(format!("    {name:<32}"), Style::default().fg(theme::FG)),
            Span::styled(
                format!(" dosage={}", entry.dosage),
                Style::default().fg(theme::WARN),
            ),
        ]));
    }
    if below > 0 {
        lines.push(Line::from(Span::styled(
            format!("    ↓ {below} more"),
            Style::default().fg(theme::MUTED),
        )));
    }
    Paragraph::new(lines).render(area, buf);
}

fn draw_column_picker(
    view: &VariantViewState,
    area: Rect,
    buf: &mut Buffer,
    picker: &ColumnsModalState,
) {
    let overlay = centered(area, 50, 70);
    Clear.render(overlay, buf);
    let inner = overlay;
    let schema = view.scroller.schema();
    let items: Vec<ListItem> = schema
        .fields()
        .iter()
        .enumerate()
        .map(|(i, f)| {
            let mark = if picker.selected[i] { "[x]" } else { "[ ]" };
            let style = if i == picker.cursor {
                Style::default().fg(theme::ACCENT).bold()
            } else {
                Style::default().fg(theme::FG)
            };
            ListItem::new(Line::from(vec![
                Span::styled(format!(" {mark} "), Style::default().fg(theme::WARN)),
                Span::styled(f.name().to_string(), style),
            ]))
        })
        .collect();
    let mut state = ListState::default();
    state.select(Some(picker.cursor));
    StatefulWidget::render(List::new(items), inner, buf, &mut state);
}

fn first_results_parquet(results_dir: &Path) -> Result<PathBuf, CohortError> {
    if let Some(p) = first_parquet_in(results_dir) {
        return Ok(p);
    }
    Err(CohortError::Input(format!(
        "no parquet files under {}",
        results_dir.display()
    )))
}

fn pick_first_chrom_parquet(set_root: &Path) -> Result<PathBuf, CohortError> {
    let mut chrom_dirs: Vec<(String, PathBuf)> = Vec::new();
    let entries = std::fs::read_dir(set_root)
        .map_err(|e| CohortError::Resource(format!("read {}: {e}", set_root.display())))?;
    for entry in entries {
        let Ok(entry) = entry else { continue };
        let name = entry.file_name().to_string_lossy().into_owned();
        let Some(chrom) = name.strip_prefix("chromosome=") else {
            continue;
        };
        let dir = entry.path();
        if !dir.is_dir() {
            continue;
        }
        chrom_dirs.push((chrom.to_string(), dir));
    }
    if chrom_dirs.is_empty() {
        return Err(CohortError::Input(format!(
            "no chromosome=* directories under {}",
            set_root.display()
        )));
    }
    chrom_dirs.sort_by_key(|(c, _)| chrom_sort_key(c));
    for (_, dir) in chrom_dirs {
        if let Some(p) = first_parquet_in(&dir) {
            return Ok(p);
        }
    }
    Err(CohortError::Input(format!(
        "no parquet files under {}",
        set_root.display()
    )))
}

fn first_parquet_in(dir: &Path) -> Option<PathBuf> {
    let mut files: Vec<PathBuf> = std::fs::read_dir(dir)
        .ok()?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .filter(|p| p.extension().map(|x| x == "parquet").unwrap_or(false))
        .collect();
    files.sort();
    files.into_iter().next()
}

fn chrom_sort_key(chrom: &str) -> (u8, u32) {
    if let Ok(n) = chrom.trim_start_matches("chr").parse::<u32>() {
        return (0, n);
    }
    let bare = chrom.trim_start_matches("chr");
    match bare {
        "X" => (1, 0),
        "Y" => (1, 1),
        "M" | "MT" => (1, 2),
        _ => (2, 0),
    }
}

fn locate_store(parquet_path: &Path) -> Option<PathBuf> {
    let mut p = parquet_path.parent()?;
    loop {
        let candidate = p.join(".genotype_store");
        if candidate.join("manifest.json").exists() {
            return Some(candidate);
        }
        match p.parent() {
            Some(parent) => p = parent,
            None => return None,
        }
    }
}

fn load_sample_names(store_dir: &Path) -> Result<Vec<String>, CohortError> {
    let path = store_dir.join("samples.txt");
    let content = std::fs::read_to_string(&path)
        .map_err(|e| CohortError::Resource(format!("read {}: {e}", path.display())))?;
    Ok(content.lines().map(|l| l.trim().to_string()).collect())
}

fn find_vid(batch: &RecordBatch, schema: &SchemaRef, row: usize) -> Option<String> {
    if let Ok(idx) = schema.index_of(Col::Vid.as_str()) {
        let arr = batch.column(idx);
        if let Some(s) = arr.as_any().downcast_ref::<StringArray>() {
            if !s.is_null(row) {
                return Some(s.value(row).to_string());
            }
        }
    }
    let chrom_idx = schema.index_of(Col::Chromosome.as_str()).ok();
    let pos_idx = schema.index_of(Col::Position.as_str()).ok();
    let ref_idx = schema.index_of(Col::RefAllele.as_str()).ok();
    let alt_idx = schema.index_of(Col::AltAllele.as_str()).ok();
    let (chrom, pos, refa, alta) = match (chrom_idx, pos_idx, ref_idx, alt_idx) {
        (Some(c), Some(p), Some(r), Some(a)) => (c, p, r, a),
        _ => return None,
    };
    let chrom = string_at(batch.column(chrom), row)?;
    let pos = int_at(batch.column(pos), row)? as u32;
    let refa = string_at(batch.column(refa), row)?;
    let alta = string_at(batch.column(alta), row)?;
    Some(crate::types::format_vid(&chrom, pos, &refa, &alta))
}

fn chrom_from_path_or_batch(
    path: &Path,
    batch: &RecordBatch,
    schema: &SchemaRef,
    row: usize,
) -> Option<String> {
    for component in path.parent()?.components() {
        let s = component.as_os_str().to_string_lossy();
        if let Some(rest) = s.strip_prefix("chromosome=") {
            return Some(rest.to_string());
        }
    }
    let idx = schema.index_of(Col::Chromosome.as_str()).ok()?;
    string_at(batch.column(idx), row)
}

fn string_at(array: &Arc<dyn Array>, row: usize) -> Option<String> {
    if let Some(s) = array.as_any().downcast_ref::<StringArray>() {
        if !s.is_null(row) {
            return Some(s.value(row).to_string());
        }
    }
    None
}

fn int_at(array: &Arc<dyn Array>, row: usize) -> Option<i64> {
    if let Some(a) = array.as_any().downcast_ref::<Int32Array>() {
        if !a.is_null(row) {
            return Some(a.value(row) as i64);
        }
    }
    if let Some(a) = array.as_any().downcast_ref::<Int64Array>() {
        if !a.is_null(row) {
            return Some(a.value(row));
        }
    }
    None
}

fn compute_widths(batch: &RecordBatch, cols: &[usize], total_width: usize) -> Vec<usize> {
    if cols.is_empty() {
        return Vec::new();
    }
    let per = (total_width / cols.len()).clamp(8, 28);
    let mut widths = vec![per; cols.len()];
    for (idx, &c) in cols.iter().enumerate() {
        let header_len = batch.schema().field(c).name().len();
        if header_len > widths[idx] {
            widths[idx] = header_len.min(28);
        }
    }
    widths
}

fn pad(s: &str, w: usize) -> String {
    let mut out = String::with_capacity(w + 1);
    let mut count = 0;
    for ch in s.chars() {
        if count + 1 >= w {
            out.push('…');
            count += 1;
            break;
        }
        out.push(ch);
        count += 1;
    }
    while count < w {
        out.push(' ');
        count += 1;
    }
    out.push(' ');
    out
}

fn format_cell(batch: &RecordBatch, col: usize, row: usize) -> String {
    let array = batch.column(col);
    if array.is_null(row) {
        return "·".to_string();
    }
    let opts = FormatOptions::default();
    match ArrayFormatter::try_new(array.as_ref(), &opts) {
        Ok(f) => f.value(row).to_string(),
        Err(_) => String::new(),
    }
}

fn centered(area: Rect, pct_w: u16, pct_h: u16) -> Rect {
    let w = area.width * pct_w / 100;
    let h = area.height * pct_h / 100;
    let x = area.x + (area.width.saturating_sub(w)) / 2;
    let y = area.y + (area.height.saturating_sub(h)) / 2;
    Rect { x, y, width: w, height: h }
}
