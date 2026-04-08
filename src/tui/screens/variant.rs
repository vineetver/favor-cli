use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::{
    Array, BooleanArray, Float32Array, Float64Array, Int32Array, Int64Array, RecordBatch,
    StringArray, UInt32Array, UInt64Array,
};
use arrow::datatypes::{DataType, SchemaRef};
use crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
use parquet::arrow::ProjectionMask;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Clear, List, ListItem, ListState, Paragraph};
use ratatui::Frame;

use crate::column::Col;
use crate::error::CohortError;
use crate::staar::carrier::reader::{CarrierList, VariantIndex};
use crate::staar::sparse_g::SparseG;
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::artifact_view::{ArtifactFilter, ArtifactOutcome, ArtifactView, Predicate};
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::state::arrow_predicate::CompiledFilter;
use crate::tui::state::parquet_scroller::{ParquetScroller, RowFilterFactory};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

pub struct VariantScreen {
    title: String,
    summary: Option<String>,
    companion: Option<PathBuf>,
    scroller: ParquetScroller,
    display_columns: Vec<usize>,
    active_filter: Option<Arc<CompiledFilter>>,
    modal: VariantModal,
    chrom_cache: HashMap<String, ChromCacheEntry>,
    sample_names: Option<Vec<String>>,
    error: Option<String>,
    pending_open: Option<(PathBuf, Option<ArtifactFilter>)>,
}

enum VariantModal {
    None,
    Filter(FilterModal),
    Columns(ColumnsModal),
    Carrier(CarrierModal),
}

struct FilterModal {
    buf: String,
    parse_error: Option<String>,
}

struct ColumnsModal {
    cursor: usize,
    selected: Vec<bool>,
}

struct CarrierModal {
    vid: String,
    carriers: Option<CarrierList>,
    error: Option<String>,
    scroll: usize,
}

struct ChromCacheEntry {
    sparse_g: SparseG,
    index: VariantIndex,
}

impl VariantScreen {
    fn open(
        path: PathBuf,
        title_prefix: Option<&str>,
        summary: Option<String>,
        companion: Option<PathBuf>,
    ) -> Result<Self, CohortError> {
        let scroller = ParquetScroller::open(&path)?;
        let display_columns: Vec<usize> = (0..scroller.schema().fields().len()).collect();
        let leaf = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());
        let title = match title_prefix {
            Some(p) => format!("{p} :: {leaf}"),
            None => leaf,
        };
        Ok(Self {
            title,
            summary,
            companion,
            scroller,
            display_columns,
            active_filter: None,
            modal: VariantModal::None,
            chrom_cache: HashMap::new(),
            sample_names: None,
            error: None,
            pending_open: None,
        })
    }

    pub fn new_for_parquet(path: PathBuf) -> Result<Self, CohortError> {
        Self::open(path, None, None, None)
    }

    pub fn new_for_annotated_set(set_path: PathBuf) -> Result<Self, CohortError> {
        let parquet = pick_first_chrom_parquet(&set_path)?;
        let prefix = set_path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        Self::open(parquet, Some(&prefix), None, None)
    }

    pub fn new_for_staar_results(
        results_path: PathBuf,
        companion_annotated_set: Option<PathBuf>,
    ) -> Result<Self, CohortError> {
        let parquet = first_results_parquet(&results_path)?;
        let prefix = results_path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        let summary = format!(
            "STAAR results · enter on a row opens {} filtered to that gene",
            companion_annotated_set
                .as_ref()
                .and_then(|p| p.file_name())
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_else(|| "companion set".into())
        );
        Self::open(parquet, Some(&prefix), Some(summary), companion_annotated_set)
    }

    pub fn with_initial_filter(mut self, text: &str) -> Result<Self, CohortError> {
        let parsed = CompiledFilter::parse(text).map_err(CohortError::Input)?;
        let arc = Arc::new(parsed);
        let factory: Arc<dyn RowFilterFactory> = arc.clone();
        self.scroller.set_filter(Some(factory))?;
        self.active_filter = Some(arc);
        self.rebuild_projection()?;
        Ok(self)
    }

    fn try_cross_open(&mut self) {
        let Some((batch, row)) = self.scroller.focused_record() else {
            return;
        };
        let schema = batch.schema();
        let Ok(idx) = schema.index_of(Col::GeneName.as_str()) else {
            self.error = Some("row has no gene_name column".into());
            return;
        };
        let arr = batch.column(idx);
        let Some(s) = arr.as_any().downcast_ref::<StringArray>() else {
            return;
        };
        if s.is_null(row) {
            return;
        }
        let gene = s.value(row).to_string();
        let Some(target) = self.companion.clone() else {
            self.error = Some("no companion artifact attached to this view".into());
            return;
        };
        self.pending_open = Some((
            target,
            Some(ArtifactFilter::new(Col::GeneName, Predicate::Eq(gene))),
        ));
    }

    fn rebuild_projection(&mut self) -> Result<(), CohortError> {
        let n_fields = self.scroller.schema().fields().len();
        let mut wanted: Vec<bool> = vec![false; n_fields];
        for &i in &self.display_columns {
            if i < n_fields {
                wanted[i] = true;
            }
        }
        if let Some(filter) = &self.active_filter {
            let schema = self.scroller.schema().clone();
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
        let mask = ProjectionMask::leaves(
            self.scroller.metadata().parquet_schema(),
            leaves,
        );
        self.scroller.set_projection(mask)
    }

    fn apply_filter_text(&mut self) {
        let VariantModal::Filter(modal) = &mut self.modal else {
            return;
        };
        let text = modal.buf.trim().to_string();
        if text.is_empty() {
            self.active_filter = None;
            if let Err(e) = self
                .scroller
                .set_filter(None)
                .and_then(|_| self.rebuild_projection())
            {
                self.error = Some(format!("{e}"));
            }
            self.modal = VariantModal::None;
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
        if let Err(e) = self.scroller.set_filter(Some(factory)) {
            if let VariantModal::Filter(m) = &mut self.modal {
                m.parse_error = Some(format!("{e}"));
            }
            let _ = self.scroller.set_filter(None);
            return;
        }
        self.active_filter = Some(arc);
        if let Err(e) = self.rebuild_projection() {
            if let VariantModal::Filter(m) = &mut self.modal {
                m.parse_error = Some(format!("{e}"));
            }
            return;
        }
        self.modal = VariantModal::None;
    }

    fn clear_filter(&mut self) {
        self.active_filter = None;
        let _ = self.scroller.set_filter(None);
        let _ = self.rebuild_projection();
    }

    fn open_filter(&mut self) {
        let buf = self
            .scroller
            .filter_text()
            .map(|s| s.to_string())
            .unwrap_or_default();
        self.modal = VariantModal::Filter(FilterModal { buf, parse_error: None });
    }

    fn open_column_picker(&mut self) {
        let n = self.scroller.schema().fields().len();
        let mut selected = vec![false; n];
        for &i in &self.display_columns {
            if i < n {
                selected[i] = true;
            }
        }
        self.modal = VariantModal::Columns(ColumnsModal { cursor: 0, selected });
    }

    fn close_column_picker(&mut self) {
        if let VariantModal::Columns(picker) = std::mem::replace(&mut self.modal, VariantModal::None) {
            self.display_columns = picker
                .selected
                .iter()
                .enumerate()
                .filter_map(|(i, &k)| if k { Some(i) } else { None })
                .collect();
            if self.display_columns.is_empty() {
                let n = self.scroller.schema().fields().len();
                self.display_columns = (0..n).collect();
            }
            if let Err(e) = self.rebuild_projection() {
                self.error = Some(format!("{e}"));
            }
        }
    }

    fn open_carrier_view(&mut self) {
        let Some((batch, row)) = self.scroller.focused_record() else {
            return;
        };
        let schema = batch.schema();
        let vid = match find_vid(batch, &schema, row) {
            Some(v) => v,
            None => {
                self.modal = VariantModal::Carrier(CarrierModal {
                    vid: String::new(),
                    carriers: None,
                    error: Some(
                        "row has no vid and no (chromosome, position, ref_allele, alt_allele) tuple"
                            .into(),
                    ),
                    scroll: 0,
                });
                return;
            }
        };
        let chrom = chrom_from_path_or_batch(self.scroller.path(), batch, &schema, row)
            .unwrap_or_default();
        let path = self.scroller.path().to_path_buf();
        let mut panel = CarrierModal {
            vid: vid.clone(),
            carriers: None,
            error: None,
            scroll: 0,
        };
        match self.load_carriers(&path, &chrom, &vid) {
            Ok(list) => panel.carriers = Some(list),
            Err(e) => panel.error = Some(format!("{e}")),
        }
        self.modal = VariantModal::Carrier(panel);
    }

    fn close_modal(&mut self) {
        if matches!(self.modal, VariantModal::Columns(_)) {
            self.close_column_picker();
            return;
        }
        self.modal = VariantModal::None;
    }

    fn load_carriers(
        &mut self,
        parquet_path: &Path,
        chrom: &str,
        vid: &str,
    ) -> Result<CarrierList, CohortError> {
        if chrom.is_empty() {
            return Err(CohortError::Input("row has no chromosome".into()));
        }
        let store_dir = locate_store(parquet_path)
            .ok_or_else(|| CohortError::DataMissing(format!(
                "no genotype store found near {}",
                parquet_path.display()
            )))?;
        if self.sample_names.is_none() {
            self.sample_names = Some(load_sample_names(&store_dir)?);
        }
        let chrom_dir = store_dir.join(format!("chromosome={chrom}"));
        if !chrom_dir.exists() {
            return Err(CohortError::DataMissing(format!(
                "store has no chromosome {chrom}"
            )));
        }
        if !self.chrom_cache.contains_key(chrom) {
            let sparse_g = SparseG::open(&chrom_dir)?;
            let index = VariantIndex::load(&chrom_dir)?;
            self.chrom_cache.insert(
                chrom.to_string(),
                ChromCacheEntry { sparse_g, index },
            );
        }
        let entry = self.chrom_cache.get(chrom).unwrap();
        let vvcf = entry
            .index
            .resolve_vid(vid)
            .ok_or_else(|| CohortError::DataMissing(format!("vid not in store: {vid}")))?;
        Ok(entry.sparse_g.load_variant(vvcf))
    }

    fn handle_modal_key(&mut self, key: KeyEvent) -> Transition {
        if matches!(key.code, KeyCode::Esc) {
            self.close_modal();
            return Transition::Stay;
        }
        match &mut self.modal {
            VariantModal::None => Transition::Stay,
            VariantModal::Filter(modal) => {
                match (key.code, key.modifiers) {
                    (KeyCode::Enter, _) => {
                        self.apply_filter_text();
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
                        if !m.contains(KeyModifiers::CONTROL)
                            && !m.contains(KeyModifiers::ALT) =>
                    {
                        modal.buf.push(c);
                        modal.parse_error = None;
                    }
                    _ => {}
                }
                Transition::Stay
            }
            VariantModal::Columns(picker) => {
                let n = picker.selected.len();
                match key.code {
                    KeyCode::Enter => {
                        self.close_column_picker();
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
                Transition::Stay
            }
            VariantModal::Carrier(panel) => {
                match key.code {
                    KeyCode::Char('j') | KeyCode::Down => {
                        if let Some(list) = &panel.carriers {
                            if panel.scroll + 1 < list.len() {
                                panel.scroll += 1;
                            }
                        }
                    }
                    KeyCode::Char('k') | KeyCode::Up => {
                        panel.scroll = panel.scroll.saturating_sub(1);
                    }
                    _ => {}
                }
                Transition::Stay
            }
        }
    }

    fn draw_table(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" {} ", self.header_text()))
            .border_style(Style::default().fg(theme::MUTED));
        let inner = block.inner(area);
        frame.render_widget(block, area);

        let Some((batch, focus)) = self.scroller.focused_record() else {
            let msg = if self.scroller.row_group_count() == 0 {
                "  empty parquet"
            } else {
                "  no rows in this row group (filter excluded all)"
            };
            frame.render_widget(
                Paragraph::new(Line::from(Span::styled(
                    msg,
                    Style::default().fg(theme::MUTED),
                ))),
                inner,
            );
            return;
        };

        let schema = batch.schema();
        let cols: Vec<usize> = self
            .display_columns
            .iter()
            .filter_map(|i| schema.index_of(self.scroller.schema().field(*i).name()).ok())
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
        frame.render_widget(Paragraph::new(lines), inner);
    }

    fn header_text(&self) -> String {
        let total = self.scroller.total_rows_estimate();
        let rg = self.scroller.row_group_count();
        let cur = self.scroller.current_row_group() + 1;
        let visible = self.scroller.current_batch_len();
        let filter = self.scroller.filter_text().unwrap_or("none");
        let base = format!(
            "{} · row group {}/{} · {} visible · {} total est · filter: {}",
            self.title, cur, rg.max(1), visible, total, filter
        );
        match <Self as ArtifactView>::header(self) {
            Some(s) => format!("{base} · {s}"),
            None => base,
        }
    }

    fn draw_detail(&self, frame: &mut Frame, area: Rect) {
        let inner = area;
        let Some((batch, row)) = self.scroller.focused_record() else {
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
        frame.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_filter_bar(&self, frame: &mut Frame, area: Rect, modal: &FilterModal) {
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
        frame.render_widget(Paragraph::new(vec![input, hint]), area);
    }

    fn draw_carrier_panel(&self, frame: &mut Frame, area: Rect, panel: &CarrierModal) {
        if let Some(err) = &panel.error {
            frame.render_widget(
                Paragraph::new(Line::from(Span::styled(
                    format!(" carriers: {err}"),
                    Style::default().fg(theme::BAD),
                ))),
                area,
            );
            return;
        }
        let Some(list) = &panel.carriers else {
            return;
        };
        let inner = area;
        let header = format!(
            "  {} carriers of {}",
            list.len(),
            if panel.vid.is_empty() { "<row>" } else { panel.vid.as_str() }
        );
        let mut lines: Vec<Line> = Vec::new();
        lines.push(Line::from(Span::styled(
            header,
            Style::default().fg(theme::ACCENT).bold(),
        )));
        let visible = inner.height.saturating_sub(1) as usize;
        let names = self.sample_names.as_deref().unwrap_or(&[]);
        if panel.scroll > 0 {
            lines.push(Line::from(Span::styled(
                format!("    ↑ {} above", panel.scroll),
                Style::default().fg(theme::MUTED),
            )));
        }
        let body_cap = visible.saturating_sub(usize::from(panel.scroll > 0));
        let total = list.entries.len();
        let end = (panel.scroll + body_cap).min(total);
        let below = total.saturating_sub(end);
        let body_take = body_cap.saturating_sub(usize::from(below > 0));
        for entry in list
            .entries
            .iter()
            .skip(panel.scroll)
            .take(body_take)
        {
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
        frame.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_column_picker(&self, frame: &mut Frame, area: Rect, picker: &ColumnsModal) {
        let overlay = centered(area, 50, 70);
        frame.render_widget(Clear, overlay);
        let inner = overlay;
        let schema = self.scroller.schema();
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
        frame.render_stateful_widget(List::new(items), inner, &mut state);
    }
}

impl Screen for VariantScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, _log: &LogTail) {
        let modal_height: u16 = match &self.modal {
            VariantModal::Filter(_) => 2,
            VariantModal::Carrier(_) => 10,
            _ => 0,
        };
        let detail_height: u16 = if area.height >= 24 { 8 } else { 0 };
        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(6),
                Constraint::Length(detail_height),
                Constraint::Length(modal_height),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        self.draw_table(frame, v[0]);
        if detail_height > 0 && matches!(self.modal, VariantModal::None) {
            self.draw_detail(frame, v[1]);
        }
        match &self.modal {
            VariantModal::Filter(m) => self.draw_filter_bar(frame, v[2], m),
            VariantModal::Carrier(p) => self.draw_carrier_panel(frame, v[2], p),
            _ => {}
        }

        let error_line = self
            .error
            .clone()
            .or_else(|| match &self.modal {
                VariantModal::Filter(m) => m.parse_error.clone().map(|e| format!("error: {e}")),
                _ => None,
            })
            .unwrap_or_default();
        frame.render_widget(
            Paragraph::new(Line::from(Span::styled(
                format!(" {error_line}"),
                Style::default().fg(theme::BAD),
            ))),
            v[3],
        );

        let hint = match &self.modal {
            VariantModal::None => {
                "/ filter  ! columns  c carriers  { } page  g G start/end  q back"
            }
            VariantModal::Filter(_) => "enter apply  esc cancel  ctrl-u clear",
            VariantModal::Columns(_) => "space toggle  enter close  esc cancel",
            VariantModal::Carrier(_) => "j k scroll  esc close",
        };
        StatusBar { title: &self.title, keys: hint }.render(frame, v[4]);

        if let VariantModal::Columns(p) = &self.modal {
            self.draw_column_picker(frame, area, p);
        }
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Variant
    }

    fn keys(&self) -> KeyMap {
        if !matches!(self.modal, VariantModal::None) {
            return KeyMap::new();
        }
        let none = KeyModifiers::NONE;
        KeyMap::new()
            .bind(KeyCode::Char('q'), none, Action::VariantClose)
            .bind(KeyCode::Esc, none, Action::VariantClose)
            .bind(KeyCode::Char('j'), none, Action::VariantScrollRowDown)
            .bind(KeyCode::Down, none, Action::VariantScrollRowDown)
            .bind(KeyCode::Char('k'), none, Action::VariantScrollRowUp)
            .bind(KeyCode::Up, none, Action::VariantScrollRowUp)
            .bind(KeyCode::Char('{'), none, Action::VariantPrevRowGroup)
            .bind(KeyCode::Char('}'), none, Action::VariantNextRowGroup)
            .bind(KeyCode::Char('g'), none, Action::VariantJumpStart)
            .bind(KeyCode::Char('G'), KeyModifiers::SHIFT, Action::VariantJumpEnd)
            .bind(KeyCode::Char('/'), none, Action::VariantOpenFilter)
            .bind(KeyCode::Char('x'), none, Action::VariantFilterClear)
            .bind(KeyCode::Char('!'), none, Action::VariantOpenColumnPicker)
            .bind(KeyCode::Char('c'), none, Action::VariantOpenCarrierView)
            .bind(KeyCode::Enter, none, Action::VariantOpenLinked)
    }

    fn handle(&mut self, event: &AppEvent) -> Transition {
        let AppEvent::Key(k) = event else {
            return Transition::Stay;
        };
        if matches!(self.modal, VariantModal::None) {
            return match self.keys().lookup(k.code, k.modifiers) {
                Some(action) => self.on_action(action),
                None => Transition::Stay,
            };
        }
        self.handle_modal_key(*k)
    }

    fn on_action(&mut self, action: Action) -> Transition {
        self.error = None;
        match action {
            Action::VariantClose => Transition::Pop,
            Action::VariantScrollRowDown => {
                self.scroller.scroll_down();
                Transition::Stay
            }
            Action::VariantScrollRowUp => {
                self.scroller.scroll_up();
                Transition::Stay
            }
            Action::VariantPrevRowGroup => {
                if let Err(e) = self.scroller.prev_row_group() {
                    self.error = Some(format!("{e}"));
                }
                Transition::Stay
            }
            Action::VariantNextRowGroup => {
                if let Err(e) = self.scroller.next_row_group() {
                    self.error = Some(format!("{e}"));
                }
                Transition::Stay
            }
            Action::VariantJumpStart => {
                if let Err(e) = self.scroller.goto_row_group(0) {
                    self.error = Some(format!("{e}"));
                }
                Transition::Stay
            }
            Action::VariantJumpEnd => {
                let last = self.scroller.row_group_count().saturating_sub(1);
                if let Err(e) = self.scroller.goto_row_group(last) {
                    self.error = Some(format!("{e}"));
                }
                Transition::Stay
            }
            Action::VariantOpenFilter => {
                self.open_filter();
                Transition::Stay
            }
            Action::VariantFilterClear => {
                self.clear_filter();
                Transition::Stay
            }
            Action::VariantOpenColumnPicker => {
                self.open_column_picker();
                Transition::Stay
            }
            Action::VariantOpenCarrierView => {
                self.open_carrier_view();
                Transition::Stay
            }
            Action::VariantOpenLinked => {
                self.try_cross_open();
                if let ArtifactOutcome::OpenArtifact { path, filter } =
                    <Self as ArtifactView>::outcome(self)
                {
                    let companion = self.companion.clone();
                    let prefix = path
                        .file_name()
                        .map(|n| n.to_string_lossy().into_owned());
                    let opened = Self::open(path, prefix.as_deref(), None, companion).and_then(
                        |s| match filter {
                            Some(f) => {
                                let text = match &f.predicate {
                                    Predicate::Eq(v) => {
                                        format!("{} = \"{}\"", f.column.as_str(), v)
                                    }
                                };
                                s.with_initial_filter(&text)
                            }
                            None => Ok(s),
                        },
                    );
                    match opened {
                        Ok(screen) => return Transition::Push(Box::new(screen)),
                        Err(e) => self.error = Some(format!("{e}")),
                    }
                }
                Transition::Stay
            }
            Action::VariantCloseCarrierView | Action::VariantColumnPickerClose => {
                self.close_modal();
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }
}

impl ArtifactView for VariantScreen {
    fn header(&self) -> Option<String> {
        self.summary.clone()
    }

    fn outcome(&mut self) -> ArtifactOutcome {
        match self.pending_open.take() {
            Some((path, filter)) => ArtifactOutcome::OpenArtifact { path, filter },
            None => ArtifactOutcome::Continue,
        }
    }
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
    let entries = std::fs::read_dir(set_root).map_err(|e| {
        CohortError::Resource(format!("read {}: {e}", set_root.display()))
    })?;
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
    match array.data_type() {
        DataType::Int32 => array
            .as_any()
            .downcast_ref::<Int32Array>()
            .map(|a| a.value(row).to_string())
            .unwrap_or_default(),
        DataType::Int64 => array
            .as_any()
            .downcast_ref::<Int64Array>()
            .map(|a| a.value(row).to_string())
            .unwrap_or_default(),
        DataType::UInt32 => array
            .as_any()
            .downcast_ref::<UInt32Array>()
            .map(|a| a.value(row).to_string())
            .unwrap_or_default(),
        DataType::UInt64 => array
            .as_any()
            .downcast_ref::<UInt64Array>()
            .map(|a| a.value(row).to_string())
            .unwrap_or_default(),
        DataType::Float32 => array
            .as_any()
            .downcast_ref::<Float32Array>()
            .map(|a| format!("{:.4}", a.value(row)))
            .unwrap_or_default(),
        DataType::Float64 => array
            .as_any()
            .downcast_ref::<Float64Array>()
            .map(|a| format!("{:.4}", a.value(row)))
            .unwrap_or_default(),
        DataType::Utf8 | DataType::LargeUtf8 | DataType::Utf8View => array
            .as_any()
            .downcast_ref::<StringArray>()
            .map(|a| a.value(row).to_string())
            .unwrap_or_default(),
        DataType::Boolean => array
            .as_any()
            .downcast_ref::<BooleanArray>()
            .map(|a| if a.value(row) { "true" } else { "false" }.to_string())
            .unwrap_or_default(),
        other => format!("{other:?}"),
    }
}

fn centered(area: Rect, pct_w: u16, pct_h: u16) -> Rect {
    let w = area.width * pct_w / 100;
    let h = area.height * pct_h / 100;
    let x = area.x + (area.width.saturating_sub(w)) / 2;
    let y = area.y + (area.height.saturating_sub(h)) / 2;
    Rect { x, y, width: w, height: h }
}
