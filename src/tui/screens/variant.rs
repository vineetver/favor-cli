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
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::state::arrow_predicate::CompiledFilter;
use crate::tui::state::parquet_scroller::{ParquetScroller, RowFilterFactory};
use crate::tui::theme;
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

pub struct VariantScreen {
    title: String,
    scroller: ParquetScroller,
    display_columns: Vec<usize>,
    filter_bar: FilterBarState,
    column_picker: Option<ColumnPickerState>,
    carrier: Option<CarrierPanelState>,
    chrom_cache: HashMap<String, ChromCacheEntry>,
    sample_names: Option<Vec<String>>,
    error: Option<String>,
    sub_mode: SubMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SubMode {
    Browse,
    FilterEdit,
    ColumnPicker,
    Carrier,
}

struct FilterBarState {
    buf: String,
    compiled: Option<Arc<CompiledFilter>>,
    parse_error: Option<String>,
}

struct ColumnPickerState {
    cursor: usize,
    selected: Vec<bool>,
}

struct CarrierPanelState {
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
    pub fn new_for_parquet(path: PathBuf) -> Result<Self, CohortError> {
        let scroller = ParquetScroller::open(&path)?;
        let display_columns: Vec<usize> = (0..scroller.schema().fields().len()).collect();
        let title = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_else(|| path.display().to_string());
        Ok(Self {
            title,
            scroller,
            display_columns,
            filter_bar: FilterBarState {
                buf: String::new(),
                compiled: None,
                parse_error: None,
            },
            column_picker: None,
            carrier: None,
            chrom_cache: HashMap::new(),
            sample_names: None,
            error: None,
            sub_mode: SubMode::Browse,
        })
    }

    pub fn new_for_annotated_set(set_path: PathBuf) -> Result<Self, CohortError> {
        let parquet = pick_first_chrom_parquet(&set_path)?;
        let mut s = Self::new_for_parquet(parquet)?;
        s.title = format!(
            "{} :: {}",
            set_path
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
            s.title
        );
        Ok(s)
    }

    fn rebuild_projection(&mut self) -> Result<(), CohortError> {
        let n_fields = self.scroller.schema().fields().len();
        let mut wanted: Vec<bool> = vec![false; n_fields];
        for &i in &self.display_columns {
            if i < n_fields {
                wanted[i] = true;
            }
        }
        if let Some(filter) = &self.filter_bar.compiled {
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
            leaves.into_iter(),
        );
        self.scroller.set_projection(mask)
    }

    fn apply_filter_text(&mut self) {
        let text = self.filter_bar.buf.trim().to_string();
        if text.is_empty() {
            self.filter_bar.compiled = None;
            self.filter_bar.parse_error = None;
            if let Err(e) = self
                .scroller
                .set_filter(None)
                .and_then(|_| self.rebuild_projection())
            {
                self.error = Some(format!("{e}"));
            }
            self.sub_mode = SubMode::Browse;
            return;
        }
        let parsed = match CompiledFilter::parse(&text) {
            Ok(c) => c,
            Err(e) => {
                self.filter_bar.parse_error = Some(e);
                return;
            }
        };
        let arc = Arc::new(parsed);
        self.filter_bar.compiled = Some(arc.clone());
        let factory: Arc<dyn RowFilterFactory> = arc;
        if let Err(e) = self.scroller.set_filter(Some(factory)) {
            self.filter_bar.parse_error = Some(format!("{e}"));
            self.filter_bar.compiled = None;
            let _ = self.scroller.set_filter(None);
            return;
        }
        if let Err(e) = self.rebuild_projection() {
            self.filter_bar.parse_error = Some(format!("{e}"));
            return;
        }
        self.filter_bar.parse_error = None;
        self.sub_mode = SubMode::Browse;
    }

    fn clear_filter(&mut self) {
        self.filter_bar.buf.clear();
        self.filter_bar.compiled = None;
        self.filter_bar.parse_error = None;
        let _ = self.scroller.set_filter(None);
        let _ = self.rebuild_projection();
    }

    fn open_column_picker(&mut self) {
        let n = self.scroller.schema().fields().len();
        let mut selected = vec![false; n];
        for &i in &self.display_columns {
            if i < n {
                selected[i] = true;
            }
        }
        self.column_picker = Some(ColumnPickerState { cursor: 0, selected });
        self.sub_mode = SubMode::ColumnPicker;
    }

    fn close_column_picker(&mut self) {
        if let Some(picker) = self.column_picker.take() {
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
        self.sub_mode = SubMode::Browse;
    }

    fn open_carrier_view(&mut self) {
        let Some((batch, row)) = self.scroller.focused_record() else {
            return;
        };
        let schema = batch.schema();
        let vid = match find_vid(batch, &schema, row) {
            Some(v) => v,
            None => {
                self.carrier = Some(CarrierPanelState {
                    vid: String::new(),
                    carriers: None,
                    error: Some(
                        "row has no vid and no (chromosome, position, ref_allele, alt_allele) tuple"
                            .into(),
                    ),
                    scroll: 0,
                });
                self.sub_mode = SubMode::Carrier;
                return;
            }
        };
        let chrom = chrom_from_path_or_batch(self.scroller.path(), batch, &schema, row)
            .unwrap_or_default();
        let path = self.scroller.path().to_path_buf();
        let mut panel = CarrierPanelState {
            vid: vid.clone(),
            carriers: None,
            error: None,
            scroll: 0,
        };
        match self.load_carriers(&path, &chrom, &vid) {
            Ok(list) => panel.carriers = Some(list),
            Err(e) => panel.error = Some(format!("{e}")),
        }
        self.carrier = Some(panel);
        self.sub_mode = SubMode::Carrier;
    }

    fn close_carrier_view(&mut self) {
        self.carrier = None;
        self.sub_mode = SubMode::Browse;
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

    fn handle_filter_edit_key(&mut self, key: KeyEvent) -> Transition {
        match (key.code, key.modifiers) {
            (KeyCode::Esc, _) => {
                self.filter_bar.parse_error = None;
                self.sub_mode = SubMode::Browse;
                Transition::Stay
            }
            (KeyCode::Enter, _) => {
                self.apply_filter_text();
                Transition::Stay
            }
            (KeyCode::Backspace, _) => {
                self.filter_bar.buf.pop();
                self.filter_bar.parse_error = None;
                Transition::Stay
            }
            (KeyCode::Char('u'), m) if m.contains(KeyModifiers::CONTROL) => {
                self.filter_bar.buf.clear();
                self.filter_bar.parse_error = None;
                Transition::Stay
            }
            (KeyCode::Char(c), m) if !m.contains(KeyModifiers::CONTROL) && !m.contains(KeyModifiers::ALT) => {
                self.filter_bar.buf.push(c);
                self.filter_bar.parse_error = None;
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }

    fn handle_column_picker_key(&mut self, key: KeyEvent) -> Transition {
        let Some(picker) = self.column_picker.as_mut() else {
            self.sub_mode = SubMode::Browse;
            return Transition::Stay;
        };
        let n = picker.selected.len();
        match key.code {
            KeyCode::Esc | KeyCode::Enter => {
                self.close_column_picker();
                Transition::Stay
            }
            KeyCode::Char('j') | KeyCode::Down => {
                if n > 0 && picker.cursor + 1 < n {
                    picker.cursor += 1;
                }
                Transition::Stay
            }
            KeyCode::Char('k') | KeyCode::Up => {
                if picker.cursor > 0 {
                    picker.cursor -= 1;
                }
                Transition::Stay
            }
            KeyCode::Char(' ') => {
                if picker.cursor < n {
                    picker.selected[picker.cursor] = !picker.selected[picker.cursor];
                }
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }

    fn handle_carrier_key(&mut self, key: KeyEvent) -> Transition {
        match key.code {
            KeyCode::Esc => {
                self.close_carrier_view();
                Transition::Stay
            }
            KeyCode::Char('j') | KeyCode::Down => {
                if let Some(panel) = self.carrier.as_mut() {
                    if let Some(list) = &panel.carriers {
                        if panel.scroll + 1 < list.len() {
                            panel.scroll += 1;
                        }
                    }
                }
                Transition::Stay
            }
            KeyCode::Char('k') | KeyCode::Up => {
                if let Some(panel) = self.carrier.as_mut() {
                    panel.scroll = panel.scroll.saturating_sub(1);
                }
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }

    fn draw_table(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(format!(" {} ", self.header_text()))
            .border_style(Style::default().fg(theme::ACCENT));
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
        let widths = compute_widths(batch, &cols, inner.width as usize);
        let header_line = Line::from(
            cols.iter()
                .zip(widths.iter())
                .map(|(c, w)| {
                    Span::styled(
                        pad(&schema.field(*c).name(), *w),
                        Style::default().fg(theme::ACCENT).bold(),
                    )
                })
                .collect::<Vec<_>>(),
        );

        let visible_rows = inner.height.saturating_sub(1) as usize;
        let total = batch.num_rows();
        let start = focus.saturating_sub(visible_rows / 2).min(total.saturating_sub(visible_rows.max(1)));
        let end = (start + visible_rows).min(total);

        let mut lines: Vec<Line> = Vec::with_capacity(end - start + 1);
        lines.push(header_line);
        for row in start..end {
            let style = if row == focus {
                Style::default().bg(theme::MUTED).fg(theme::FG).bold()
            } else {
                Style::default().fg(theme::FG)
            };
            let spans: Vec<Span> = cols
                .iter()
                .zip(widths.iter())
                .map(|(c, w)| Span::styled(pad(&format_cell(batch, *c, row), *w), style))
                .collect();
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
        format!(
            "{} · row group {}/{} · {} visible · {} total est · filter: {}",
            self.title, cur, rg.max(1), visible, total, filter
        )
    }

    fn draw_detail(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(" Detail ")
            .border_style(Style::default().fg(theme::MUTED));
        let inner = block.inner(area);
        frame.render_widget(block, area);

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

        let half = (pairs.len() + 1) / 2;
        let max_rows = inner.height as usize;
        let mut lines: Vec<Line> = Vec::with_capacity(max_rows);
        for i in 0..half.min(max_rows) {
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
        frame.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_filter_bar(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(" Filter ")
            .border_style(Style::default().fg(theme::ACCENT));
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let mut lines = vec![Line::from(vec![
            Span::styled("/ ", Style::default().fg(theme::WARN).bold()),
            Span::styled(self.filter_bar.buf.as_str(), Style::default().fg(theme::FG)),
            Span::styled("_", Style::default().fg(theme::ACCENT)),
        ])];
        if let Some(err) = &self.filter_bar.parse_error {
            lines.push(Line::from(Span::styled(
                format!("  error: {err}"),
                Style::default().fg(theme::BAD),
            )));
        }
        frame.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_carrier_panel(&self, frame: &mut Frame, area: Rect) {
        let block = Block::default()
            .borders(Borders::ALL)
            .title(" Carriers ")
            .border_style(Style::default().fg(theme::ACCENT));
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let Some(panel) = &self.carrier else {
            return;
        };
        if let Some(err) = &panel.error {
            frame.render_widget(
                Paragraph::new(Line::from(Span::styled(
                    format!("  {err}"),
                    Style::default().fg(theme::BAD),
                ))),
                inner,
            );
            return;
        }
        let Some(list) = &panel.carriers else {
            return;
        };
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
        for entry in list
            .entries
            .iter()
            .skip(panel.scroll)
            .take(visible)
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
        frame.render_widget(Paragraph::new(lines), inner);
    }

    fn draw_column_picker(&self, frame: &mut Frame, area: Rect) {
        let Some(picker) = &self.column_picker else {
            return;
        };
        let overlay = centered(area, 50, 70);
        frame.render_widget(Clear, overlay);
        let block = Block::default()
            .borders(Borders::ALL)
            .title(" Columns — space toggle, enter close ")
            .border_style(Style::default().fg(theme::ACCENT));
        let inner = block.inner(overlay);
        frame.render_widget(block, overlay);
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
        let bottom_height: u16 = match self.sub_mode {
            SubMode::FilterEdit => 4,
            SubMode::Carrier => 10,
            _ => 1,
        };
        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(6),
                Constraint::Length(8),
                Constraint::Length(bottom_height),
                Constraint::Length(1),
            ])
            .split(area);

        self.draw_table(frame, v[0]);
        self.draw_detail(frame, v[1]);
        match self.sub_mode {
            SubMode::FilterEdit => self.draw_filter_bar(frame, v[2]),
            SubMode::Carrier => self.draw_carrier_panel(frame, v[2]),
            _ => {}
        }

        let status = self
            .error
            .clone()
            .unwrap_or_else(|| "/ filter  ! columns  c carriers  { } page  g G start/end  q back".into());
        StatusBar { title: &self.title, keys: &status }.render(frame, v[3]);

        if self.sub_mode == SubMode::ColumnPicker {
            self.draw_column_picker(frame, area);
        }
    }

    fn scope(&self) -> ActionScope {
        ActionScope::Variant
    }

    fn keys(&self) -> KeyMap {
        if self.sub_mode != SubMode::Browse {
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
    }

    /// Sub-modes (filter edit, column picker, carrier panel) bypass the
    /// scoped KeyMap so raw keystrokes go directly to the active sub-widget.
    /// Same correctness invariant as SetupScreen::handle.
    fn handle(&mut self, event: &AppEvent) -> Transition {
        let AppEvent::Key(k) = event else {
            return Transition::Stay;
        };
        match self.sub_mode {
            SubMode::FilterEdit => self.handle_filter_edit_key(*k),
            SubMode::ColumnPicker => self.handle_column_picker_key(*k),
            SubMode::Carrier => self.handle_carrier_key(*k),
            SubMode::Browse => match self.keys().lookup(k.code, k.modifiers) {
                Some(action) => self.on_action(action),
                None => Transition::Stay,
            },
        }
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
                self.sub_mode = SubMode::FilterEdit;
                self.filter_bar.parse_error = None;
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
            Action::VariantCloseCarrierView => {
                self.close_carrier_view();
                Transition::Stay
            }
            Action::VariantColumnPickerClose => {
                self.close_column_picker();
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }
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
