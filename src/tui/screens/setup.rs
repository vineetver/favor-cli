use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph};
use ratatui::Frame;

use crate::config::{Config, Environment, ResourceConfig, Tier};
use crate::data::Pack;
use crate::resource::Resources;
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::theme::{self, Tone, FOCUS_GLYPH};
use crate::tui::widgets::file_picker::{self, tab_complete, DirBrowserState};
use crate::tui::widgets::log_tail::LogTail;

const TIERS: [(Tier, &str); 2] = [
    (Tier::Base, "curated annotations for most analyses"),
    (Tier::Full, "complete FAVOR database, all annotations"),
];

const MEMORY_PRESETS: &[(&str, u64)] = &[
    ("8GB", 8 * 1024 * 1024 * 1024),
    ("16GB", 16 * 1024 * 1024 * 1024),
    ("32GB", 32 * 1024 * 1024 * 1024),
    ("64GB", 64 * 1024 * 1024 * 1024),
    ("128GB", 128 * 1024 * 1024 * 1024),
    ("256GB", 256 * 1024 * 1024 * 1024),
];

#[derive(Debug, Clone)]
pub struct SetupOutcome {
    pub tier: Tier,
    pub root: PathBuf,
    pub packs: Vec<String>,
    pub environment: Option<Environment>,
    pub memory_budget: Option<String>,
}

pub type OutcomeSink = Arc<Mutex<Option<SetupOutcome>>>;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Field {
    Tier,
    Root,
    Advanced,
    Env,
    Memory,
    Packs,
    Save,
}

enum Mode {
    Form,
    BrowseRoot(DirBrowserState),
    EditMemory(String),
    EditPacks(usize, Vec<bool>),
}

pub struct SetupScreen {
    mode: Mode,
    focus: Field,
    advanced_open: bool,
    tier: Tier,
    root: PathBuf,
    packs: Vec<String>,
    environment: Option<Environment>,
    memory_budget: Option<String>,
    resources: Resources,
    sink: Option<OutcomeSink>,
    error: Option<String>,
}

impl SetupScreen {
    pub fn new() -> Self {
        let cfg = Config::load().unwrap_or_default();
        let root = if cfg.data.root_dir.is_empty() {
            Config::default_root_dir()
        } else {
            PathBuf::from(&cfg.data.root_dir)
        };
        Self {
            mode: Mode::Form,
            focus: Field::Tier,
            advanced_open: false,
            tier: cfg.data.tier,
            root,
            packs: Vec::new(),
            environment: cfg.resources.environment,
            memory_budget: cfg.resources.memory_budget.clone(),
            resources: Resources::detect(),
            sink: None,
            error: None,
        }
    }

    pub fn with_sink(
        env: Option<Environment>,
        memory_budget: Option<String>,
        sink: OutcomeSink,
    ) -> Self {
        let mut s = Self::new();
        if env.is_some() {
            s.environment = env;
        }
        if memory_budget.is_some() {
            s.memory_budget = memory_budget;
        }
        s.sink = Some(sink);
        s
    }

    fn visible_fields(&self) -> Vec<Field> {
        let mut v = vec![Field::Tier, Field::Root, Field::Advanced];
        if self.advanced_open {
            v.push(Field::Env);
            v.push(Field::Memory);
            v.push(Field::Packs);
        }
        v.push(Field::Save);
        v
    }

    fn focus_index(&self) -> usize {
        self.visible_fields()
            .iter()
            .position(|f| *f == self.focus)
            .unwrap_or(0)
    }

    fn focus_next(&mut self) {
        let v = self.visible_fields();
        let i = (self.focus_index() + 1) % v.len();
        self.focus = v[i];
    }

    fn focus_prev(&mut self) {
        let v = self.visible_fields();
        let n = v.len();
        let i = (self.focus_index() + n - 1) % n;
        self.focus = v[i];
    }

    fn save(&mut self) -> Transition {
        if !self.root.is_dir() {
            self.error = Some(format!("data root not a directory: {}", self.root.display()));
            return Transition::Stay;
        }
        if let Some(sink) = self.sink.as_ref() {
            *sink.lock().unwrap() = Some(SetupOutcome {
                tier: self.tier,
                root: self.root.clone(),
                packs: self.packs.clone(),
                environment: self.environment,
                memory_budget: self.memory_budget.clone(),
            });
        }
        Transition::Pop
    }

    fn open_browse(&mut self) {
        let start = if self.root.is_dir() {
            self.root.clone()
        } else {
            std::env::current_dir().unwrap_or_else(|_| Config::default_root_dir())
        };
        self.mode = Mode::BrowseRoot(DirBrowserState::new("Select FAVOR data root", &start));
    }

    fn open_edit_packs(&mut self) {
        let optional = Pack::optional();
        let checked: Vec<bool> = optional
            .iter()
            .map(|p| self.packs.iter().any(|id| id == p.id))
            .collect();
        self.mode = Mode::EditPacks(0, checked);
    }

    fn open_edit_memory(&mut self) {
        let buf = self.memory_budget.clone().unwrap_or_default();
        self.mode = Mode::EditMemory(buf);
    }

    fn cycle_tier(&mut self) {
        self.tier = match self.tier {
            Tier::Base => Tier::Full,
            Tier::Full => Tier::Base,
        };
    }

    fn cycle_env(&mut self) {
        self.environment = match self.environment {
            None => Some(Environment::Hpc),
            Some(Environment::Hpc) => Some(Environment::Workstation),
            Some(Environment::Workstation) => None,
        };
    }

    fn cycle_memory(&mut self) {
        let cur = self.memory_budget.as_deref();
        let idx = MEMORY_PRESETS.iter().position(|(v, _)| Some(*v) == cur);
        let next = match idx {
            None => 0,
            Some(i) if i + 1 < MEMORY_PRESETS.len() => i + 1,
            Some(_) => {
                self.memory_budget = None;
                return;
            }
        };
        self.memory_budget = Some(MEMORY_PRESETS[next].0.to_string());
    }

    fn handle_form(&mut self, code: KeyCode) -> Transition {
        self.error = None;
        match code {
            KeyCode::Esc | KeyCode::Char('q') => return Transition::Pop,
            KeyCode::Up | KeyCode::Char('k') => self.focus_prev(),
            KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => self.focus_next(),
            KeyCode::Left | KeyCode::Char('h') => match self.focus {
                Field::Tier => self.cycle_tier(),
                Field::Env => self.cycle_env(),
                Field::Memory => self.cycle_memory(),
                _ => {}
            },
            KeyCode::Right | KeyCode::Char('l') => match self.focus {
                Field::Tier => self.cycle_tier(),
                Field::Env => self.cycle_env(),
                Field::Memory => self.cycle_memory(),
                _ => {}
            },
            KeyCode::Char(' ') => match self.focus {
                Field::Tier => self.cycle_tier(),
                Field::Env => self.cycle_env(),
                Field::Memory => self.cycle_memory(),
                Field::Advanced => self.advanced_open = !self.advanced_open,
                Field::Root => self.open_browse(),
                Field::Packs => self.open_edit_packs(),
                Field::Save => return self.save(),
            },
            KeyCode::Enter => match self.focus {
                Field::Tier => self.cycle_tier(),
                Field::Root => self.open_browse(),
                Field::Advanced => self.advanced_open = !self.advanced_open,
                Field::Env => self.cycle_env(),
                Field::Memory => self.open_edit_memory(),
                Field::Packs => self.open_edit_packs(),
                Field::Save => return self.save(),
            },
            _ => {}
        }
        Transition::Stay
    }

    fn handle_browse(&mut self, code: KeyCode) -> Transition {
        let Mode::BrowseRoot(state) = &mut self.mode else {
            return Transition::Stay;
        };
        if state.typing_path {
            match code {
                KeyCode::Enter => {
                    let p = PathBuf::from(&state.input_buf);
                    if p.is_dir() {
                        state.navigate_to(p);
                        state.typing_path = false;
                    }
                }
                KeyCode::Esc => {
                    state.input_buf = state.current_dir.to_string_lossy().to_string();
                    state.typing_path = false;
                }
                KeyCode::Backspace => {
                    state.input_buf.pop();
                }
                KeyCode::Char(c) => state.input_buf.push(c),
                KeyCode::Tab => {
                    if let Some(c) = tab_complete(&state.input_buf) {
                        state.input_buf = c;
                    }
                }
                _ => {}
            }
            let typed = std::path::Path::new(&state.input_buf);
            if typed.is_dir() {
                state.probe = crate::config::DirProbe::scan(typed);
            }
            return Transition::Stay;
        }
        match code {
            KeyCode::Up | KeyCode::Char('k') => state.select_up(),
            KeyCode::Down | KeyCode::Char('j') => state.select_down(),
            KeyCode::Right => {
                if let Some(sel) = state.enter_selected() {
                    self.root = sel;
                    self.mode = Mode::Form;
                }
            }
            KeyCode::Left | KeyCode::Backspace => state.go_parent(),
            KeyCode::Enter | KeyCode::Char(' ') => {
                self.root = state.current_dir.clone();
                self.mode = Mode::Form;
            }
            KeyCode::Char('c') => {
                let _ = std::fs::create_dir_all(&state.current_dir);
                self.root = state.current_dir.clone();
                self.mode = Mode::Form;
            }
            KeyCode::Char('/') | KeyCode::Char('g') => {
                state.typing_path = true;
                state.input_buf = state.current_dir.to_string_lossy().to_string();
            }
            KeyCode::Esc | KeyCode::Char('q') => self.mode = Mode::Form,
            _ => {}
        }
        Transition::Stay
    }

    fn handle_memory_edit(&mut self, code: KeyCode) -> Transition {
        let Mode::EditMemory(buf) = &mut self.mode else {
            return Transition::Stay;
        };
        match code {
            KeyCode::Enter => {
                if buf.is_empty() {
                    self.memory_budget = None;
                    self.mode = Mode::Form;
                } else if ResourceConfig::parse_memory_bytes(buf).is_some() {
                    self.memory_budget = Some(buf.clone());
                    self.mode = Mode::Form;
                } else {
                    self.error = Some(format!("invalid memory: {buf} (e.g. 48GB, 12288MB)"));
                }
            }
            KeyCode::Esc => self.mode = Mode::Form,
            KeyCode::Backspace => {
                buf.pop();
            }
            KeyCode::Char(c) => buf.push(c),
            _ => {}
        }
        Transition::Stay
    }

    fn handle_packs_edit(&mut self, code: KeyCode) -> Transition {
        let Mode::EditPacks(cursor, checked) = &mut self.mode else {
            return Transition::Stay;
        };
        let n = Pack::optional().len();
        match code {
            KeyCode::Up | KeyCode::Char('k') => *cursor = cursor.saturating_sub(1),
            KeyCode::Down | KeyCode::Char('j') => {
                if *cursor + 1 < n {
                    *cursor += 1;
                }
            }
            KeyCode::Char(' ') => {
                checked[*cursor] = !checked[*cursor];
            }
            KeyCode::Char('a') => {
                let all = checked.iter().all(|&c| c);
                for c in checked.iter_mut() {
                    *c = !all;
                }
            }
            KeyCode::Enter => {
                self.packs = Pack::optional()
                    .into_iter()
                    .zip(checked.iter())
                    .filter(|(_, &on)| on)
                    .map(|(p, _)| p.id.to_string())
                    .collect();
                self.mode = Mode::Form;
            }
            KeyCode::Esc | KeyCode::Char('q') => self.mode = Mode::Form,
            _ => {}
        }
        Transition::Stay
    }

    fn glyph(&self, f: Field) -> &'static str {
        if self.focus == f {
            FOCUS_GLYPH
        } else {
            " "
        }
    }

    fn field_line(&self, f: Field) -> Line<'static> {
        let g = self.glyph(f);
        let label_tone = if self.focus == f { Tone::Focus } else { Tone::Normal };
        match f {
            Field::Tier => {
                let (_, summary) = TIERS.iter().find(|(t, _)| *t == self.tier).unwrap();
                Line::from(vec![
                    Span::styled(format!(" {g} "), label_tone.style()),
                    Span::styled("tier      ", label_tone.style()),
                    Span::styled(
                        format!("{} {}", self.tier.as_str(), self.tier.size_human()),
                        Tone::Warn.style(),
                    ),
                    Span::styled(format!("  {summary}"), Tone::Muted.style()),
                ])
            }
            Field::Root => {
                let exists = self.root.is_dir();
                let path_tone = if exists { Tone::Normal } else { Tone::Bad };
                Line::from(vec![
                    Span::styled(format!(" {g} "), label_tone.style()),
                    Span::styled("data root ", label_tone.style()),
                    Span::styled(self.root.display().to_string(), path_tone.style()),
                ])
            }
            Field::Advanced => {
                let mark = if self.advanced_open { "[-]" } else { "[+]" };
                Line::from(vec![
                    Span::styled(format!(" {g} "), label_tone.style()),
                    Span::styled(format!("{mark} advanced"), label_tone.style()),
                ])
            }
            Field::Env => {
                let val = match self.environment {
                    None => "auto".to_string(),
                    Some(Environment::Hpc) => "hpc".to_string(),
                    Some(Environment::Workstation) => "workstation".to_string(),
                };
                Line::from(vec![
                    Span::styled(format!("   {g} "), label_tone.style()),
                    Span::styled("env       ", label_tone.style()),
                    Span::styled(val, Tone::Warn.style()),
                ])
            }
            Field::Memory => {
                let val = self
                    .memory_budget
                    .clone()
                    .unwrap_or_else(|| format!("auto ({})", self.resources.memory_human()));
                Line::from(vec![
                    Span::styled(format!("   {g} "), label_tone.style()),
                    Span::styled("memory    ", label_tone.style()),
                    Span::styled(val, Tone::Warn.style()),
                ])
            }
            Field::Packs => {
                let val = if self.packs.is_empty() {
                    "none".to_string()
                } else {
                    self.packs.join(", ")
                };
                Line::from(vec![
                    Span::styled(format!("   {g} "), label_tone.style()),
                    Span::styled("packs     ", label_tone.style()),
                    Span::styled(val, Tone::Warn.style()),
                ])
            }
            Field::Save => Line::from(""),
        }
    }

    fn draw_form(&self, frame: &mut Frame, area: Rect) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(2),
                Constraint::Min(6),
                Constraint::Length(3),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let title = Paragraph::new(Line::from(Span::styled(
            "  cohort setup",
            Style::default().fg(theme::ACCENT).bold(),
        )));
        frame.render_widget(title, layout[0]);

        let lines: Vec<Line> = self
            .visible_fields()
            .into_iter()
            .filter(|f| *f != Field::Save)
            .map(|f| self.field_line(f))
            .collect();
        frame.render_widget(Paragraph::new(lines), layout[1]);

        let save_focused = self.focus == Field::Save;
        let save_label = if save_focused {
            format!(" {FOCUS_GLYPH} Save ")
        } else {
            "   Save ".to_string()
        };
        let save_block = Block::default()
            .borders(Borders::ALL)
            .border_style(Style::default().fg(if save_focused { theme::ACCENT } else { theme::MUTED }));
        let save = Paragraph::new(Line::from(Span::styled(
            save_label,
            if save_focused {
                Tone::Focus.style()
            } else {
                Tone::Normal.style()
            },
        )))
        .block(save_block);
        let save_area = Rect {
            x: layout[2].x + 2,
            y: layout[2].y,
            width: 12,
            height: 3,
        };
        frame.render_widget(save, save_area);

        let err_text = self.error.clone().unwrap_or_default();
        let err = Paragraph::new(format!("  {err_text}")).style(theme::error_slot_style());
        frame.render_widget(err, layout[3]);

        let hint = match self.focus {
            Field::Tier => "  ←/→ cycle tier   ↑/↓ move   enter activate   esc cancel",
            Field::Root => "  enter browse   ↑/↓ move   esc cancel",
            Field::Advanced => "  enter toggle   ↑/↓ move   esc cancel",
            Field::Env => "  ←/→ cycle env   ↑/↓ move   esc cancel",
            Field::Memory => "  ←/→ cycle preset   enter type custom   esc cancel",
            Field::Packs => "  enter edit packs   ↑/↓ move   esc cancel",
            Field::Save => "  enter save and exit   esc cancel",
        };
        frame.render_widget(Paragraph::new(hint).style(theme::hint_bar_style()), layout[4]);
    }

    fn draw_memory_edit(&self, frame: &mut Frame, area: Rect, buf: &str) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(2),
                Constraint::Length(3),
                Constraint::Min(0),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);
        let title = Paragraph::new(Line::from(Span::styled(
            "  cohort setup — memory budget",
            Style::default().fg(theme::ACCENT).bold(),
        )));
        frame.render_widget(title, layout[0]);

        let valid = buf.is_empty() || ResourceConfig::parse_memory_bytes(buf).is_some();
        let color = if valid { theme::FG } else { theme::BAD };
        let input = Paragraph::new(Line::from(vec![
            Span::styled("  amount: ", Tone::Muted.style()),
            Span::styled(buf.to_string(), Style::default().fg(color)),
            Span::styled("_", Style::default().fg(theme::ACCENT)),
        ]))
        .block(
            Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(input, layout[1]);

        let err = Paragraph::new(format!("  {}", self.error.clone().unwrap_or_default()))
            .style(theme::error_slot_style());
        frame.render_widget(err, layout[3]);

        let hint = "  type amount (e.g. 48GB)   enter confirm   esc cancel";
        frame.render_widget(Paragraph::new(hint).style(theme::hint_bar_style()), layout[4]);
    }

    fn draw_packs_edit(&self, frame: &mut Frame, area: Rect, cursor: usize, checked: &[bool]) {
        let packs = Pack::optional();
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(2),
                Constraint::Min(4),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);
        let title = Paragraph::new(Line::from(Span::styled(
            "  cohort setup — add-on packs",
            Style::default().fg(theme::ACCENT).bold(),
        )));
        frame.render_widget(title, layout[0]);

        let lines: Vec<Line> = packs
            .iter()
            .enumerate()
            .map(|(i, p)| {
                let g = if i == cursor { FOCUS_GLYPH } else { " " };
                let mark = if checked[i] { "[x]" } else { "[ ]" };
                let tone = if i == cursor { Tone::Focus } else { Tone::Normal };
                Line::from(vec![
                    Span::styled(format!(" {g} {mark} "), tone.style()),
                    Span::styled(format!("{:<14}", p.id), tone.style()),
                    Span::styled(format!("{:>7}  ", p.size_human), Tone::Warn.style()),
                    Span::styled(p.description.to_string(), Tone::Muted.style()),
                ])
            })
            .collect();
        frame.render_widget(Paragraph::new(lines), layout[1]);

        let err = Paragraph::new(format!("  {}", self.error.clone().unwrap_or_default()))
            .style(theme::error_slot_style());
        frame.render_widget(err, layout[2]);

        let hint = "  space toggle   a toggle all   enter done   esc cancel";
        frame.render_widget(Paragraph::new(hint).style(theme::hint_bar_style()), layout[3]);
    }
}

impl Screen for SetupScreen {
    fn title(&self) -> &str {
        "Setup"
    }

    fn scope(&self) -> ActionScope {
        match &self.mode {
            Mode::BrowseRoot(_) => ActionScope::FilePicker,
            _ => ActionScope::Setup,
        }
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        let mut map = KeyMap::new()
            .bind(KeyCode::Esc, none, Action::SetupCancel)
            .bind(KeyCode::Char('q'), none, Action::SetupCancel);
        match &self.mode {
            Mode::Form => {
                map = map
                    .bind(KeyCode::Up, none, Action::SetupPrev)
                    .bind(KeyCode::Char('k'), none, Action::SetupPrev)
                    .bind(KeyCode::Down, none, Action::SetupNext)
                    .bind(KeyCode::Char('j'), none, Action::SetupNext)
                    .bind(KeyCode::Tab, none, Action::SetupNext)
                    .bind(KeyCode::Enter, none, Action::SetupConfirm)
                    .bind(KeyCode::Char(' '), none, Action::SetupToggle);
            }
            Mode::EditPacks(_, _) => {
                map = map
                    .bind(KeyCode::Up, none, Action::SetupPrev)
                    .bind(KeyCode::Char('k'), none, Action::SetupPrev)
                    .bind(KeyCode::Down, none, Action::SetupNext)
                    .bind(KeyCode::Char('j'), none, Action::SetupNext)
                    .bind(KeyCode::Char(' '), none, Action::SetupToggle)
                    .bind(KeyCode::Char('a'), none, Action::SetupToggleAll)
                    .bind(KeyCode::Enter, none, Action::SetupConfirm);
            }
            Mode::EditMemory(_) => {
                map = map.bind(KeyCode::Enter, none, Action::SetupConfirm);
            }
            Mode::BrowseRoot(_) => {
                map = KeyMap::new()
                    .bind(KeyCode::Esc, none, Action::PickerCancel)
                    .bind(KeyCode::Char('q'), none, Action::PickerCancel)
                    .bind(KeyCode::Up, none, Action::PickerUp)
                    .bind(KeyCode::Char('k'), none, Action::PickerUp)
                    .bind(KeyCode::Down, none, Action::PickerDown)
                    .bind(KeyCode::Char('j'), none, Action::PickerDown)
                    .bind(KeyCode::Right, none, Action::PickerInto)
                    .bind(KeyCode::Enter, none, Action::PickerSelect)
                    .bind(KeyCode::Left, none, Action::PickerParent)
                    .bind(KeyCode::Backspace, none, Action::PickerParent)
                    .bind(KeyCode::Char(' '), none, Action::PickerSelect)
                    .bind(KeyCode::Char('c'), none, Action::PickerCreate)
                    .bind(KeyCode::Char('/'), none, Action::PickerTypePath)
                    .bind(KeyCode::Char('g'), none, Action::PickerTypePath);
            }
        }
        map
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, _log: &LogTail) {
        match &mut self.mode {
            Mode::Form => self.draw_form(frame, area),
            Mode::BrowseRoot(state) => file_picker::draw(frame, area, state),
            Mode::EditMemory(buf) => {
                let buf = buf.clone();
                self.draw_memory_edit(frame, area, &buf);
            }
            Mode::EditPacks(cursor, checked) => {
                let cursor = *cursor;
                let checked = checked.clone();
                self.draw_packs_edit(frame, area, cursor, &checked);
            }
        }
    }

    fn handle(&mut self, event: &AppEvent) -> Transition {
        let AppEvent::Key(k) = event else {
            return Transition::Stay;
        };
        match self.mode {
            Mode::Form => self.handle_form(k.code),
            Mode::BrowseRoot(_) => self.handle_browse(k.code),
            Mode::EditMemory(_) => self.handle_memory_edit(k.code),
            Mode::EditPacks(_, _) => self.handle_packs_edit(k.code),
        }
    }
}
