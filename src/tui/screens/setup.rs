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
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::stages::types::{FormField, FormSchema, PathKind};
use crate::tui::theme::{self, Tone, FOCUS_GLYPH};
use crate::tui::widgets::file_picker::{self, tab_complete, DirBrowserState};
use crate::tui::widgets::form::{Form, FormOutcome};

const MEMORY_PRESETS: &[&str] = &["auto", "8GB", "16GB", "32GB", "64GB", "128GB", "256GB"];

#[derive(Debug, Clone)]
pub struct SetupOutcome {
    pub tier: Tier,
    pub root: PathBuf,
    pub packs: Vec<String>,
    pub environment: Option<Environment>,
    pub memory_budget: Option<String>,
}

pub type OutcomeSink = Arc<Mutex<Option<SetupOutcome>>>;

enum Mode {
    Form,
    BrowseRoot(DirBrowserState),
    EditMemory(String),
    EditPacks(usize, Vec<bool>),
}

pub struct SetupScreen {
    mode: Mode,
    form: Form,
    sink: Option<OutcomeSink>,
    error: Option<String>,
}

fn build_schema(cfg: &Config) -> FormSchema {
    FormSchema {
        fields: vec![
            FormField::Choice {
                id: "tier",
                label: "tier",
                options: &["base", "full"],
                default: Some(match cfg.data.tier {
                    Tier::Base => "base",
                    Tier::Full => "full",
                }),
            },
            FormField::Path {
                id: "root",
                label: "data root",
                kind: PathKind::Dir,
                default: None,
            },
        ],
        advanced: vec![
            FormField::Choice {
                id: "env",
                label: "env",
                options: &["auto", "hpc", "workstation"],
                default: Some(match cfg.resources.environment {
                    None => "auto",
                    Some(Environment::Hpc) => "hpc",
                    Some(Environment::Workstation) => "workstation",
                }),
            },
            FormField::Choice {
                id: "memory",
                label: "memory",
                options: MEMORY_PRESETS,
                default: Some("auto"),
            },
            FormField::MultiSelect {
                id: "packs",
                label: "packs",
                default: &[],
            },
        ],
    }
}

impl SetupScreen {
    pub fn new() -> Self {
        let cfg = Config::load().unwrap_or_default();
        let root = if cfg.data.root_dir.is_empty() {
            Config::default_root_dir()
        } else {
            PathBuf::from(&cfg.data.root_dir)
        };
        let schema = build_schema(&cfg);
        let mut form = Form::new(schema, "Save");
        form.set_path("root", Some(root));
        if let Some(mem) = cfg.resources.memory_budget.as_deref() {
            if MEMORY_PRESETS.contains(&mem) {
                form.set_text("memory", mem.to_string());
            }
        }
        Self {
            mode: Mode::Form,
            form,
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
        if let Some(e) = env {
            let opt = match e {
                Environment::Hpc => "hpc",
                Environment::Workstation => "workstation",
            };
            s.form.set_text("env", opt.to_string());
        }
        if let Some(mb) = memory_budget {
            s.form.set_text("memory", mb);
        }
        s.sink = Some(sink);
        s
    }

    fn current_tier(&self) -> Tier {
        match self.form.values().choice("tier") {
            Some("full") => Tier::Full,
            _ => Tier::Base,
        }
    }

    fn current_root(&self) -> PathBuf {
        self.form
            .values()
            .path("root")
            .cloned()
            .unwrap_or_else(Config::default_root_dir)
    }

    fn current_env(&self) -> Option<Environment> {
        match self.form.values().choice("env") {
            Some("hpc") => Some(Environment::Hpc),
            Some("workstation") => Some(Environment::Workstation),
            _ => None,
        }
    }

    fn current_memory(&self) -> Option<String> {
        match self.form.values().choice("memory") {
            Some(v) if v != "auto" && !v.is_empty() => Some(v.to_string()),
            _ => None,
        }
    }

    fn current_packs(&self) -> Vec<String> {
        self.form
            .values()
            .multi("packs")
            .cloned()
            .unwrap_or_default()
    }

    fn save(&mut self) -> Transition {
        let root = self.current_root();
        if !root.is_dir() {
            self.error = Some(format!("data root not a directory: {}", root.display()));
            return Transition::Stay;
        }
        if let Some(sink) = self.sink.as_ref() {
            *sink.lock().unwrap() = Some(SetupOutcome {
                tier: self.current_tier(),
                root,
                packs: self.current_packs(),
                environment: self.current_env(),
                memory_budget: self.current_memory(),
            });
        }
        Transition::Pop
    }

    fn open_browse(&mut self) {
        let root = self.current_root();
        let start = if root.is_dir() {
            root
        } else {
            std::env::current_dir().unwrap_or_else(|_| Config::default_root_dir())
        };
        self.mode = Mode::BrowseRoot(DirBrowserState::new("Select FAVOR data root", &start));
    }

    fn open_edit_packs(&mut self) {
        let optional = Pack::optional();
        let current = self.current_packs();
        let checked: Vec<bool> = optional
            .iter()
            .map(|p| current.iter().any(|id| id == p.id))
            .collect();
        self.mode = Mode::EditPacks(0, checked);
    }

    fn open_edit_memory(&mut self) {
        let buf = self.current_memory().unwrap_or_default();
        self.mode = Mode::EditMemory(buf);
    }

    fn handle_form(&mut self, code: KeyCode) -> Transition {
        self.error = None;
        match self.form.handle(code) {
            FormOutcome::Continue | FormOutcome::OpenAdvanced => Transition::Stay,
            FormOutcome::Cancel => Transition::Pop,
            FormOutcome::Submit => self.save(),
            FormOutcome::RequestEdit(id) => {
                match id {
                    "root" => self.open_browse(),
                    "memory" => self.open_edit_memory(),
                    "packs" => self.open_edit_packs(),
                    _ => {}
                }
                Transition::Stay
            }
        }
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
                    self.form.set_path("root", Some(sel));
                    self.mode = Mode::Form;
                }
            }
            KeyCode::Left | KeyCode::Backspace => state.go_parent(),
            KeyCode::Enter | KeyCode::Char(' ') => {
                let dir = state.current_dir.clone();
                self.form.set_path("root", Some(dir));
                self.mode = Mode::Form;
            }
            KeyCode::Char('c') => {
                let _ = std::fs::create_dir_all(&state.current_dir);
                let dir = state.current_dir.clone();
                self.form.set_path("root", Some(dir));
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
                    self.form.set_text("memory", "auto".to_string());
                    self.mode = Mode::Form;
                } else if ResourceConfig::parse_memory_bytes(buf).is_some() {
                    self.form.set_text("memory", buf.clone());
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
                let packs: Vec<String> = Pack::optional()
                    .into_iter()
                    .zip(checked.iter())
                    .filter(|(_, &on)| on)
                    .map(|(p, _)| p.id.to_string())
                    .collect();
                self.form.set_multi("packs", packs);
                self.mode = Mode::Form;
            }
            KeyCode::Esc | KeyCode::Char('q') => self.mode = Mode::Form,
            _ => {}
        }
        Transition::Stay
    }

    fn draw_form(&self, frame: &mut Frame, area: Rect) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(2),
                Constraint::Min(6),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let title = Paragraph::new(Line::from(Span::styled(
            "  cohort setup",
            Style::default().fg(theme::ACCENT).bold(),
        )));
        frame.render_widget(title, layout[0]);

        let body = layout[1];
        self.form.render(body, frame.buffer_mut());

        let err_text = self.error.clone().unwrap_or_default();
        let err = Paragraph::new(format!("  {err_text}")).style(theme::error_slot_style());
        frame.render_widget(err, layout[2]);

        let hint = match self.form.focused_field_id() {
            Some("tier") => "  ←/→ cycle tier   ↑/↓ move   enter activate   esc cancel",
            Some("root") => "  enter browse   ↑/↓ move   esc cancel",
            Some("env") => "  ←/→ cycle env   ↑/↓ move   esc cancel",
            Some("memory") => "  ←/→ cycle preset   enter type custom   esc cancel",
            Some("packs") => "  enter edit packs   ↑/↓ move   esc cancel",
            _ => "  enter save and exit   esc cancel",
        };
        frame.render_widget(Paragraph::new(hint).style(theme::hint_bar_style()), layout[3]);
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

    fn draw(&mut self, frame: &mut Frame, area: Rect) {
        if area.width < 80 || area.height < 24 {
            let msg = Paragraph::new(Line::from(Span::styled(
                "terminal too small (need 80x24)",
                Style::default().fg(theme::WARN),
            )));
            frame.render_widget(msg, area);
            return;
        }
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
