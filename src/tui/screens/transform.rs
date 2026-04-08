use std::path::PathBuf;

use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, Paragraph};
use ratatui::Frame;

use crate::cli::GenomeBuild;
use crate::commands::{AnnotateConfig, IngestConfig};
use crate::config::{Config, Tier};
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::screen::{RunRequest, Screen, Transition};
use crate::tui::state::artifacts::{Artifact, ArtifactKind};
use crate::tui::theme;
use crate::tui::widgets::file_picker::{self, DirBrowserState};
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

pub struct IngestForm {
    inputs: Vec<PathBuf>,
    output: Option<PathBuf>,
    emit_sql: bool,
    build: Option<GenomeBuild>,
}

impl IngestForm {
    fn from_artifact(art: Option<&Artifact>) -> Self {
        let inputs = match art {
            Some(a) if matches!(a.kind, ArtifactKind::RawVcf) => vec![a.path.clone()],
            _ => Vec::new(),
        };
        Self {
            inputs,
            output: None,
            emit_sql: false,
            build: None,
        }
    }

    fn field_count(&self) -> usize {
        5
    }

    fn field_label(&self, idx: usize) -> &'static str {
        match idx {
            0 => "inputs",
            1 => "output",
            2 => "emit SQL",
            3 => "build",
            4 => "[ Run ]",
            _ => "",
        }
    }

    fn field_value(&self, idx: usize) -> String {
        match idx {
            0 => {
                if self.inputs.is_empty() {
                    "(none — press Enter to add)".into()
                } else if self.inputs.len() == 1 {
                    self.inputs[0].display().to_string()
                } else {
                    format!("{} files", self.inputs.len())
                }
            }
            1 => match &self.output {
                Some(p) => p.display().to_string(),
                None => "(derived from first input)".into(),
            },
            2 => if self.emit_sql { "yes" } else { "no" }.into(),
            3 => match self.build {
                Some(GenomeBuild::Hg38) => "hg38".into(),
                Some(GenomeBuild::Hg19) => "hg19".into(),
                None => "(auto-detect)".into(),
            },
            _ => String::new(),
        }
    }

    fn command_preview(&self) -> String {
        let mut parts = vec!["> cohort ingest".to_string()];
        for p in &self.inputs {
            parts.push(p.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default());
        }
        if let Some(o) = &self.output {
            parts.push(format!("-o {}", o.display()));
        }
        if self.emit_sql {
            parts.push("--emit-sql".into());
        }
        match self.build {
            Some(GenomeBuild::Hg38) => parts.push("--build hg38".into()),
            Some(GenomeBuild::Hg19) => parts.push("--build hg19".into()),
            None => {}
        }
        parts.join(" ")
    }

    fn to_config(&self) -> Result<IngestConfig, String> {
        if self.inputs.is_empty() {
            return Err("At least one input file is required.".into());
        }
        for p in &self.inputs {
            if !p.exists() {
                return Err(format!("Input not found: {}", p.display()));
            }
        }
        let output = self.output.clone().unwrap_or_else(|| default_ingest_output(&self.inputs[0]));
        Ok(IngestConfig {
            inputs: self.inputs.clone(),
            output,
            emit_sql: self.emit_sql,
            build_override: self.build.clone(),
        })
    }
}

fn default_ingest_output(first: &std::path::Path) -> PathBuf {
    let stem = first.file_stem().unwrap_or_default().to_string_lossy().into_owned();
    let stem = stem
        .strip_suffix(".vcf")
        .or_else(|| stem.strip_suffix(".tsv"))
        .or_else(|| stem.strip_suffix(".csv"))
        .unwrap_or(&stem);
    let stem = stem
        .split("_b0_")
        .next()
        .or_else(|| stem.split("_b0.").next())
        .unwrap_or(stem);
    first
        .parent()
        .unwrap_or(first)
        .join(format!("{stem}.ingested"))
}

pub struct AnnotateForm {
    input: PathBuf,
    output: Option<PathBuf>,
    tier: Tier,
    data_root: PathBuf,
}

impl AnnotateForm {
    fn from_artifact(art: &Artifact) -> Self {
        let cfg = Config::load().ok();
        let tier = cfg.as_ref().map(|c| c.data.tier).unwrap_or(Tier::Base);
        let data_root = cfg.map(|c| c.root_dir()).unwrap_or_else(PathBuf::new);
        Self {
            input: art.path.clone(),
            output: None,
            tier,
            data_root,
        }
    }

    fn field_count(&self) -> usize {
        5
    }

    fn field_label(&self, idx: usize) -> &'static str {
        match idx {
            0 => "input",
            1 => "output",
            2 => "tier",
            3 => "data root",
            4 => "[ Run ]",
            _ => "",
        }
    }

    fn field_value(&self, idx: usize) -> String {
        match idx {
            0 => self.input.display().to_string(),
            1 => match &self.output {
                Some(p) => p.display().to_string(),
                None => "(derived from input)".into(),
            },
            2 => self.tier.as_str().into(),
            3 => {
                if self.data_root.as_os_str().is_empty() {
                    "(not configured — run setup)".into()
                } else {
                    self.data_root.display().to_string()
                }
            }
            _ => String::new(),
        }
    }

    fn command_preview(&self) -> String {
        let mut parts = vec![
            "> cohort annotate".to_string(),
            self.input.file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_default(),
        ];
        if matches!(self.tier, Tier::Full) {
            parts.push("--full".into());
        }
        if let Some(o) = &self.output {
            parts.push(format!("-o {}", o.display()));
        }
        parts.join(" ")
    }

    fn to_config(&self) -> Result<AnnotateConfig, String> {
        if !self.input.exists() {
            return Err(format!("Input not found: {}", self.input.display()));
        }
        if !self.input.join("meta.json").exists() {
            return Err(format!(
                "Input '{}' is not an ingested set (missing meta.json).",
                self.input.display()
            ));
        }
        if self.data_root.as_os_str().is_empty() {
            return Err("Data root not configured. Run `s` from workspace to open setup.".into());
        }
        if !self.data_root.exists() {
            return Err(format!("Data root not found: {}", self.data_root.display()));
        }
        let output = self.output.clone().unwrap_or_else(|| default_annotate_output(&self.input));
        Ok(AnnotateConfig {
            input: self.input.clone(),
            output,
            tier: self.tier,
            data_root: self.data_root.clone(),
        })
    }
}

fn default_annotate_output(input: &std::path::Path) -> PathBuf {
    let name = input.file_name().unwrap_or_default().to_string_lossy().into_owned();
    let stem = name.strip_suffix(".ingested").unwrap_or(&name);
    input
        .parent()
        .unwrap_or(input)
        .join(format!("{stem}.annotated"))
}

enum FormState {
    Ingest(IngestForm),
    Annotate(AnnotateForm),
}

enum PickerTarget {
    IngestAddInput,
    IngestOutput,
    AnnotateOutput,
    AnnotateDataRoot,
}

pub struct TransformScreen {
    title: String,
    form: FormState,
    focus: usize,
    picker: Option<(DirBrowserState, PickerTarget)>,
    error: Option<String>,
}

impl TransformScreen {
    pub fn new_ingest(focused: Option<&Artifact>) -> Self {
        Self {
            title: "Transform: ingest".into(),
            form: FormState::Ingest(IngestForm::from_artifact(focused)),
            focus: 0,
            picker: None,
            error: None,
        }
    }

    pub fn new_annotate(art: &Artifact) -> Self {
        Self {
            title: "Transform: annotate".into(),
            form: FormState::Annotate(AnnotateForm::from_artifact(art)),
            focus: 0,
            picker: None,
            error: None,
        }
    }

    fn field_count(&self) -> usize {
        match &self.form {
            FormState::Ingest(f) => f.field_count(),
            FormState::Annotate(f) => f.field_count(),
        }
    }

    fn run_field_idx(&self) -> usize {
        self.field_count() - 1
    }

    fn cycle_focus(&mut self, delta: isize) {
        let n = self.field_count() as isize;
        let next = (self.focus as isize + delta).rem_euclid(n);
        self.focus = next as usize;
    }

    fn try_run(&mut self) -> Transition {
        let result: Result<RunRequest, String> = match &self.form {
            FormState::Ingest(f) => f.to_config().map(RunRequest::Ingest),
            FormState::Annotate(f) => f.to_config().map(RunRequest::Annotate),
        };
        match result {
            Ok(req) => Transition::Run(req),
            Err(msg) => {
                self.error = Some(msg);
                Transition::Stay
            }
        }
    }

    fn open_picker_for_focus(&mut self) {
        let (target, start, prompt, show_files) = match (&self.form, self.focus) {
            (FormState::Ingest(_), 0) => (
                PickerTarget::IngestAddInput,
                std::env::current_dir().unwrap_or_else(|_| PathBuf::from(".")),
                "Select input file",
                true,
            ),
            (FormState::Ingest(f), 1) => (
                PickerTarget::IngestOutput,
                f.output
                    .as_deref()
                    .and_then(|p| p.parent())
                    .map(PathBuf::from)
                    .or_else(|| std::env::current_dir().ok())
                    .unwrap_or_else(|| PathBuf::from(".")),
                "Select output directory",
                false,
            ),
            (FormState::Annotate(f), 1) => (
                PickerTarget::AnnotateOutput,
                f.output
                    .as_deref()
                    .and_then(|p| p.parent())
                    .map(PathBuf::from)
                    .or_else(|| std::env::current_dir().ok())
                    .unwrap_or_else(|| PathBuf::from(".")),
                "Select output directory",
                false,
            ),
            (FormState::Annotate(f), 3) => (
                PickerTarget::AnnotateDataRoot,
                if f.data_root.as_os_str().is_empty() || !f.data_root.exists() {
                    std::env::current_dir().unwrap_or_else(|_| PathBuf::from("/"))
                } else {
                    f.data_root.clone()
                },
                "Select data root directory",
                false,
            ),
            _ => return,
        };
        self.picker = Some((DirBrowserState::with_files(prompt, &start, show_files), target));
    }

    fn apply_picker_choice(&mut self, target: PickerTarget, chosen: PathBuf) {
        match (&mut self.form, target) {
            (FormState::Ingest(f), PickerTarget::IngestAddInput) => {
                f.inputs.push(chosen);
            }
            (FormState::Ingest(f), PickerTarget::IngestOutput) => {
                f.output = Some(chosen);
            }
            (FormState::Annotate(f), PickerTarget::AnnotateOutput) => {
                f.output = Some(chosen);
            }
            (FormState::Annotate(f), PickerTarget::AnnotateDataRoot) => {
                f.data_root = chosen;
            }
            _ => {}
        }
        self.error = None;
    }

    fn activate_field(&mut self) -> Transition {
        if self.focus == self.run_field_idx() {
            return self.try_run();
        }
        match (&mut self.form, self.focus) {
            (FormState::Ingest(f), 2) => {
                f.emit_sql = !f.emit_sql;
                Transition::Stay
            }
            (FormState::Ingest(f), 3) => {
                f.build = match f.build {
                    None => Some(GenomeBuild::Hg38),
                    Some(GenomeBuild::Hg38) => Some(GenomeBuild::Hg19),
                    Some(GenomeBuild::Hg19) => None,
                };
                Transition::Stay
            }
            (FormState::Annotate(f), 2) => {
                f.tier = match f.tier {
                    Tier::Base => Tier::Full,
                    Tier::Full => Tier::Base,
                };
                Transition::Stay
            }
            _ => {
                self.open_picker_for_focus();
                Transition::Stay
            }
        }
    }

    fn clear_input_at_focus(&mut self) {
        match (&mut self.form, self.focus) {
            (FormState::Ingest(f), 0) => {
                f.inputs.clear();
            }
            (FormState::Ingest(f), 1) => {
                f.output = None;
            }
            (FormState::Annotate(f), 1) => {
                f.output = None;
            }
            _ => {}
        }
    }

    fn command_preview(&self) -> String {
        match &self.form {
            FormState::Ingest(f) => f.command_preview(),
            FormState::Annotate(f) => f.command_preview(),
        }
    }
}

impl Screen for TransformScreen {
    fn title(&self) -> &str {
        &self.title
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail) {
        if let Some((picker, _)) = self.picker.as_mut() {
            file_picker::draw(frame, area, picker);
            return;
        }

        let v = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Min(6),
                Constraint::Length(1),
                Constraint::Length(1),
                Constraint::Length(6),
                Constraint::Length(1),
            ])
            .split(area);

        let items: Vec<ListItem> = (0..self.field_count())
            .map(|i| {
                let is_focus = i == self.focus;
                let (label, value) = match &self.form {
                    FormState::Ingest(f) => (f.field_label(i), f.field_value(i)),
                    FormState::Annotate(f) => (f.field_label(i), f.field_value(i)),
                };
                let label_style = if is_focus {
                    Style::default().fg(theme::ACCENT).bold()
                } else {
                    Style::default().fg(theme::MUTED)
                };
                let value_style = if is_focus {
                    Style::default().fg(theme::FG).bold()
                } else {
                    Style::default().fg(theme::FG)
                };
                let marker = if is_focus { " > " } else { "   " };
                ListItem::new(Line::from(vec![
                    Span::raw(marker),
                    Span::styled(format!("{label:<14}"), label_style),
                    Span::styled(value, value_style),
                ]))
            })
            .collect();

        let form_title = format!(" {} ", self.title);
        let list = List::new(items).block(
            Block::default()
                .borders(Borders::ALL)
                .title(form_title)
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(list, v[0]);

        let preview = Paragraph::new(Line::from(Span::styled(
            format!("  {}", self.command_preview()),
            Style::default().fg(theme::MUTED),
        )));
        frame.render_widget(preview, v[1]);

        let error_line = match &self.error {
            Some(msg) => Line::from(Span::styled(
                format!("  {msg}"),
                Style::default().fg(theme::BAD),
            )),
            None => Line::from(""),
        };
        frame.render_widget(Paragraph::new(error_line), v[2]);

        log.draw(frame, v[3], "Log");

        StatusBar {
            title: &self.title,
            keys: "tab/shift-tab move  enter activate  backspace clear  esc back",
        }
        .render(frame, v[4]);
    }

    fn scope(&self) -> ActionScope {
        if self.picker.is_some() {
            ActionScope::FilePicker
        } else {
            ActionScope::Transform
        }
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        if self.picker.is_some() {
            KeyMap::new()
                .bind(KeyCode::Esc, none, Action::PickerCancel)
                .bind(KeyCode::Up, none, Action::PickerUp)
                .bind(KeyCode::Char('k'), none, Action::PickerUp)
                .bind(KeyCode::Down, none, Action::PickerDown)
                .bind(KeyCode::Char('j'), none, Action::PickerDown)
                .bind(KeyCode::Left, none, Action::PickerParent)
                .bind(KeyCode::Backspace, none, Action::PickerParent)
                .bind(KeyCode::Right, none, Action::PickerInto)
                .bind(KeyCode::Enter, none, Action::PickerInto)
        } else {
            KeyMap::new()
                .bind(KeyCode::Esc, none, Action::TransformCancel)
                .bind(KeyCode::Tab, none, Action::TransformNextField)
                .bind(KeyCode::Down, none, Action::TransformNextField)
                .bind(KeyCode::BackTab, KeyModifiers::SHIFT, Action::TransformPrevField)
                .bind(KeyCode::Up, none, Action::TransformPrevField)
                .bind(KeyCode::Enter, none, Action::TransformActivate)
                .bind(KeyCode::Char(' '), none, Action::TransformToggleBool)
                .bind(KeyCode::Backspace, none, Action::TransformClearField)
        }
    }

    fn on_action(&mut self, action: Action) -> Transition {
        if let Some((picker, _)) = self.picker.as_mut() {
            match action {
                Action::PickerCancel => {
                    self.picker = None;
                }
                Action::PickerUp => picker.select_up(),
                Action::PickerDown => picker.select_down(),
                Action::PickerParent => picker.go_parent(),
                Action::PickerInto => {
                    if let Some(chosen) = picker.enter_selected() {
                        let (_, target) = self.picker.take().unwrap();
                        self.apply_picker_choice(target, chosen);
                    }
                }
                _ => {}
            }
            return Transition::Stay;
        }

        match action {
            Action::TransformCancel => Transition::Pop,
            Action::TransformNextField => {
                self.cycle_focus(1);
                Transition::Stay
            }
            Action::TransformPrevField => {
                self.cycle_focus(-1);
                Transition::Stay
            }
            Action::TransformActivate | Action::TransformToggleBool => self.activate_field(),
            Action::TransformClearField => {
                self.clear_input_at_focus();
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }
}
