use std::path::PathBuf;

use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::cli::GenomeBuild;
use crate::commands::{AnnotateConfig, IngestConfig};
use crate::config::{Config, Tier};
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::screen::{RunRequest, Screen, Transition};
use crate::tui::stages::types::{FormField, FormSchema, PathKind};
use crate::tui::state::artifacts::{Artifact, ArtifactKind};
use crate::tui::state::{BuildTag, SessionState, TransformSnapshot};
use crate::tui::theme;
use crate::tui::widgets::file_picker::{self, DirBrowserState};
use crate::tui::widgets::form::{Form, FormOutcome};
use crate::tui::widgets::log_tail::LogTail;
use crate::tui::widgets::status_bar::StatusBar;

enum Kind {
    Ingest,
    Annotate,
}

enum PickerTarget {
    IngestAddInput,
    IngestOutput,
    AnnotateOutput,
    AnnotateDataRoot,
}

pub struct TransformScreen {
    title: String,
    kind: Kind,
    form: Form,
    inputs: Vec<PathBuf>,
    annotate_input: PathBuf,
    picker: Option<(DirBrowserState, PickerTarget)>,
    error: Option<String>,
}

fn ingest_schema() -> FormSchema {
    FormSchema {
        fields: vec![
            FormField::Path {
                id: "inputs",
                label: "inputs",
                kind: PathKind::File,
                default: None,
            },
            FormField::Path {
                id: "output",
                label: "output",
                kind: PathKind::Dir,
                default: None,
            },
        ],
        advanced: vec![
            FormField::Toggle {
                id: "emit_sql",
                label: "emit SQL",
                default: false,
            },
            FormField::Choice {
                id: "build",
                label: "build",
                options: &["auto", "hg38", "hg19"],
                default: Some("auto"),
            },
        ],
    }
}

fn annotate_schema(default_tier: Tier) -> FormSchema {
    FormSchema {
        fields: vec![
            FormField::Path {
                id: "input",
                label: "input",
                kind: PathKind::Dir,
                default: None,
            },
            FormField::Path {
                id: "output",
                label: "output",
                kind: PathKind::Dir,
                default: None,
            },
            FormField::Choice {
                id: "tier",
                label: "tier",
                options: &["base", "full"],
                default: Some(default_tier.as_str()),
            },
            FormField::Path {
                id: "data_root",
                label: "data root",
                kind: PathKind::Dir,
                default: None,
            },
        ],
        advanced: vec![],
    }
}

impl TransformScreen {
    pub fn new_ingest(focused: Option<&Artifact>) -> Self {
        let inputs: Vec<PathBuf> = match focused {
            Some(a) if matches!(a.kind, ArtifactKind::RawVcf) => vec![a.path.clone()],
            _ => Vec::new(),
        };
        let mut form = Form::new(ingest_schema(), "Run");
        if let Some(first) = inputs.first() {
            form.set_path("inputs", Some(first.clone()));
        }
        Self {
            title: "Transform: ingest".into(),
            kind: Kind::Ingest,
            form,
            inputs,
            annotate_input: PathBuf::new(),
            picker: None,
            error: None,
        }
    }

    pub fn new_annotate(art: &Artifact) -> Self {
        let cfg = Config::load().ok();
        let tier = cfg.as_ref().map(|c| c.data.tier).unwrap_or(Tier::Base);
        let data_root = cfg.map(|c| c.root_dir()).unwrap_or_default();
        let mut form = Form::new(annotate_schema(tier), "Run");
        form.set_path("input", Some(art.path.clone()));
        if !data_root.as_os_str().is_empty() {
            form.set_path("data_root", Some(data_root));
        }
        Self {
            title: "Transform: annotate".into(),
            kind: Kind::Annotate,
            form,
            inputs: Vec::new(),
            annotate_input: art.path.clone(),
            picker: None,
            error: None,
        }
    }

    fn current_build(&self) -> Option<GenomeBuild> {
        match self.form.values().choice("build") {
            Some("hg38") => Some(GenomeBuild::Hg38),
            Some("hg19") => Some(GenomeBuild::Hg19),
            _ => None,
        }
    }

    fn current_tier(&self) -> Tier {
        match self.form.values().choice("tier") {
            Some("full") => Tier::Full,
            _ => Tier::Base,
        }
    }

    fn ingest_config(&self) -> Result<IngestConfig, String> {
        if self.inputs.is_empty() {
            return Err("At least one input file is required.".into());
        }
        for p in &self.inputs {
            if !p.exists() {
                return Err(format!("Input not found: {}", p.display()));
            }
        }
        let output = self
            .form
            .values()
            .path("output")
            .cloned()
            .unwrap_or_else(|| default_ingest_output(&self.inputs[0]));
        Ok(IngestConfig {
            inputs: self.inputs.clone(),
            output,
            emit_sql: self.form.values().toggle("emit_sql").unwrap_or(false),
            build_override: self.current_build(),
        })
    }

    fn annotate_config(&self) -> Result<AnnotateConfig, String> {
        let input = self.annotate_input.clone();
        if !input.exists() {
            return Err(format!("Input not found: {}", input.display()));
        }
        if !input.join("meta.json").exists() {
            return Err(format!(
                "Input '{}' is not an ingested set (missing meta.json).",
                input.display()
            ));
        }
        let data_root = self
            .form
            .values()
            .path("data_root")
            .cloned()
            .unwrap_or_default();
        if data_root.as_os_str().is_empty() {
            return Err("Data root not configured. Run `s` from workspace to open setup.".into());
        }
        if !data_root.exists() {
            return Err(format!("Data root not found: {}", data_root.display()));
        }
        let output = self
            .form
            .values()
            .path("output")
            .cloned()
            .unwrap_or_else(|| default_annotate_output(&input));
        Ok(AnnotateConfig {
            input,
            output,
            tier: self.current_tier(),
            data_root,
        })
    }

    fn try_run(&mut self) -> Transition {
        let result = match self.kind {
            Kind::Ingest => self.ingest_config().map(RunRequest::Ingest),
            Kind::Annotate => self.annotate_config().map(RunRequest::Annotate),
        };
        match result {
            Ok(req) => Transition::Run(req),
            Err(msg) => {
                self.error = Some(msg);
                Transition::Stay
            }
        }
    }

    fn open_picker(&mut self, id: &'static str) {
        let cwd = || std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        let (target, start, prompt, show_files) = match (&self.kind, id) {
            (Kind::Ingest, "inputs") => {
                (PickerTarget::IngestAddInput, cwd(), "Select input file", true)
            }
            (Kind::Ingest, "output") => {
                let start = self
                    .form
                    .values()
                    .path("output")
                    .and_then(|p| p.parent().map(PathBuf::from))
                    .unwrap_or_else(cwd);
                (PickerTarget::IngestOutput, start, "Select output directory", false)
            }
            (Kind::Annotate, "output") => {
                let start = self
                    .form
                    .values()
                    .path("output")
                    .and_then(|p| p.parent().map(PathBuf::from))
                    .unwrap_or_else(cwd);
                (PickerTarget::AnnotateOutput, start, "Select output directory", false)
            }
            (Kind::Annotate, "data_root") => {
                let cur = self.form.values().path("data_root").cloned();
                let start = match cur {
                    Some(p) if p.exists() => p,
                    _ => cwd(),
                };
                (PickerTarget::AnnotateDataRoot, start, "Select data root directory", false)
            }
            _ => return,
        };
        self.picker = Some((DirBrowserState::with_files(prompt, &start, show_files), target));
    }

    fn apply_picker(&mut self, target: PickerTarget, chosen: PathBuf) {
        match target {
            PickerTarget::IngestAddInput => {
                self.inputs.push(chosen.clone());
                self.form.set_path("inputs", Some(chosen));
            }
            PickerTarget::IngestOutput => self.form.set_path("output", Some(chosen)),
            PickerTarget::AnnotateOutput => self.form.set_path("output", Some(chosen)),
            PickerTarget::AnnotateDataRoot => self.form.set_path("data_root", Some(chosen)),
        }
        self.error = None;
    }

    fn clear_focused(&mut self) {
        let id = match self.form.focused_field_id() {
            Some(i) => i,
            None => return,
        };
        match (&self.kind, id) {
            (Kind::Ingest, "inputs") => {
                self.inputs.clear();
                self.form.set_path("inputs", None);
            }
            (_, "output") => self.form.set_path("output", None),
            _ => {}
        }
    }

    fn ingest_preview(&self) -> String {
        let mut parts = vec!["> cohort ingest".to_string()];
        for p in &self.inputs {
            parts.push(
                p.file_name()
                    .map(|n| n.to_string_lossy().into_owned())
                    .unwrap_or_default(),
            );
        }
        if let Some(o) = self.form.values().path("output") {
            parts.push(format!("-o {}", o.display()));
        }
        if self.form.values().toggle("emit_sql").unwrap_or(false) {
            parts.push("--emit-sql".into());
        }
        match self.current_build() {
            Some(GenomeBuild::Hg38) => parts.push("--build hg38".into()),
            Some(GenomeBuild::Hg19) => parts.push("--build hg19".into()),
            None => {}
        }
        parts.join(" ")
    }

    fn annotate_preview(&self) -> String {
        let mut parts = vec![
            "> cohort annotate".to_string(),
            self.annotate_input
                .file_name()
                .map(|n| n.to_string_lossy().into_owned())
                .unwrap_or_default(),
        ];
        if matches!(self.current_tier(), Tier::Full) {
            parts.push("--full".into());
        }
        if let Some(o) = self.form.values().path("output") {
            parts.push(format!("-o {}", o.display()));
        }
        parts.join(" ")
    }

    fn command_preview(&self) -> String {
        match self.kind {
            Kind::Ingest => self.ingest_preview(),
            Kind::Annotate => self.annotate_preview(),
        }
    }

    fn drive_form(&mut self, code: KeyCode) -> Transition {
        self.error = None;
        match self.form.handle(code) {
            FormOutcome::Continue | FormOutcome::OpenAdvanced => Transition::Stay,
            FormOutcome::Cancel => Transition::Pop,
            FormOutcome::Submit => self.try_run(),
            FormOutcome::RequestEdit(id) => {
                self.open_picker(id);
                Transition::Stay
            }
        }
    }
}

fn default_ingest_output(first: &std::path::Path) -> PathBuf {
    let stem = first
        .file_stem()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
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

fn default_annotate_output(input: &std::path::Path) -> PathBuf {
    let name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    let stem = name.strip_suffix(".ingested").unwrap_or(&name);
    input
        .parent()
        .unwrap_or(input)
        .join(format!("{stem}.annotated"))
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
                Constraint::Length(2),
                Constraint::Min(6),
                Constraint::Length(1),
                Constraint::Length(4),
                Constraint::Length(1),
                Constraint::Length(1),
            ])
            .split(area);

        let title = Paragraph::new(Line::from(Span::styled(
            format!("  {}", self.title),
            Style::default().fg(theme::ACCENT).bold(),
        )));
        frame.render_widget(title, v[0]);

        self.form.render(v[1], frame.buffer_mut());

        let preview = Paragraph::new(Line::from(Span::styled(
            format!("  {}", self.command_preview()),
            Style::default().fg(theme::MUTED),
        )));
        frame.render_widget(preview, v[2]);

        log.draw(frame, v[3], "Log");

        let error_line = match &self.error {
            Some(msg) => Line::from(Span::styled(
                format!("  {msg}"),
                Style::default().fg(theme::BAD),
            )),
            None => Line::from(""),
        };
        frame.render_widget(Paragraph::new(error_line), v[4]);

        StatusBar {
            title: &self.title,
            keys: "tab/shift-tab move  enter activate  backspace clear  esc back",
        }
        .render(frame, v[5]);
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
                        self.apply_picker(target, chosen);
                    }
                }
                _ => {}
            }
            return Transition::Stay;
        }

        match action {
            Action::TransformCancel => self.drive_form(KeyCode::Esc),
            Action::TransformNextField => self.drive_form(KeyCode::Down),
            Action::TransformPrevField => self.drive_form(KeyCode::Up),
            Action::TransformActivate => self.drive_form(KeyCode::Enter),
            Action::TransformToggleBool => self.drive_form(KeyCode::Char(' ')),
            Action::TransformClearField => {
                self.clear_focused();
                Transition::Stay
            }
            _ => Transition::Stay,
        }
    }

    fn contribute_session(&self, state: &mut SessionState) {
        state.transform = Some(match self.kind {
            Kind::Ingest => TransformSnapshot::Ingest {
                inputs: self.inputs.clone(),
                output: self.form.values().path("output").cloned(),
                emit_sql: self.form.values().toggle("emit_sql").unwrap_or(false),
                build: self.current_build().as_ref().map(BuildTag::from_build),
            },
            Kind::Annotate => TransformSnapshot::Annotate {
                input: self.annotate_input.clone(),
                output: self.form.values().path("output").cloned(),
                tier: self.current_tier(),
                data_root: self
                    .form
                    .values()
                    .path("data_root")
                    .cloned()
                    .unwrap_or_default(),
            },
        });
    }
}
