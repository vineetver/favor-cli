use std::path::PathBuf;

use crossterm::event::KeyCode;
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::widgets::Widget;
use ratatui::Frame;

use crate::tui::action::{Action, ActionScope};
use crate::tui::screen::{Screen, Transition};
use crate::tui::shell::{ErrorMessage, ScreenChrome, Shell};
use crate::tui::stages::types::{FormField, PathKind, SessionCtx};
use crate::tui::stages::Stage;
use crate::tui::widgets::file_picker::{self, DirBrowserState};
use crate::tui::widgets::form::{Form, FormOutcome};

pub struct StageView {
    stage: &'static dyn Stage,
    title: String,
    form: Form,
    picker: Option<(DirBrowserState, &'static str)>,
    error: Option<ErrorMessage>,
}

impl StageView {
    pub fn new(stage: &'static dyn Stage, ctx: SessionCtx<'_>) -> Self {
        let schema = stage.form_schema(&ctx);
        let form = Form::new(schema, "Run");
        let title = format!("Stage: {}", stage.label());
        Self {
            stage,
            title,
            form,
            picker: None,
            error: None,
        }
    }

    fn open_picker(&mut self, id: &'static str) {
        let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));
        let kind = match self.form.field(id) {
            Some(FormField::Path { kind, .. }) => *kind,
            _ => return,
        };
        let start = self
            .form
            .values()
            .path(id)
            .and_then(|p| {
                if p.is_dir() {
                    Some(p.clone())
                } else {
                    p.parent().map(PathBuf::from)
                }
            })
            .unwrap_or(cwd);
        let show_files = !matches!(kind, PathKind::Dir);
        self.picker = Some((
            DirBrowserState::with_files("Select path", &start, show_files),
            id,
        ));
    }

    fn apply_picker(&mut self, id: &'static str, chosen: PathBuf) {
        self.form.set_path(id, Some(chosen));
        self.error = None;
    }

    fn try_run(&mut self) -> Transition {
        match self.stage.build_command(self.form.values()) {
            Ok(req) => Transition::Run(req),
            Err(err) => {
                self.error = Some(ErrorMessage { text: form_error_text(&err) });
                Transition::Stay
            }
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

    fn clear_focused(&mut self) {
        let Some(id) = self.form.focused_field_id() else {
            return;
        };
        if matches!(self.form.field(id), Some(FormField::Path { .. })) {
            self.form.set_path(id, None);
        }
    }
}

fn form_error_text(err: &crate::tui::stages::types::FormError) -> String {
    use crate::tui::stages::types::FormError;
    match err {
        FormError::Missing(field) => format!("missing field: {field}"),
        FormError::Invalid { field, reason } => format!("{field}: {reason}"),
    }
}

impl Screen for StageView {
    fn title(&self) -> &str {
        &self.title
    }

    fn draw(&mut self, frame: &mut Frame, area: Rect) {
        if let Some((picker, _)) = self.picker.as_mut() {
            file_picker::draw(frame, area, picker);
            return;
        }
        let chrome = ScreenChrome {
            title: &self.title,
            status: None,
            error: self.error.as_ref(),
            scope: ActionScope::Transform,
            graph: None,
        };
        let form = &self.form;
        let body = |inner: Rect, buf: &mut Buffer| {
            form.render(inner, buf);
        };
        Shell::new(chrome, body).render(area, frame.buffer_mut());
    }

    fn scope(&self) -> ActionScope {
        if self.picker.is_some() {
            ActionScope::FilePicker
        } else {
            ActionScope::Transform
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
                        let (_, id) = self.picker.take().unwrap();
                        self.apply_picker(id, chosen);
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
}
