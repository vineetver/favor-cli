use std::path::PathBuf;

use crossterm::event::KeyCode;
use ratatui::buffer::Buffer;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph, Widget};
use ratatui::Frame;

use crate::tui::action::{Action, ActionScope};
use crate::tui::screen::{Screen, Transition};
use crate::tui::shell::{ErrorMessage, ScreenChrome, Shell};
use crate::tui::stages::types::{FormField, PathKind, SessionCtx};
use crate::tui::stages::Stage;
use crate::tui::theme::{self, Tone, FOCUS_GLYPH};
use crate::tui::widgets::file_picker::{self, DirBrowserState};
use crate::tui::widgets::form::{Form, FormOutcome};

struct MultiPickerState {
    field_id: &'static str,
    label: &'static str,
    options: &'static [&'static str],
    cursor: usize,
    checked: Vec<bool>,
}

enum Editor {
    None,
    Path(DirBrowserState, &'static str),
    Multi(MultiPickerState),
}

pub struct StageView {
    stage: &'static dyn Stage,
    title: String,
    form: Form,
    editor: Editor,
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
            editor: Editor::None,
            error: None,
        }
    }

    fn open_path_picker(&mut self, id: &'static str) {
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
        self.editor = Editor::Path(
            DirBrowserState::with_files("Select path", &start, show_files),
            id,
        );
    }

    fn open_multi_picker(&mut self, id: &'static str) {
        let (label, options) = match self.form.field(id) {
            Some(FormField::MultiSelect { label, options, .. }) => (*label, *options),
            _ => return,
        };
        let current = self.form.values().multi(id).cloned().unwrap_or_default();
        let checked: Vec<bool> = options.iter().map(|o| current.iter().any(|c| c == o)).collect();
        self.editor = Editor::Multi(MultiPickerState {
            field_id: id,
            label,
            options,
            cursor: 0,
            checked,
        });
    }

    fn apply_path(&mut self, id: &'static str, chosen: PathBuf) {
        self.form.set_path(id, Some(chosen));
        self.error = None;
    }

    fn commit_multi(&mut self) {
        if let Editor::Multi(state) = &self.editor {
            let id = state.field_id;
            let picks: Vec<String> = state
                .options
                .iter()
                .zip(state.checked.iter())
                .filter(|(_, &on)| on)
                .map(|(o, _)| (*o).to_string())
                .collect();
            self.form.set_multi(id, picks);
        }
        self.editor = Editor::None;
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
                match self.form.field(id) {
                    Some(FormField::Path { .. }) => self.open_path_picker(id),
                    Some(FormField::MultiSelect { .. }) => self.open_multi_picker(id),
                    _ => {}
                }
                Transition::Stay
            }
        }
    }

    fn draw_multi(frame: &mut Frame, area: Rect, state: &MultiPickerState) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(2),
                Constraint::Min(4),
                Constraint::Length(1),
            ])
            .split(area);
        let title = Paragraph::new(Line::from(Span::styled(
            format!("  {}", state.label),
            Style::default().fg(theme::ACCENT).bold(),
        )));
        frame.render_widget(title, layout[0]);

        let lines: Vec<Line> = state
            .options
            .iter()
            .enumerate()
            .map(|(i, opt)| {
                let g = if i == state.cursor { FOCUS_GLYPH } else { " " };
                let mark = if state.checked[i] { "[x]" } else { "[ ]" };
                let tone = if i == state.cursor { Tone::Focus } else { Tone::Normal };
                Line::from(vec![
                    Span::styled(format!(" {g} {mark} "), tone.style()),
                    Span::styled((*opt).to_string(), tone.style()),
                ])
            })
            .collect();
        let body = Paragraph::new(lines).block(
            Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(body, layout[1]);

        let hint = "  space toggle   a toggle all   enter done   esc cancel";
        frame.render_widget(Paragraph::new(hint).style(theme::hint_bar_style()), layout[2]);
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
        match &mut self.editor {
            Editor::Path(picker, _) => {
                file_picker::draw(frame, area, picker);
                return;
            }
            Editor::Multi(state) => {
                Self::draw_multi(frame, area, state);
                return;
            }
            Editor::None => {}
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
        match &self.editor {
            Editor::Path(_, _) => ActionScope::FilePicker,
            _ => ActionScope::Transform,
        }
    }

    fn on_action(&mut self, action: Action) -> Transition {
        if let Editor::Path(picker, _) = &mut self.editor {
            match action {
                Action::PickerCancel => {
                    self.editor = Editor::None;
                }
                Action::PickerUp => picker.select_up(),
                Action::PickerDown => picker.select_down(),
                Action::PickerParent => picker.go_parent(),
                Action::PickerInto => {
                    if let Some(chosen) = picker.enter_selected() {
                        if let Editor::Path(_, id) = std::mem::replace(&mut self.editor, Editor::None) {
                            self.apply_path(id, chosen);
                        }
                    }
                }
                _ => {}
            }
            return Transition::Stay;
        }
        if let Editor::Multi(state) = &mut self.editor {
            match action {
                Action::TransformCancel => self.editor = Editor::None,
                Action::TransformPrevField => {
                    if state.cursor > 0 {
                        state.cursor -= 1;
                    }
                }
                Action::TransformNextField => {
                    if state.cursor + 1 < state.options.len() {
                        state.cursor += 1;
                    }
                }
                Action::TransformToggleBool => {
                    if let Some(c) = state.checked.get_mut(state.cursor) {
                        *c = !*c;
                    }
                }
                Action::TransformActivate => self.commit_multi(),
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
