use std::path::PathBuf;

use crossterm::event::KeyCode;
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph, Widget};

use crate::tui::stages::types::{
    FieldData, FieldValue, FormField, FormSchema, FormValues, PathKind,
};
use crate::tui::theme::{self, Tone, FOCUS_GLYPH};

pub enum FormOutcome {
    Continue,
    Submit,
    Cancel,
    OpenAdvanced,
    RequestEdit(&'static str),
}

pub struct Form {
    schema: FormSchema,
    values: FormValues,
    advanced_open: bool,
    focus: usize,
    save_label: &'static str,
}

impl Form {
    pub fn new(schema: FormSchema, save_label: &'static str) -> Self {
        let mut values = FormValues::default();
        seed_defaults(&schema.fields, &mut values);
        seed_defaults(&schema.advanced, &mut values);
        Self {
            schema,
            values,
            advanced_open: false,
            focus: 0,
            save_label,
        }
    }

    pub fn values(&self) -> &FormValues {
        &self.values
    }

    pub fn set_path(&mut self, id: &'static str, path: Option<PathBuf>) {
        self.values
            .values
            .insert(id, FieldData::Path(FieldValue::Edited(path)));
    }

    pub fn set_text(&mut self, id: &'static str, text: String) {
        self.values
            .values
            .insert(id, FieldData::Text(FieldValue::Edited(text)));
    }

    pub fn set_multi(&mut self, id: &'static str, items: Vec<String>) {
        self.values
            .values
            .insert(id, FieldData::Multi(FieldValue::Edited(items)));
    }

    fn rows(&self) -> Vec<Row<'_>> {
        let mut rows: Vec<Row> = self
            .schema
            .fields
            .iter()
            .map(|f| Row::Field { field: f, advanced: false })
            .collect();
        if !self.schema.advanced.is_empty() {
            rows.push(Row::AdvancedExpander);
            if self.advanced_open {
                for f in &self.schema.advanced {
                    rows.push(Row::Field { field: f, advanced: true });
                }
            }
        }
        rows.push(Row::Save);
        rows
    }

    pub fn focused_field_id(&self) -> Option<&'static str> {
        match self.rows().get(self.focus)? {
            Row::Field { field, .. } => Some(field.id()),
            _ => None,
        }
    }

    pub fn handle(&mut self, code: KeyCode) -> FormOutcome {
        let row_count = self.rows().len();
        if row_count == 0 {
            return FormOutcome::Continue;
        }
        if self.focus >= row_count {
            self.focus = 0;
        }
        match code {
            KeyCode::Esc | KeyCode::Char('q') => FormOutcome::Cancel,
            KeyCode::Up | KeyCode::Char('k') => {
                self.focus = (self.focus + row_count - 1) % row_count;
                FormOutcome::Continue
            }
            KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                self.focus = (self.focus + 1) % row_count;
                FormOutcome::Continue
            }
            KeyCode::Left | KeyCode::Char('h') => {
                if let Some(field) = self.focused_field_owned() {
                    self.cycle(&field, false);
                }
                FormOutcome::Continue
            }
            KeyCode::Right | KeyCode::Char('l') => {
                if let Some(field) = self.focused_field_owned() {
                    self.cycle(&field, true);
                }
                FormOutcome::Continue
            }
            KeyCode::Char(' ') | KeyCode::Enter => self.activate(),
            _ => FormOutcome::Continue,
        }
    }

    fn focused_field_owned(&self) -> Option<FormField> {
        match self.rows().get(self.focus)? {
            Row::Field { field, .. } => Some((*field).clone()),
            _ => None,
        }
    }

    fn activate(&mut self) -> FormOutcome {
        let row_kind = match self.rows().get(self.focus) {
            Some(Row::Save) => ActivateKind::Save,
            Some(Row::AdvancedExpander) => ActivateKind::Advanced,
            Some(Row::Field { field, .. }) => ActivateKind::Field((*field).clone()),
            None => return FormOutcome::Continue,
        };
        match row_kind {
            ActivateKind::Save => FormOutcome::Submit,
            ActivateKind::Advanced => {
                self.advanced_open = !self.advanced_open;
                FormOutcome::OpenAdvanced
            }
            ActivateKind::Field(field) => match field {
                FormField::Choice { .. } | FormField::Toggle { .. } => {
                    self.cycle(&field, true);
                    FormOutcome::Continue
                }
                FormField::Path { id, .. }
                | FormField::Text { id, .. }
                | FormField::Number { id, .. }
                | FormField::MultiSelect { id, .. } => FormOutcome::RequestEdit(id),
            },
        }
    }

    fn cycle(&mut self, field: &FormField, forward: bool) {
        match field {
            FormField::Choice { id, options, .. } => {
                let cur = self.values.choice(id).map(str::to_string);
                let idx = options.iter().position(|o| Some(*o) == cur.as_deref());
                let next = match (idx, forward) {
                    (None, _) => 0,
                    (Some(i), true) => (i + 1) % options.len(),
                    (Some(i), false) => (i + options.len() - 1) % options.len(),
                };
                self.values.values.insert(
                    id,
                    FieldData::Choice(FieldValue::Edited(Some(options[next].to_string()))),
                );
            }
            FormField::Toggle { id, .. } => {
                let cur = self.values.toggle(id).unwrap_or(false);
                self.values
                    .values
                    .insert(id, FieldData::Toggle(FieldValue::Edited(!cur)));
            }
            _ => {}
        }
    }
}

#[derive(Clone, Copy)]
enum Row<'a> {
    Field { field: &'a FormField, advanced: bool },
    AdvancedExpander,
    Save,
}

enum ActivateKind {
    Save,
    Advanced,
    Field(FormField),
}

fn seed_defaults(fields: &[FormField], values: &mut FormValues) {
    for f in fields {
        match f {
            FormField::Path { id, default, .. } => {
                values
                    .values
                    .insert(id, FieldData::Path(FieldValue::Inferred(default.clone())));
            }
            FormField::Choice { id, default, .. } => {
                values.values.insert(
                    id,
                    FieldData::Choice(FieldValue::Inferred(default.map(str::to_string))),
                );
            }
            FormField::Toggle { id, default, .. } => {
                values
                    .values
                    .insert(id, FieldData::Toggle(FieldValue::Inferred(*default)));
            }
            FormField::Text { id, default, .. } => {
                values.values.insert(
                    id,
                    FieldData::Text(FieldValue::Inferred(default.clone().unwrap_or_default())),
                );
            }
            FormField::Number { id, default, .. } => {
                values.values.insert(
                    id,
                    FieldData::Number(FieldValue::Inferred(default.unwrap_or(0.0))),
                );
            }
            FormField::MultiSelect { id, default, .. } => {
                values.values.insert(
                    id,
                    FieldData::Multi(FieldValue::Inferred(
                        default.iter().map(|s| s.to_string()).collect(),
                    )),
                );
            }
        }
    }
}

impl Form {
    pub fn render(&self, area: Rect, buf: &mut Buffer) {
        let rows = self.rows();
        let mut lines: Vec<Line> = Vec::with_capacity(rows.len());
        for (i, row) in rows.iter().enumerate() {
            if matches!(row, Row::Save) {
                continue;
            }
            let focused = i == self.focus;
            lines.push(render_row(row, focused, &self.values));
        }
        let body_height = area.height.saturating_sub(3);
        let body_area = Rect { height: body_height, ..area };
        Paragraph::new(lines).render(body_area, buf);

        let save_focused = matches!(rows.get(self.focus), Some(Row::Save));
        let save_y = area.y + body_height;
        if save_y + 3 <= area.y + area.height {
            let save_area = Rect {
                x: area.x + 2,
                y: save_y,
                width: (self.save_label.len() as u16 + 6).min(area.width.saturating_sub(2)),
                height: 3,
            };
            let glyph = if save_focused { FOCUS_GLYPH } else { " " };
            let label = format!(" {glyph} {} ", self.save_label);
            let border_color = if save_focused { theme::ACCENT } else { theme::MUTED };
            let block = Block::default()
                .borders(Borders::ALL)
                .border_style(Style::default().fg(border_color));
            let tone = if save_focused { Tone::Focus } else { Tone::Normal };
            Paragraph::new(Line::from(Span::styled(label, tone.style())))
                .block(block)
                .render(save_area, buf);
        }
    }
}

fn render_row(row: &Row<'_>, focused: bool, values: &FormValues) -> Line<'static> {
    let label_tone = if focused { Tone::Focus } else { Tone::Normal };
    let glyph = if focused { FOCUS_GLYPH } else { " " };
    match row {
        Row::Save => Line::from(""),
        Row::AdvancedExpander => Line::from(vec![
            Span::styled(format!(" {glyph} "), label_tone.style()),
            Span::styled("[+] advanced", label_tone.style()),
        ]),
        Row::Field { field, advanced } => {
            let indent = if *advanced { "   " } else { "" };
            let label = field_label(field);
            let (value_text, edited) = value_display(field, values);
            let value_style = if edited {
                Style::default().fg(theme::WARN)
            } else {
                Style::default().fg(theme::MUTED).add_modifier(Modifier::DIM)
            };
            Line::from(vec![
                Span::styled(format!("{indent} {glyph} "), label_tone.style()),
                Span::styled(format!("{label:<10}"), label_tone.style()),
                Span::styled(value_text, value_style),
            ])
        }
    }
}

fn field_label(f: &FormField) -> &'static str {
    match f {
        FormField::Path { label, .. }
        | FormField::Choice { label, .. }
        | FormField::Toggle { label, .. }
        | FormField::Text { label, .. }
        | FormField::Number { label, .. }
        | FormField::MultiSelect { label, .. } => label,
    }
}

fn value_display(f: &FormField, values: &FormValues) -> (String, bool) {
    match f {
        FormField::Path { id, kind, .. } => {
            let edited = matches!(
                values.values.get(id),
                Some(FieldData::Path(FieldValue::Edited(_)))
            );
            let text = values
                .path(id)
                .map(|p| p.display().to_string())
                .unwrap_or_else(|| match kind {
                    PathKind::Dir => "<dir>".to_string(),
                    PathKind::File => "<file>".to_string(),
                    PathKind::Any => "<path>".to_string(),
                });
            (text, edited)
        }
        FormField::Choice { id, .. } => {
            let edited = matches!(
                values.values.get(id),
                Some(FieldData::Choice(FieldValue::Edited(_)))
            );
            let text = values.choice(id).unwrap_or("<unset>").to_string();
            (text, edited)
        }
        FormField::Toggle { id, .. } => {
            let edited = matches!(
                values.values.get(id),
                Some(FieldData::Toggle(FieldValue::Edited(_)))
            );
            let text = if values.toggle(id).unwrap_or(false) {
                "on".to_string()
            } else {
                "off".to_string()
            };
            (text, edited)
        }
        FormField::Text { id, .. } => {
            let edited = matches!(
                values.values.get(id),
                Some(FieldData::Text(FieldValue::Edited(_)))
            );
            let text = values.text(id).unwrap_or("").to_string();
            let text = if text.is_empty() { "<empty>".to_string() } else { text };
            (text, edited)
        }
        FormField::Number { id, .. } => {
            let edited = matches!(
                values.values.get(id),
                Some(FieldData::Number(FieldValue::Edited(_)))
            );
            let text = values
                .number(id)
                .map(|n| format!("{n}"))
                .unwrap_or_else(|| "0".to_string());
            (text, edited)
        }
        FormField::MultiSelect { id, .. } => {
            let edited = matches!(
                values.values.get(id),
                Some(FieldData::Multi(FieldValue::Edited(_)))
            );
            let text = values
                .multi(id)
                .filter(|v| !v.is_empty())
                .map(|v| v.join(", "))
                .unwrap_or_else(|| "none".to_string());
            (text, edited)
        }
    }
}
