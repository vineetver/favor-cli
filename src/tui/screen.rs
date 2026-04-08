use ratatui::layout::Rect;
use ratatui::Frame;

use crate::commands::{AnnotateConfig, IngestConfig};

use super::action::{Action, ActionScope, KeyMap};
use super::event::AppEvent;
use super::state::SessionState;
use super::widgets::log_tail::LogTail;

pub trait Screen {
    fn title(&self) -> &str;
    fn draw(&mut self, frame: &mut Frame, area: Rect, log: &LogTail);
    fn scope(&self) -> ActionScope;
    fn keys(&self) -> KeyMap;

    fn handle(&mut self, event: &AppEvent) -> Transition {
        match event {
            AppEvent::Key(k) => match self.keys().lookup(k.code, k.modifiers) {
                Some(action) => self.on_action(action),
                None => Transition::Stay,
            },
            other => self.on_other_event(other),
        }
    }

    fn on_action(&mut self, _action: Action) -> Transition {
        Transition::Stay
    }

    fn on_other_event(&mut self, _event: &AppEvent) -> Transition {
        Transition::Stay
    }

    fn on_focus(&mut self) {}

    fn contribute_session(&self, _state: &mut SessionState) {}

    fn restore_session(&mut self, _state: &SessionState) {}

    fn set_session_error(&mut self, _msg: String) {}
}

pub enum Transition {
    Stay,
    Quit,
    Push(Box<dyn Screen>),
    Pop,
    #[allow(dead_code)]
    Replace(Box<dyn Screen>),
    Run(RunRequest),
}

pub enum RunRequest {
    Ingest(IngestConfig),
    Annotate(AnnotateConfig),
}

impl RunRequest {
    pub fn description(&self) -> String {
        match self {
            RunRequest::Ingest(cfg) => {
                let first = cfg
                    .inputs
                    .first()
                    .map(|p| p.file_name().unwrap_or_default().to_string_lossy().into_owned())
                    .unwrap_or_default();
                if cfg.inputs.len() > 1 {
                    format!("Ingesting {first} (+{} more)", cfg.inputs.len() - 1)
                } else {
                    format!("Ingesting {first}")
                }
            }
            RunRequest::Annotate(cfg) => {
                let name = cfg
                    .input
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy();
                format!("Annotating {name} ({} tier)", cfg.tier.as_str())
            }
        }
    }
}
