use ratatui::layout::Rect;
use ratatui::Frame;

use super::action::{Action, ActionScope, KeyMap};
use super::event::AppEvent;
use super::stages::types::RunRequest;
use super::state::SessionState;

pub trait Screen {
    fn title(&self) -> &str;
    fn draw(&mut self, frame: &mut Frame, area: Rect);
    fn scope(&self) -> ActionScope;

    fn keys(&self) -> KeyMap {
        KeyMap::for_scope(self.scope())
    }

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
    Run(RunRequest),
}
