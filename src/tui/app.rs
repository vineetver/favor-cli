use std::io::{self, Stdout};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::{Arc, Mutex};
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind, KeyModifiers};
use crossterm::terminal::SetTitle;
use crossterm::ExecutableCommand;
use indicatif::ProgressBar;
use ratatui::backend::CrosstermBackend;
use ratatui::Terminal;

use crate::commands;
use crate::error::CohortError;
use crate::output::Output;

use super::action::{Action, KeyMap};
use super::event::AppEvent;
use super::output::{BarRegistry, LogLine, ProgressSnapshot, TuiOutput};
use super::screens::{help, stage_view, variant, workspace};
use super::stages::types::{RunRequest, SetupConfig};
use super::state::app::{AppState, Modal, Outcome, View};
use super::state::{SessionId, SessionState};
use super::widgets::palette::{Palette, PaletteOutcome};
use super::widgets::run_overlay::{dim_area, RunOverlay};
use crate::staar::pipeline::StaarPipeline;

struct BarEntry {
    arc: Arc<ProgressBar>,
}

pub struct App {
    state: AppState,
    log_rx: Receiver<LogLine>,
    shared_bars: BarRegistry,
    bars: Vec<BarEntry>,
    tui_out: Arc<TuiOutput>,
    cmd_tx: Sender<Result<(), CohortError>>,
    cmd_rx: Receiver<Result<(), CohortError>>,
    global_keys: KeyMap,
}

impl App {
    pub fn new(
        mut state: AppState,
        log_rx: Receiver<LogLine>,
        shared_bars: BarRegistry,
        tui_out: Arc<TuiOutput>,
    ) -> Self {
        let probe = SessionState {
            cwd: state.workspace.cwd.clone(),
            ..Default::default()
        };
        if let Some(store) = state.session_store.as_ref() {
            if !probe.cwd.as_os_str().is_empty() {
                let id = SessionId::for_cwd(&probe.cwd);
                match store.load(&id) {
                    Ok(Some(saved)) => state.restore_session(&saved),
                    Ok(None) => {}
                    Err(e) => {
                        state.error = Some(super::shell::ErrorMessage {
                            text: format!("session load failed: {e}"),
                        });
                    }
                }
            }
        }
        let (cmd_tx, cmd_rx) = mpsc::channel();
        let none = KeyModifiers::NONE;
        let global_keys = KeyMap::new()
            .bind(KeyCode::Char('c'), KeyModifiers::CONTROL, Action::Quit)
            .bind(KeyCode::Char('p'), KeyModifiers::CONTROL, Action::OpenPalette)
            .bind(KeyCode::Char('?'), none, Action::OpenHelp);
        Self {
            state,
            log_rx,
            shared_bars,
            bars: Vec::new(),
            tui_out,
            cmd_tx,
            cmd_rx,
            global_keys,
        }
    }

    pub fn with_setup_sink(mut self, sink: Arc<Mutex<Option<SetupConfig>>>) -> Self {
        self.state.setup_sink = Some(sink);
        self
    }

    pub fn run(
        mut self,
        terminal: &mut Terminal<CrosstermBackend<Stdout>>,
    ) -> Result<(), CohortError> {
        let tick = Duration::from_millis(50);
        let mut last_title = String::new();
        loop {
            let current_title = title_for(&self.state).to_string();
            if current_title != last_title {
                let _ = io::stdout().execute(SetTitle(format!("cohort — {current_title}")));
                last_title = current_title;
            }
            terminal
                .draw(|f| {
                    let area = f.area();
                    render(&mut self.state, f, area);
                    if let Some(Modal::Run(_)) = &self.state.modal {
                        let buf = f.buffer_mut();
                        dim_area(area, buf);
                    }
                    render_modal(&mut self.state, f, area);
                })
                .map_err(|e| CohortError::Internal(e.into()))?;

            let mut pending: Vec<AppEvent> = Vec::new();
            if event::poll(tick).map_err(|e| CohortError::Internal(e.into()))? {
                if let Event::Key(k) = event::read().map_err(|e| CohortError::Internal(e.into()))? {
                    if k.kind == KeyEventKind::Press {
                        if k.code == KeyCode::Char('c')
                            && k.modifiers.contains(KeyModifiers::CONTROL)
                        {
                            self.state.save_session();
                            return Ok(());
                        }
                        if let Some(Modal::Run(o)) = &mut self.state.modal {
                            if o.handle(k).is_some() {
                                self.state.modal = None;
                            }
                            continue;
                        }
                        if let Some(Modal::Palette(p)) = &mut self.state.modal {
                            match p.handle_key(k) {
                                PaletteOutcome::Stay => {}
                                PaletteOutcome::Close => {
                                    self.state.modal = None;
                                }
                                PaletteOutcome::Run(action) => {
                                    self.state.modal = None;
                                    let outcome = self.dispatch_palette_action(action);
                                    if self.apply_outcome(outcome) {
                                        return Ok(());
                                    }
                                }
                            }
                            continue;
                        }
                        if let Some(Modal::Help(h)) = &mut self.state.modal {
                            match k.code {
                                KeyCode::Esc => self.state.modal = None,
                                KeyCode::Tab => h.cycle(),
                                _ => {}
                            }
                            continue;
                        }
                        pending.push(AppEvent::Key(k));
                    }
                }
            } else {
                pending.push(AppEvent::Tick);
            }
            pending.extend(self.sample_bars());

            let mut outcomes: Vec<Outcome> = Vec::new();
            for ev in pending {
                self.apply_to_log(&ev);
                outcomes.push(self.dispatch_event(ev));
            }
            for outcome in outcomes {
                if self.apply_outcome(outcome) {
                    return Ok(());
                }
            }
        }
    }

    fn dispatch_event(&mut self, ev: AppEvent) -> Outcome {
        match ev {
            AppEvent::Tick => {
                if matches!(self.state.view, View::Workspace) {
                    workspace::on_tick(&mut self.state);
                }
                Outcome::Stay
            }
            AppEvent::Key(k) => {
                if let View::Variant(v) = &self.state.view {
                    if !matches!(v.modal, variant::VariantModal::None) {
                        return variant::handle_modal_key(&mut self.state, k);
                    }
                }
                let scope = self.state.active_scope();
                let scope_keys = KeyMap::for_scope(scope);
                if let Some(action) = scope_keys.lookup(k.code, k.modifiers) {
                    return self.dispatch_action(action);
                }
                if let Some(global) = self.global_keys.lookup(k.code, k.modifiers) {
                    return self.handle_global_action(global);
                }
                Outcome::Stay
            }
            _ => Outcome::Stay,
        }
    }

    fn dispatch_action(&mut self, action: Action) -> Outcome {
        match &self.state.view {
            View::Workspace => workspace::handle_action(&mut self.state, action),
            View::Stage(_) => stage_view::handle_action(&mut self.state, action),
            View::Variant(_) => variant::handle_action(&mut self.state, action),
        }
    }

    fn handle_global_action(&mut self, action: Action) -> Outcome {
        match action {
            Action::OpenPalette => {
                let scope = self.state.active_scope();
                self.state.modal = Some(Modal::Palette(Palette::open(scope)));
                Outcome::Stay
            }
            Action::OpenHelp => {
                let scope = self.state.active_scope();
                self.state.modal = Some(Modal::Help(help::HelpState::new(scope)));
                Outcome::Stay
            }
            Action::Quit => Outcome::Quit,
            _ => Outcome::Stay,
        }
    }

    fn dispatch_palette_action(&mut self, action: Action) -> Outcome {
        match action {
            Action::Quit => Outcome::Quit,
            Action::OpenHelp => {
                let scope = self.state.active_scope();
                self.state.modal = Some(Modal::Help(help::HelpState::new(scope)));
                Outcome::Stay
            }
            Action::OpenPalette | Action::ClosePalette | Action::ClosePaletteAndRun => {
                Outcome::Stay
            }
            other => self.dispatch_action(other),
        }
    }

    fn apply_outcome(&mut self, outcome: Outcome) -> bool {
        match outcome {
            Outcome::Stay => false,
            Outcome::Quit => {
                self.state.save_session();
                true
            }
            Outcome::Run(req) => {
                if let RunRequest::Setup(cfg) = &req {
                    if let Some(sink) = self.state.setup_sink.as_ref() {
                        *sink.lock().unwrap() = Some(cfg.clone());
                    }
                }
                let cancel = self.tui_out.arm_cancel();
                self.state.modal = Some(Modal::Run(RunOverlay::launch(&req, cancel)));
                self.spawn_command(req);
                false
            }
        }
    }

    fn spawn_command(&self, req: RunRequest) {
        let out = Arc::clone(&self.tui_out);
        let tx = self.cmd_tx.clone();
        std::thread::spawn(move || {
            let res = dispatch(req, &*out);
            let _ = tx.send(res);
        });
    }

    fn apply_to_log(&mut self, ev: &AppEvent) {
        match ev {
            AppEvent::Log(line) => {
                if let Some(Modal::Run(o)) = &mut self.state.modal {
                    o.push_log(line.clone());
                }
            }
            AppEvent::ProgressUpdate(snap) => {
                if let Some(Modal::Run(o)) = &mut self.state.modal {
                    o.push_progress(snap.clone());
                }
            }
            AppEvent::CommandDone(res) => {
                if let Some(Modal::Run(o)) = &mut self.state.modal {
                    let mirror = match res {
                        Ok(()) => Ok(()),
                        Err(e) => Err(e.to_string()),
                    };
                    o.on_command_done(mirror);
                }
            }
            _ => {}
        }
    }

    fn sample_bars(&mut self) -> Vec<AppEvent> {
        let mut events = Vec::new();
        while let Ok(line) = self.log_rx.try_recv() {
            events.push(AppEvent::Log(line));
        }
        while let Ok(res) = self.cmd_rx.try_recv() {
            events.push(AppEvent::CommandDone(res));
        }
        {
            let mut shared = self.shared_bars.lock().unwrap();
            for arc in shared.drain(..) {
                self.bars.push(BarEntry { arc });
            }
        }
        for entry in &self.bars {
            let pb = &entry.arc;
            let length_opt = pb.length();
            let length = length_opt.unwrap_or(u64::MAX);
            let indeterminate = length_opt.map(|l| l == u64::MAX).unwrap_or(true);
            events.push(AppEvent::ProgressUpdate(ProgressSnapshot {
                position: pb.position(),
                length,
                message: pb.message().to_string(),
                indeterminate,
            }));
        }
        self.bars.retain(|entry| !entry.arc.is_finished());
        events
    }
}

fn render(state: &mut AppState, frame: &mut ratatui::Frame, area: ratatui::layout::Rect) {
    match &state.view {
        View::Workspace => workspace::render(state, frame, area),
        View::Stage(_) => stage_view::render(state, frame, area),
        View::Variant(_) => variant::render(state, frame, area),
    }
}

fn render_modal(state: &mut AppState, frame: &mut ratatui::Frame, area: ratatui::layout::Rect) {
    let Some(modal) = &state.modal else {
        return;
    };
    match modal {
        Modal::Help(h) => help::render(frame, area, h),
        Modal::Palette(p) => p.draw(frame, area),
        Modal::Run(o) => o.render(area, frame.buffer_mut()),
    }
}

fn title_for(state: &AppState) -> &str {
    match &state.view {
        View::Workspace => "Workspace",
        View::Stage(s) => s.title.as_str(),
        View::Variant(v) => v.title.as_str(),
    }
}

fn dispatch(req: RunRequest, out: &dyn Output) -> Result<(), CohortError> {
    match req {
        RunRequest::Ingest(cfg) => commands::ingest::run_ingest(&cfg, out),
        RunRequest::Annotate(cfg) => commands::annotate::run_annotate(&cfg, out),
        RunRequest::Staar(cfg) => StaarPipeline::new(*cfg, out)?.run(),
        RunRequest::MetaStaar(cfg) => commands::meta_staar::run_meta_staar(&cfg, out),
        RunRequest::Setup(cfg) => apply_setup(cfg, out),
    }
}

fn apply_setup(cfg: SetupConfig, out: &dyn Output) -> Result<(), CohortError> {
    let mut existing = crate::config::Config::load().unwrap_or_default();
    existing.data.tier = cfg.tier;
    existing.data.root_dir = cfg.root_dir.to_string_lossy().to_string();
    existing.data.packs = cfg.packs.clone();
    if cfg.environment.is_some() {
        existing.resources.environment = cfg.environment;
    }
    if cfg.memory_budget.is_some() {
        existing.resources.memory_budget = cfg.memory_budget.clone();
    }
    existing.save()?;
    out.success(&format!(
        "saved configuration to {}",
        crate::config::Config::config_path().display()
    ));
    Ok(())
}

pub fn make_terminal() -> io::Result<Terminal<CrosstermBackend<Stdout>>> {
    let backend = CrosstermBackend::new(io::stdout());
    Terminal::new(backend)
}
