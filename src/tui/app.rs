use std::io::{self, Stdout};
use std::sync::mpsc::{self, Receiver, Sender};
use std::sync::Arc;
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
use super::screen::{RunRequest, Screen, Transition};
use super::screens::help::HelpScreen;
use super::screens::run::{RunScreen, RunWiring};
use super::widgets::log_tail::LogTail;
use super::widgets::palette::{Palette, PaletteOutcome};

struct BarEntry {
    id: u64,
    arc: Arc<ProgressBar>,
    finished_once: bool,
}

pub struct App {
    screens: Vec<Box<dyn Screen>>,
    log_rx: Receiver<LogLine>,
    shared_bars: BarRegistry,
    bars: Vec<BarEntry>,
    next_bar_id: u64,
    log: LogTail,
    tui_out: Arc<TuiOutput>,
    cmd_tx: Sender<Result<(), CohortError>>,
    cmd_rx: Receiver<Result<(), CohortError>>,
    global_keys: KeyMap,
    palette: Option<Palette>,
}

impl App {
    pub fn new(
        initial: Box<dyn Screen>,
        log_rx: Receiver<LogLine>,
        shared_bars: BarRegistry,
        tui_out: Arc<TuiOutput>,
    ) -> Self {
        let (cmd_tx, cmd_rx) = mpsc::channel();
        let none = KeyModifiers::NONE;
        let global_keys = KeyMap::new()
            .bind(KeyCode::Char('c'), KeyModifiers::CONTROL, Action::Quit)
            .bind(KeyCode::Char('p'), KeyModifiers::CONTROL, Action::OpenPalette)
            .bind(KeyCode::Char('?'), none, Action::OpenHelp);
        Self {
            screens: vec![initial],
            log_rx,
            shared_bars,
            bars: Vec::new(),
            next_bar_id: 0,
            log: LogTail::new(256),
            tui_out,
            cmd_tx,
            cmd_rx,
            global_keys,
            palette: None,
        }
    }

    pub fn run(mut self, terminal: &mut Terminal<CrosstermBackend<Stdout>>) -> Result<(), CohortError> {
        let tick = Duration::from_millis(50);
        let mut last_title = String::new();
        while let Some(top_idx) = self.screens.len().checked_sub(1) {
            let current_title = self.screens[top_idx].title().to_string();
            if current_title != last_title {
                let _ = io::stdout().execute(SetTitle(format!("cohort — {current_title}")));
                last_title = current_title;
            }
            {
                let log_ref: &LogTail = &self.log;
                let palette_ref = self.palette.as_ref();
                let screen = &mut self.screens[top_idx];
                terminal
                    .draw(|f| {
                        let area = f.area();
                        screen.draw(f, area, log_ref);
                        if let Some(p) = palette_ref {
                            p.draw(f, area);
                        }
                    })
                    .map_err(|e| CohortError::Internal(e.into()))?;
            }

            let mut pending: Vec<AppEvent> = Vec::new();
            if event::poll(tick).map_err(|e| CohortError::Internal(e.into()))? {
                if let Event::Key(k) = event::read().map_err(|e| CohortError::Internal(e.into()))? {
                    if k.kind == KeyEventKind::Press {
                        if k.code == KeyCode::Char('c')
                            && k.modifiers.contains(KeyModifiers::CONTROL)
                        {
                            return Ok(());
                        }
                        if let Some(palette) = self.palette.as_mut() {
                            match palette.handle_key(k) {
                                PaletteOutcome::Stay => {}
                                PaletteOutcome::Close => {
                                    self.palette = None;
                                }
                                PaletteOutcome::Run(action) => {
                                    self.palette = None;
                                    let transition =
                                        self.dispatch_palette_action(top_idx, action);
                                    if self.apply_transition(top_idx, transition) {
                                        continue;
                                    }
                                }
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

            let mut break_loop = false;
            for ev in pending {
                self.apply_to_log(&ev);
                let transition = self.dispatch_event(top_idx, ev);
                if self.apply_transition(top_idx, transition) {
                    break_loop = true;
                    break;
                }
            }
            if break_loop {
                continue;
            }
        }
        Ok(())
    }

    fn dispatch_event(&mut self, top_idx: usize, ev: AppEvent) -> Transition {
        let AppEvent::Key(k) = &ev else {
            return self.screens[top_idx].handle(&ev);
        };
        let screen_bound = self.screens[top_idx]
            .keys()
            .lookup(k.code, k.modifiers)
            .is_some();
        if screen_bound {
            return self.screens[top_idx].handle(&ev);
        }
        if let Some(global_action) = self.global_keys.lookup(k.code, k.modifiers) {
            return self.handle_global_action(top_idx, global_action);
        }
        self.screens[top_idx].handle(&ev)
    }

    fn handle_global_action(&mut self, top_idx: usize, action: Action) -> Transition {
        match action {
            Action::OpenPalette => {
                let scope = self.screens[top_idx].scope();
                self.palette = Some(Palette::open(scope));
                Transition::Stay
            }
            Action::OpenHelp => {
                let scope = self.screens[top_idx].scope();
                Transition::Push(Box::new(HelpScreen::new(scope)))
            }
            Action::Quit => Transition::Quit,
            _ => Transition::Stay,
        }
    }

    fn dispatch_palette_action(&mut self, top_idx: usize, action: Action) -> Transition {
        match action {
            Action::Quit => Transition::Quit,
            Action::OpenHelp => {
                let scope = self.screens[top_idx].scope();
                Transition::Push(Box::new(HelpScreen::new(scope)))
            }
            Action::OpenPalette
            | Action::ClosePalette
            | Action::ClosePaletteAndRun => Transition::Stay,
            other => self.screens[top_idx].on_action(other),
        }
    }

    fn apply_transition(&mut self, _top_idx: usize, transition: Transition) -> bool {
        match transition {
            Transition::Stay => false,
            Transition::Quit => {
                self.screens.clear();
                true
            }
            Transition::Pop => {
                self.screens.pop();
                if let Some(top) = self.screens.last_mut() {
                    top.on_focus();
                }
                true
            }
            Transition::Push(s) => {
                self.screens.push(s);
                true
            }
            Transition::Replace(s) => {
                self.screens.pop();
                self.screens.push(s);
                true
            }
            Transition::Run(req) => {
                let description = req.description();
                let artifact = req.expected_artifact();
                let cancel = self.tui_out.arm_cancel();
                let wiring = RunWiring::Wired { cancel };
                self.spawn_command(req);
                self.screens.push(Box::new(RunScreen::new(description, wiring, artifact)));
                true
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
            AppEvent::Log(line) => self.log.push_log(line.clone()),
            AppEvent::ProgressUpdate(snap) => self.log.push_snapshot(snap.clone()),
            AppEvent::ProgressSweep(id) => self.log.sweep(*id),
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
                let id = self.next_bar_id;
                self.next_bar_id += 1;
                self.bars.push(BarEntry { id, arc, finished_once: false });
            }
        }
        for entry in &self.bars {
            let pb = &entry.arc;
            let length_opt = pb.length();
            let length = length_opt.unwrap_or(u64::MAX);
            let indeterminate = length_opt.map(|l| l == u64::MAX).unwrap_or(true);
            events.push(AppEvent::ProgressUpdate(ProgressSnapshot {
                id: entry.id,
                position: pb.position(),
                length,
                message: pb.message().to_string(),
                indeterminate,
                finished: pb.is_finished(),
            }));
        }
        let mut keep = Vec::with_capacity(self.bars.len());
        for mut entry in self.bars.drain(..) {
            if entry.arc.is_finished() && entry.finished_once {
                events.push(AppEvent::ProgressSweep(entry.id));
            } else {
                if entry.arc.is_finished() {
                    entry.finished_once = true;
                }
                keep.push(entry);
            }
        }
        self.bars = keep;
        events
    }
}

fn dispatch(req: RunRequest, out: &dyn Output) -> Result<(), CohortError> {
    match req {
        RunRequest::Ingest(cfg) => commands::ingest::run_ingest(&cfg, out),
        RunRequest::Annotate(cfg) => commands::annotate::run_annotate(&cfg, out),
    }
}

pub fn make_terminal() -> io::Result<Terminal<CrosstermBackend<Stdout>>> {
    let backend = CrosstermBackend::new(io::stdout());
    Terminal::new(backend)
}
