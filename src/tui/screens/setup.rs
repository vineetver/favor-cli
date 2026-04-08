use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use crossterm::event::{KeyCode, KeyModifiers};
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Style, Stylize};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, List, ListItem, Padding, Paragraph};
use ratatui::Frame;

use crate::config::{Config, Environment, ResourceConfig, Tier};
use crate::data::Pack;
use crate::resource::Resources;
use crate::tui::action::{Action, ActionScope, KeyMap};
use crate::tui::event::AppEvent;
use crate::tui::screen::{Screen, Transition};
use crate::tui::theme;
use crate::tui::widgets::file_picker::{self, tab_complete, DirBrowserState};
use crate::tui::widgets::log_tail::LogTail;

struct TierOption {
    tier: Tier,
    summary: &'static str,
    details: &'static [&'static str],
}

const TIERS: &[TierOption] = &[
    TierOption {
        tier: Tier::Base,
        summary: "Curated annotations for most analyses",
        details: &[
            "Genes        GENCODE consequence, transcripts, hgvsc/hgvsp",
            "Clinical     ClinVar significance, disease, review status",
            "Frequency    gnomAD genome+exome (AF only), TOPMed, 1000G 6 pops",
            "Coding       CADD, REVEL, SpliceAI, AlphaMissense, MaveDB",
            "Noncoding    LINSIGHT, FATHMM-XF, GPN-MSA, JARVIS, ReMM, ncER",
            "Conservation PhyloP, PhastCons (3 levels each), GERP, B-statistic",
            "Constraint   gnomAD constraint score + phred",
            "Integrative  13 aPC annotation principal component scores",
            "Regulatory   MACIE, cV2F, cCRE, GeneHancer, CAGE, super enhancers",
            "dbNSFP       15 of 30 predictors (REVEL, MetaSVM, BayesDel, ...)",
            "Other        Mutation rate, distance to TSS/TSE",
            "",
            "Not in base: UCSC, RefSeq, COSMIC, ChromHMM, ENCODE histone marks,",
            "  full gnomAD pops, 15 extra dbNSFP, mappability, variant density",
        ],
    },
    TierOption {
        tier: Tier::Full,
        summary: "Complete FAVOR database - all annotations",
        details: &[
            "Genes        GENCODE + UCSC + RefSeq",
            "Clinical     ClinVar, COSMIC cancer mutations",
            "Frequency    gnomAD 9 pops + sex-stratified + FAF, TOPMed, 1000G",
            "Coding       CADD, REVEL, SpliceAI, AlphaMissense, MaveDB",
            "Noncoding    LINSIGHT, FATHMM-XF, GPN-MSA, JARVIS, ReMM, ncER,",
            "             ncBoost, PGBoost, FunSeq, ALoFT",
            "Conservation PhyloP, PhastCons (3 levels each), GERP, B-statistic",
            "Constraint   gnomAD constraint score + phred",
            "Integrative  13 aPC annotation principal component scores",
            "Regulatory   MACIE, cV2F, cCRE, GeneHancer, CAGE, super enhancers",
            "Epigenomics  ChromHMM 25 states, ENCODE 13 histone marks (raw+phred)",
            "dbNSFP       All 30 predictors",
            "Sequence     Variant density, ReMap, GC/CpG, mappability (4 levels)",
            "Other        Mutation rate, distance, recombination, nucleotide div",
        ],
    },
];

struct EnvOption {
    env: Environment,
    label: &'static str,
    description: &'static str,
}

const ENV_OPTIONS: &[EnvOption] = &[
    EnvOption {
        env: Environment::Hpc,
        label: "HPC cluster",
        description: "Shared cluster with SLURM (srun/sbatch). Memory budget is your \
                       default; srun allocations automatically override it.",
    },
    EnvOption {
        env: Environment::Workstation,
        label: "Workstation",
        description: "Personal machine or dedicated server. Memory budget is the \
                       hard limit for all operations.",
    },
];

struct MemoryPreset {
    label: &'static str,
    value: &'static str,
    bytes: u64,
}

const MEMORY_PRESETS: &[MemoryPreset] = &[
    MemoryPreset { label: "8 GB", value: "8GB", bytes: 8 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "16 GB", value: "16GB", bytes: 16 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "32 GB", value: "32GB", bytes: 32 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "64 GB", value: "64GB", bytes: 64 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "128 GB", value: "128GB", bytes: 128 * 1024 * 1024 * 1024 },
    MemoryPreset { label: "256 GB", value: "256GB", bytes: 256 * 1024 * 1024 * 1024 },
];

#[derive(Debug, Clone)]
pub struct SetupOutcome {
    pub tier: Tier,
    pub root: PathBuf,
    pub packs: Vec<String>,
    pub environment: Option<Environment>,
    pub memory_budget: Option<String>,
}

enum Stage {
    Tier {
        selected: usize,
    },
    Dir(DirBrowserState),
    Packs {
        cursor: usize,
        checked: Vec<bool>,
        installed: Vec<String>,
    },
    Env {
        selected: usize,
    },
    Memory {
        selected: usize,
        custom_mode: bool,
        custom_buf: String,
    },
}

pub type OutcomeSink = Arc<Mutex<Option<SetupOutcome>>>;

pub struct SetupScreen {
    stage: Stage,
    tier: Option<Tier>,
    root: Option<PathBuf>,
    packs: Option<Vec<String>>,
    environment: Option<Environment>,
    memory_budget: Option<String>,
    resources: Resources,
    sink: Option<OutcomeSink>,
}

impl SetupScreen {
    pub fn new() -> Self {
        Self {
            stage: Stage::Tier { selected: 0 },
            tier: None,
            root: None,
            packs: None,
            environment: None,
            memory_budget: None,
            resources: Resources::detect(),
            sink: None,
        }
    }

    pub fn with_sink(
        env: Option<Environment>,
        memory_budget: Option<String>,
        sink: OutcomeSink,
    ) -> Self {
        let mut s = Self::new();
        s.environment = env;
        s.memory_budget = memory_budget;
        s.sink = Some(sink);
        s
    }

    fn finish(&mut self) {
        let Some(sink) = self.sink.as_ref() else {
            return;
        };
        if let (Some(tier), Some(root)) = (self.tier, self.root.clone()) {
            *sink.lock().unwrap() = Some(SetupOutcome {
                tier,
                root,
                packs: self.packs.clone().unwrap_or_default(),
                environment: self.environment,
                memory_budget: self.memory_budget.clone(),
            });
        }
    }

    fn enter_dir_stage(&mut self) {
        let cwd = std::env::current_dir().unwrap_or_else(|_| Config::default_root_dir());
        self.stage = Stage::Dir(DirBrowserState::new(
            "Select FAVOR data root directory",
            &cwd,
        ));
    }

    fn enter_packs_stage(&mut self) {
        let root = self.root.clone().unwrap_or_else(Config::default_root_dir);
        let optional = Pack::optional();
        let installed: Vec<String> = optional
            .iter()
            .filter(|p| p.tables.iter().any(|t| p.local_dir(&root).join(t).is_dir()))
            .map(|p| p.id.to_string())
            .collect();
        let checked: Vec<bool> = optional
            .iter()
            .map(|p| installed.iter().any(|id| id == p.id))
            .collect();
        self.stage = Stage::Packs { cursor: 0, checked, installed };
    }

    fn enter_env_stage(&mut self) {
        if self.environment.is_some() {
            self.enter_memory_stage();
            return;
        }
        self.stage = Stage::Env { selected: 0 };
    }

    fn enter_memory_stage(&mut self) {
        if self.memory_budget.is_some() {
            self.finish();
            return;
        }
        let detected = self.resources.memory_bytes;
        let selected = MEMORY_PRESETS
            .iter()
            .rposition(|p| p.bytes <= detected * 100 / 80)
            .unwrap_or(1);
        self.stage = Stage::Memory {
            selected,
            custom_mode: false,
            custom_buf: String::new(),
        };
    }

    fn handle_tier(&mut self, code: KeyCode) -> Transition {
        let Stage::Tier { selected } = &mut self.stage else {
            return Transition::Stay;
        };
        match code {
            KeyCode::Up | KeyCode::Char('k') => {
                *selected = selected.saturating_sub(1);
            }
            KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                if *selected < TIERS.len() - 1 {
                    *selected += 1;
                }
            }
            KeyCode::Enter => {
                self.tier = Some(TIERS[*selected].tier);
                self.enter_dir_stage();
            }
            KeyCode::Esc | KeyCode::Char('q') => {
                return Transition::Pop;
            }
            _ => {}
        }
        Transition::Stay
    }

    fn handle_dir(&mut self, code: KeyCode) -> Transition {
        let Stage::Dir(state) = &mut self.stage else {
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
                KeyCode::Char(c) => {
                    state.input_buf.push(c);
                }
                KeyCode::Tab => {
                    if let Some(completed) = tab_complete(&state.input_buf) {
                        state.input_buf = completed;
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
            KeyCode::Right | KeyCode::Enter => {
                if let Some(selected) = state.enter_selected() {
                    self.root = Some(selected);
                    self.enter_packs_stage();
                }
            }
            KeyCode::Left | KeyCode::Backspace => state.go_parent(),
            KeyCode::Char(' ') => {
                self.root = Some(state.current_dir.clone());
                self.enter_packs_stage();
            }
            KeyCode::Char('c') => {
                let _ = std::fs::create_dir_all(&state.current_dir);
                self.root = Some(state.current_dir.clone());
                self.enter_packs_stage();
            }
            KeyCode::Char('/') | KeyCode::Char('g') => {
                state.typing_path = true;
                state.input_buf = state.current_dir.to_string_lossy().to_string();
            }
            KeyCode::Esc | KeyCode::Char('q') => {
                return Transition::Pop;
            }
            _ => {}
        }
        Transition::Stay
    }

    fn handle_packs(&mut self, code: KeyCode) -> Transition {
        let Stage::Packs { cursor, checked, .. } = &mut self.stage else {
            return Transition::Stay;
        };
        let optional_count = Pack::optional().len();
        match code {
            KeyCode::Up | KeyCode::Char('k') => {
                *cursor = cursor.saturating_sub(1);
            }
            KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                if *cursor + 1 < optional_count {
                    *cursor += 1;
                }
            }
            KeyCode::Char(' ') => {
                checked[*cursor] = !checked[*cursor];
            }
            KeyCode::Char('a') => {
                let all_on = checked.iter().all(|&c| c);
                for c in checked.iter_mut() {
                    *c = !all_on;
                }
            }
            KeyCode::Enter => {
                let selected: Vec<String> = Pack::optional()
                    .into_iter()
                    .zip(checked.iter())
                    .filter(|(_, &on)| on)
                    .map(|(p, _)| p.id.to_string())
                    .collect();
                self.packs = Some(selected);
                self.enter_env_stage();
            }
            KeyCode::Esc | KeyCode::Char('q') => {
                self.packs = Some(Vec::new());
                self.enter_env_stage();
            }
            _ => {}
        }
        Transition::Stay
    }

    fn handle_env(&mut self, code: KeyCode) -> Transition {
        let Stage::Env { selected } = &mut self.stage else {
            return Transition::Stay;
        };
        match code {
            KeyCode::Up | KeyCode::Char('k') => {
                *selected = selected.saturating_sub(1);
            }
            KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                if *selected < ENV_OPTIONS.len() - 1 {
                    *selected += 1;
                }
            }
            KeyCode::Enter => {
                self.environment = Some(ENV_OPTIONS[*selected].env);
                self.enter_memory_stage();
            }
            KeyCode::Esc | KeyCode::Char('q') => {
                self.environment = None;
                self.enter_memory_stage();
            }
            _ => {}
        }
        Transition::Stay
    }

    fn handle_memory(&mut self, code: KeyCode) -> Transition {
        let Stage::Memory { selected, custom_mode, custom_buf } = &mut self.stage else {
            return Transition::Stay;
        };
        let custom_idx = MEMORY_PRESETS.len();

        if *custom_mode {
            match code {
                KeyCode::Enter => {
                    if ResourceConfig::parse_memory_bytes(custom_buf).is_some() {
                        self.memory_budget = Some(custom_buf.clone());
                        self.finish();
                        return Transition::Pop;
                    }
                }
                KeyCode::Esc => {
                    *custom_mode = false;
                    custom_buf.clear();
                }
                KeyCode::Backspace => {
                    custom_buf.pop();
                }
                KeyCode::Char(c) => {
                    custom_buf.push(c);
                }
                _ => {}
            }
            return Transition::Stay;
        }

        match code {
            KeyCode::Up | KeyCode::Char('k') => {
                *selected = selected.saturating_sub(1);
            }
            KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                if *selected < custom_idx {
                    *selected += 1;
                }
            }
            KeyCode::Enter => {
                if *selected == custom_idx {
                    *custom_mode = true;
                    custom_buf.clear();
                } else {
                    self.memory_budget = Some(MEMORY_PRESETS[*selected].value.to_string());
                    self.finish();
                    return Transition::Pop;
                }
            }
            KeyCode::Esc | KeyCode::Char('q') => {
                self.memory_budget = None;
                self.finish();
                return Transition::Pop;
            }
            _ => {}
        }
        Transition::Stay
    }

    fn draw_tier(frame: &mut Frame, area: Rect, selected: usize) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Length(5),
                Constraint::Min(10),
                Constraint::Length(2),
            ])
            .split(area);

        let title = Paragraph::new("  COHORT Setup — Annotation Tier")
            .style(Style::default().fg(theme::ACCENT).bold())
            .block(Block::default().padding(Padding::top(1)));
        frame.render_widget(title, layout[0]);

        let mut selector_lines = Vec::new();
        for (i, tier_opt) in TIERS.iter().enumerate() {
            let marker = if i == selected { " > " } else { "   " };
            let style = if i == selected {
                Style::default().fg(theme::ACCENT).bold()
            } else {
                Style::default().fg(theme::MUTED)
            };
            selector_lines.push(Line::from(vec![
                Span::raw(marker),
                Span::styled(format!("{:<6}", tier_opt.tier.as_str()), style),
                Span::styled(
                    format!("{:<12}", tier_opt.tier.size_human()),
                    if i == selected {
                        Style::default().fg(theme::WARN)
                    } else {
                        Style::default().fg(theme::MUTED)
                    },
                ),
                Span::styled(
                    tier_opt.summary,
                    if i == selected {
                        Style::default().fg(theme::FG)
                    } else {
                        Style::default().fg(theme::MUTED)
                    },
                ),
            ]));
        }
        let selector = Paragraph::new(selector_lines).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Select annotation tier ")
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(selector, layout[1]);

        let tier_opt = &TIERS[selected];
        let detail_lines: Vec<Line> = tier_opt
            .details
            .iter()
            .map(|l| {
                if l.is_empty() {
                    Line::from("")
                } else if l.starts_with("Not in") {
                    Line::from(Span::styled(
                        format!("  {l}"),
                        Style::default().fg(theme::MUTED).italic(),
                    ))
                } else {
                    let trimmed = l.trim_start();
                    if let Some(pos) = trimmed.find("  ") {
                        let (label, rest) = trimmed.split_at(pos);
                        Line::from(vec![
                            Span::styled(
                                format!("  {:<13}", label),
                                Style::default().fg(theme::WARN),
                            ),
                            Span::styled(rest.trim_start(), Style::default().fg(theme::FG)),
                        ])
                    } else {
                        Line::from(Span::styled(
                            format!("  {trimmed}"),
                            Style::default().fg(theme::FG),
                        ))
                    }
                }
            })
            .collect();

        let detail_title = format!(
            " {} - {} ",
            tier_opt.tier.as_str(),
            tier_opt.tier.size_human()
        );
        let detail = Paragraph::new(detail_lines).block(
            Block::default()
                .borders(Borders::ALL)
                .title(detail_title)
                .border_style(Style::default().fg(theme::MUTED)),
        );
        frame.render_widget(detail, layout[2]);

        let help = Paragraph::new("  up/down navigate    enter select    esc cancel")
            .style(Style::default().fg(theme::MUTED));
        frame.render_widget(help, layout[3]);
    }

    fn draw_packs(
        frame: &mut Frame,
        area: Rect,
        cursor: usize,
        checked: &[bool],
        installed: &[String],
    ) {
        let packs: Vec<&Pack> = Pack::optional();
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Length(3),
                Constraint::Min(8),
                Constraint::Length(1),
                Constraint::Length(2),
            ])
            .split(area);

        let title = Paragraph::new("  COHORT Setup — Add-on Packs")
            .style(Style::default().fg(theme::ACCENT).bold())
            .block(Block::default().padding(Padding::top(1)));
        frame.render_widget(title, layout[0]);

        let always: Vec<String> = Pack::required()
            .iter()
            .map(|p| format!("{} ({})", p.name, p.size_human))
            .collect();
        let info = Paragraph::new(Line::from(vec![
            Span::styled("  Always installed: ", Style::default().fg(theme::MUTED)),
            Span::styled(always.join(", "), Style::default().fg(theme::MUTED)),
        ]))
        .block(Block::default().padding(Padding::top(1)));
        frame.render_widget(info, layout[1]);

        let items: Vec<ListItem> = packs
            .iter()
            .enumerate()
            .map(|(i, pack)| {
                let is_cursor = i == cursor;
                let is_installed = installed.iter().any(|id| id == pack.id);
                let mark = if checked[i] { "[x]" } else { "[ ]" };

                let mark_style = if checked[i] {
                    Style::default().fg(theme::OK).bold()
                } else if is_cursor {
                    Style::default().fg(theme::FG)
                } else {
                    Style::default().fg(theme::MUTED)
                };

                let name_style = if is_cursor {
                    Style::default().fg(theme::ACCENT).bold()
                } else if checked[i] {
                    Style::default().fg(theme::FG)
                } else {
                    Style::default().fg(theme::MUTED)
                };

                let size_style = if is_cursor {
                    Style::default().fg(theme::WARN)
                } else {
                    Style::default().fg(theme::MUTED)
                };

                let status_tag = if is_installed { " installed " } else { "" };
                let status_style = Style::default().fg(theme::OK);

                let line = Line::from(vec![
                    Span::styled(format!("  {mark} "), mark_style),
                    Span::styled(format!("{:<16}", pack.id), name_style),
                    Span::styled(format!("{:>6}  ", pack.size_human), size_style),
                    Span::styled(status_tag, status_style),
                    Span::styled(pack.description, name_style),
                ]);
                ListItem::new(line)
            })
            .collect();

        let list = List::new(items).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Select packs (space toggle) ")
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(list, layout[2]);

        let selected_count = checked.iter().filter(|&&c| c).count();
        let selected_bytes: u64 = packs
            .iter()
            .zip(checked)
            .filter(|(_, &on)| on)
            .map(|(p, _)| p.size_bytes)
            .sum();
        let total_gb = selected_bytes as f64 / (1024.0 * 1024.0 * 1024.0);
        let summary_text = if selected_count > 0 {
            format!("  Selected: {selected_count} packs ({total_gb:.0} GB)")
        } else {
            "  No packs selected (you can add them later with `cohort data pull --pack <name>`)"
                .to_string()
        };
        let summary = Paragraph::new(summary_text).style(Style::default().fg(theme::WARN));
        frame.render_widget(summary, layout[3]);

        let help = Paragraph::new(
            "  up/down navigate    space toggle    a toggle all    enter done    esc skip",
        )
        .style(Style::default().fg(theme::MUTED));
        frame.render_widget(help, layout[4]);
    }

    fn draw_env(frame: &mut Frame, area: Rect, selected: usize) {
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Length(6),
                Constraint::Min(6),
                Constraint::Length(2),
            ])
            .split(area);

        let title = Paragraph::new("  COHORT Setup — Environment")
            .style(Style::default().fg(theme::ACCENT).bold())
            .block(Block::default().padding(Padding::top(1)));
        frame.render_widget(title, layout[0]);

        let mut selector_lines = Vec::new();
        for (i, opt) in ENV_OPTIONS.iter().enumerate() {
            let marker = if i == selected { " > " } else { "   " };
            let style = if i == selected {
                Style::default().fg(theme::ACCENT).bold()
            } else {
                Style::default().fg(theme::MUTED)
            };
            selector_lines.push(Line::from(vec![
                Span::raw(marker),
                Span::styled(opt.label, style),
            ]));
        }
        let selector = Paragraph::new(selector_lines).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Are you on an HPC cluster or a workstation? ")
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(selector, layout[1]);

        let desc = Paragraph::new(format!("  {}", ENV_OPTIONS[selected].description))
            .wrap(ratatui::widgets::Wrap { trim: true })
            .block(
                Block::default()
                    .borders(Borders::ALL)
                    .title(format!(" {} ", ENV_OPTIONS[selected].label))
                    .border_style(Style::default().fg(theme::MUTED)),
            );
        frame.render_widget(desc, layout[2]);

        let help = Paragraph::new("  up/down navigate    enter select    esc skip")
            .style(Style::default().fg(theme::MUTED));
        frame.render_widget(help, layout[3]);
    }

    fn draw_memory(
        frame: &mut Frame,
        area: Rect,
        resources: &Resources,
        selected: usize,
        custom_mode: bool,
        custom_buf: &str,
    ) {
        let custom_idx = MEMORY_PRESETS.len();
        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints([
                Constraint::Length(3),
                Constraint::Length(3),
                Constraint::Min(10),
                Constraint::Length(2),
            ])
            .split(area);

        let title = Paragraph::new("  COHORT Setup — Memory Budget")
            .style(Style::default().fg(theme::ACCENT).bold())
            .block(Block::default().padding(Padding::top(1)));
        frame.render_widget(title, layout[0]);

        let detected_text = format!(
            "  Detected: {}    Environment: {}",
            resources.memory_human(),
            resources.environment(),
        );
        let info = Paragraph::new(Line::from(vec![Span::styled(
            &detected_text,
            Style::default().fg(theme::WARN),
        )]))
        .block(Block::default().padding(Padding::top(1)));
        frame.render_widget(info, layout[1]);

        let mut items: Vec<ListItem> = MEMORY_PRESETS
            .iter()
            .enumerate()
            .map(|(i, preset)| {
                let is_sel = i == selected && !custom_mode;
                let marker = if is_sel { " > " } else { "   " };
                let style = if is_sel {
                    Style::default().fg(theme::ACCENT).bold()
                } else {
                    Style::default().fg(theme::FG)
                };
                let note = if preset.bytes <= resources.memory_bytes * 100 / 80 {
                    Span::styled("  (within detected)", Style::default().fg(theme::OK))
                } else {
                    Span::styled("  (exceeds detected)", Style::default().fg(theme::WARN))
                };
                ListItem::new(Line::from(vec![
                    Span::raw(marker),
                    Span::styled(format!("{:<10}", preset.label), style),
                    note,
                ]))
            })
            .collect();

        let is_custom_sel = selected == custom_idx && !custom_mode;
        if custom_mode {
            let valid = ResourceConfig::parse_memory_bytes(custom_buf).is_some();
            let color = if valid { theme::OK } else { theme::BAD };
            items.push(ListItem::new(Line::from(vec![
                Span::styled(" > Custom: ", Style::default().fg(theme::ACCENT).bold()),
                Span::styled(custom_buf, Style::default().fg(color)),
                Span::styled("_", Style::default().fg(theme::ACCENT)),
                if valid {
                    Span::styled("  (valid)", Style::default().fg(theme::OK))
                } else {
                    Span::styled("  (e.g. 48GB, 12288MB)", Style::default().fg(theme::MUTED))
                },
            ])));
        } else {
            let style = if is_custom_sel {
                Style::default().fg(theme::ACCENT).bold()
            } else {
                Style::default().fg(theme::FG)
            };
            let marker = if is_custom_sel { " > " } else { "   " };
            items.push(ListItem::new(Line::from(vec![
                Span::raw(marker),
                Span::styled("Custom...", style),
            ])));
        }

        let list = List::new(items).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Default memory budget (srun allocations override this) ")
                .border_style(Style::default().fg(theme::ACCENT)),
        );
        frame.render_widget(list, layout[2]);

        let help_text = if custom_mode {
            "  type amount (e.g. 48GB)    enter confirm    esc cancel"
        } else {
            "  up/down navigate    enter select    esc skip"
        };
        let help = Paragraph::new(help_text).style(Style::default().fg(theme::MUTED));
        frame.render_widget(help, layout[3]);
    }
}

impl Screen for SetupScreen {
    fn title(&self) -> &str {
        "Setup"
    }

    fn scope(&self) -> ActionScope {
        match &self.stage {
            Stage::Dir(_) => ActionScope::FilePicker,
            _ => ActionScope::Setup,
        }
    }

    fn keys(&self) -> KeyMap {
        let none = KeyModifiers::NONE;
        let mut map = KeyMap::new()
            .bind(KeyCode::Esc, none, Action::SetupCancel)
            .bind(KeyCode::Char('q'), none, Action::SetupCancel);
        match &self.stage {
            Stage::Tier { .. } | Stage::Env { .. } => {
                map = map
                    .bind(KeyCode::Up, none, Action::SetupPrev)
                    .bind(KeyCode::Char('k'), none, Action::SetupPrev)
                    .bind(KeyCode::Down, none, Action::SetupNext)
                    .bind(KeyCode::Char('j'), none, Action::SetupNext)
                    .bind(KeyCode::Tab, none, Action::SetupNext)
                    .bind(KeyCode::Enter, none, Action::SetupConfirm);
            }
            Stage::Packs { .. } => {
                map = map
                    .bind(KeyCode::Up, none, Action::SetupPrev)
                    .bind(KeyCode::Char('k'), none, Action::SetupPrev)
                    .bind(KeyCode::Down, none, Action::SetupNext)
                    .bind(KeyCode::Char('j'), none, Action::SetupNext)
                    .bind(KeyCode::Tab, none, Action::SetupNext)
                    .bind(KeyCode::Char(' '), none, Action::SetupToggle)
                    .bind(KeyCode::Char('a'), none, Action::SetupToggleAll)
                    .bind(KeyCode::Enter, none, Action::SetupConfirm);
            }
            Stage::Memory { .. } => {
                map = map
                    .bind(KeyCode::Up, none, Action::SetupPrev)
                    .bind(KeyCode::Char('k'), none, Action::SetupPrev)
                    .bind(KeyCode::Down, none, Action::SetupNext)
                    .bind(KeyCode::Char('j'), none, Action::SetupNext)
                    .bind(KeyCode::Tab, none, Action::SetupNext)
                    .bind(KeyCode::Enter, none, Action::SetupConfirm);
            }
            Stage::Dir(_) => {
                map = KeyMap::new()
                    .bind(KeyCode::Esc, none, Action::PickerCancel)
                    .bind(KeyCode::Char('q'), none, Action::PickerCancel)
                    .bind(KeyCode::Up, none, Action::PickerUp)
                    .bind(KeyCode::Char('k'), none, Action::PickerUp)
                    .bind(KeyCode::Down, none, Action::PickerDown)
                    .bind(KeyCode::Char('j'), none, Action::PickerDown)
                    .bind(KeyCode::Right, none, Action::PickerInto)
                    .bind(KeyCode::Enter, none, Action::PickerInto)
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

    fn draw(&mut self, frame: &mut Frame, area: Rect, _log: &LogTail) {
        match &mut self.stage {
            Stage::Tier { selected } => Self::draw_tier(frame, area, *selected),
            Stage::Dir(state) => file_picker::draw(frame, area, state),
            Stage::Packs { cursor, checked, installed } => {
                Self::draw_packs(frame, area, *cursor, checked, installed)
            }
            Stage::Env { selected } => Self::draw_env(frame, area, *selected),
            Stage::Memory { selected, custom_mode, custom_buf } => Self::draw_memory(
                frame,
                area,
                &self.resources,
                *selected,
                *custom_mode,
                custom_buf,
            ),
        }
    }

    // SetupScreen owns key dispatch directly because the Dir stage's path-typing
    // mode and the Memory stage's custom-entry mode swallow every character into a
    // text buffer; routing those through KeyMap would require modeling each
    // character as an Action. keys() above is informational only — used by help
    // and the command palette to enumerate the bindings for the active stage.
    fn handle(&mut self, event: &AppEvent) -> Transition {
        let AppEvent::Key(k) = event else {
            return Transition::Stay;
        };
        match self.stage {
            Stage::Tier { .. } => self.handle_tier(k.code),
            Stage::Dir(_) => self.handle_dir(k.code),
            Stage::Packs { .. } => self.handle_packs(k.code),
            Stage::Env { .. } => self.handle_env(k.code),
            Stage::Memory { .. } => self.handle_memory(k.code),
        }
    }
}
