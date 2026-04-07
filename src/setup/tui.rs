//! TUI widgets for the interactive setup wizard.

use std::io;
use std::path::{Path, PathBuf};

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use crossterm::ExecutableCommand;
use ratatui::prelude::*;
use ratatui::widgets::{Block, Borders, List, ListItem, ListState, Padding, Paragraph};

use crate::config::{DirProbe, Environment, ProbeStatus, ResourceConfig, Tier};
use crate::data::Pack;
use crate::resource::Resources;

/// RAII guard for raw mode + alternate screen.
/// Terminal state is always restored, even on panic.
struct TermGuard;

impl TermGuard {
    fn enter() -> io::Result<Self> {
        enable_raw_mode()?;
        io::stdout().execute(EnterAlternateScreen)?;
        Ok(TermGuard)
    }
}

impl Drop for TermGuard {
    fn drop(&mut self) {
        let _ = io::stdout().execute(LeaveAlternateScreen);
        let _ = disable_raw_mode();
    }
}

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

/// Interactive tier selector with live preview panel.
/// Returns the selected Tier, or None if user cancelled.
pub fn select_tier() -> io::Result<Option<Tier>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;
    let mut selected: usize = 0;

    loop {
        terminal.draw(|frame| draw_tier(frame, selected))?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press {
                continue;
            }
            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    selected = selected.saturating_sub(1);
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if selected < TIERS.len() - 1 {
                        selected += 1;
                    }
                }
                KeyCode::Enter => return Ok(Some(TIERS[selected].tier)),
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_tier(frame: &mut Frame, selected: usize) {
    let area = frame.area();

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
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    let mut selector_lines = Vec::new();
    for (i, tier_opt) in TIERS.iter().enumerate() {
        let marker = if i == selected { " > " } else { "   " };
        let style = if i == selected {
            Style::default().fg(Color::Cyan).bold()
        } else {
            Style::default().fg(Color::DarkGray)
        };
        selector_lines.push(Line::from(vec![
            Span::raw(marker),
            Span::styled(format!("{:<6}", tier_opt.tier.as_str()), style),
            Span::styled(
                format!("{:<12}", tier_opt.tier.size_human()),
                if i == selected {
                    Style::default().fg(Color::Yellow)
                } else {
                    Style::default().fg(Color::DarkGray)
                },
            ),
            Span::styled(
                tier_opt.summary,
                if i == selected {
                    Style::default().fg(Color::White)
                } else {
                    Style::default().fg(Color::DarkGray)
                },
            ),
        ]));
    }
    let selector = Paragraph::new(selector_lines).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Select annotation tier ")
            .border_style(Style::default().fg(Color::Cyan)),
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
                    Style::default().fg(Color::DarkGray).italic(),
                ))
            } else {
                let trimmed = l.trim_start();
                if let Some(pos) = trimmed.find("  ") {
                    let (label, rest) = trimmed.split_at(pos);
                    Line::from(vec![
                        Span::styled(
                            format!("  {:<13}", label),
                            Style::default().fg(Color::Yellow),
                        ),
                        Span::styled(rest.trim_start(), Style::default().fg(Color::White)),
                    ])
                } else {
                    Line::from(Span::styled(
                        format!("  {trimmed}"),
                        Style::default().fg(Color::White),
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
            .border_style(Style::default().fg(Color::DarkGray)),
    );
    frame.render_widget(detail, layout[2]);

    let help = Paragraph::new("  up/down navigate    enter select    esc cancel")
        .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[3]);
}

const SELECT_ENTRY: &str = "[ Use this directory ]";

struct DirBrowserState {
    prompt: String,
    current_dir: PathBuf,
    entries: Vec<String>,
    list_state: ListState,
    typing_path: bool,
    input_buf: String,
    probe: DirProbe,
}

impl DirBrowserState {
    fn new(prompt: &str, start: &Path) -> Self {
        let current_dir = if start.is_dir() {
            start.to_path_buf()
        } else {
            start
                .parent()
                .map(|p| p.to_path_buf())
                .unwrap_or_else(|| PathBuf::from("/"))
        };
        let entries = list_dirs(&current_dir);
        let mut list_state = ListState::default();
        if !entries.is_empty() {
            list_state.select(Some(0));
        }
        let probe = DirProbe::scan(&current_dir);

        Self {
            prompt: prompt.to_string(),
            current_dir,
            entries,
            list_state,
            typing_path: false,
            input_buf: String::new(),
            probe,
        }
    }

    fn navigate_to(&mut self, dir: PathBuf) {
        self.current_dir = dir;
        self.entries = list_dirs(&self.current_dir);
        self.list_state.select(if self.entries.is_empty() {
            None
        } else {
            Some(0)
        });
        self.input_buf = self.current_dir.to_string_lossy().to_string();
        self.probe = DirProbe::scan(&self.current_dir);
    }

    fn go_parent(&mut self) {
        if let Some(parent) = self.current_dir.parent() {
            self.navigate_to(parent.to_path_buf());
        }
    }

    fn enter_selected(&mut self) -> Option<PathBuf> {
        let i = self.list_state.selected()?;
        let name = &self.entries[i];

        if name == SELECT_ENTRY {
            return Some(self.current_dir.clone());
        }
        if name == ".." {
            self.go_parent();
            return None;
        }

        let target = self.current_dir.join(name);
        if target.is_dir() {
            self.navigate_to(target);
        }
        None
    }
}

/// Interactive directory browser with live FAVOR data probe panel.
/// Returns the selected path, or None if cancelled.
pub fn select_directory(prompt: &str, default_path: &Path) -> io::Result<Option<PathBuf>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;
    let mut state = DirBrowserState::new(prompt, default_path);

    loop {
        terminal.draw(|frame| draw_dir_browser(frame, &mut state))?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press {
                continue;
            }

            if state.typing_path {
                match key.code {
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
                let typed = Path::new(&state.input_buf);
                if typed.is_dir() {
                    state.probe = DirProbe::scan(typed);
                }
                continue;
            }

            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    if let Some(i) = state.list_state.selected() {
                        if i > 0 {
                            state.list_state.select(Some(i - 1));
                        }
                    }
                }
                KeyCode::Down | KeyCode::Char('j') => {
                    if let Some(i) = state.list_state.selected() {
                        if i + 1 < state.entries.len() {
                            state.list_state.select(Some(i + 1));
                        }
                    }
                }
                KeyCode::Right | KeyCode::Enter => {
                    if let Some(selected) = state.enter_selected() {
                        return Ok(Some(selected));
                    }
                }
                KeyCode::Left | KeyCode::Backspace => {
                    state.go_parent();
                }
                KeyCode::Char(' ') => {
                    return Ok(Some(state.current_dir.clone()));
                }
                KeyCode::Char('c') => {
                    let _ = std::fs::create_dir_all(&state.current_dir);
                    return Ok(Some(state.current_dir.clone()));
                }
                KeyCode::Char('/') | KeyCode::Char('g') => {
                    state.typing_path = true;
                    state.input_buf = state.current_dir.to_string_lossy().to_string();
                }
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn tab_complete(partial: &str) -> Option<String> {
    let path = Path::new(partial);
    let (parent, prefix) = if partial.ends_with('/') {
        if path.is_dir() {
            return list_first_child(path);
        }
        return None;
    } else {
        let parent = path.parent()?;
        let prefix = path.file_name()?.to_string_lossy().to_string();
        (parent, prefix)
    };

    if !parent.is_dir() {
        return None;
    }

    let matches: Vec<String> = std::fs::read_dir(parent)
        .ok()?
        .flatten()
        .filter(|e| e.path().is_dir())
        .filter_map(|e| {
            let name = e.file_name().into_string().ok()?;
            if name.starts_with(&prefix) {
                Some(name)
            } else {
                None
            }
        })
        .collect();

    match matches.len() {
        0 => None,
        1 => {
            let completed = parent.join(&matches[0]);
            Some(format!("{}/", completed.to_string_lossy()))
        }
        _ => {
            let lcp = longest_common_prefix(&matches);
            let completed = parent.join(&lcp);
            Some(completed.to_string_lossy().to_string())
        }
    }
}

fn list_first_child(dir: &Path) -> Option<String> {
    let mut entries: Vec<String> = std::fs::read_dir(dir)
        .ok()?
        .flatten()
        .filter(|e| e.path().is_dir())
        .filter_map(|e| e.file_name().into_string().ok())
        .collect();
    entries.sort();
    if entries.len() == 1 {
        Some(format!("{}/{}/", dir.to_string_lossy(), entries[0]))
    } else {
        None
    }
}

fn longest_common_prefix(strings: &[String]) -> String {
    if strings.is_empty() {
        return String::new();
    }
    let first = &strings[0];
    let mut len = first.len();
    for s in &strings[1..] {
        len = len.min(s.len());
        for (i, (a, b)) in first.bytes().zip(s.bytes()).enumerate() {
            if a != b {
                len = len.min(i);
                break;
            }
        }
    }
    first[..len].to_string()
}

fn list_dirs(dir: &Path) -> Vec<String> {
    let mut dirs = vec![SELECT_ENTRY.to_string(), "..".to_string()];
    if let Ok(read) = std::fs::read_dir(dir) {
        let mut entries: Vec<String> = read
            .flatten()
            .filter(|e| e.path().is_dir())
            .filter_map(|e| e.file_name().into_string().ok())
            .collect();
        entries.sort();
        dirs.extend(entries);
    }
    dirs
}

fn draw_dir_browser(frame: &mut Frame, state: &mut DirBrowserState) {
    let area = frame.area();

    let h_split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(60), Constraint::Percentage(40)])
        .split(area);

    let left_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Length(3),
            Constraint::Min(6),
            Constraint::Length(2),
        ])
        .split(h_split[0]);

    let title = Paragraph::new(format!("  {}", state.prompt))
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, left_layout[0]);

    if state.typing_path {
        let path_valid = Path::new(&state.input_buf).is_dir();
        let path_color = if path_valid { Color::Green } else { Color::Red };
        let path_line = Line::from(vec![
            Span::styled(" > ", Style::default().fg(Color::Yellow)),
            Span::styled(&state.input_buf, Style::default().fg(path_color)),
            Span::styled("_", Style::default().fg(Color::Cyan)),
        ]);
        let hint = if path_valid {
            " Valid directory — enter to go "
        } else {
            " Not a directory "
        };
        let path_block = Paragraph::new(path_line).block(
            Block::default()
                .borders(Borders::ALL)
                .title(hint)
                .border_style(Style::default().fg(if path_valid {
                    Color::Green
                } else {
                    Color::Yellow
                })),
        );
        frame.render_widget(path_block, left_layout[1]);
    } else {
        let path_line = Line::from(vec![
            Span::styled(" ", Style::default()),
            Span::styled(
                state.current_dir.to_string_lossy().to_string(),
                Style::default().fg(Color::White).bold(),
            ),
        ]);
        let path_block = Paragraph::new(path_line).block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Current directory ")
                .border_style(Style::default().fg(Color::Cyan)),
        );
        frame.render_widget(path_block, left_layout[1]);
    }

    let items: Vec<ListItem> = state
        .entries
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let is_selected = state.list_state.selected() == Some(i);
            let (prefix, style) = match name.as_str() {
                s if s == SELECT_ENTRY => (
                    " ",
                    if is_selected {
                        Style::default().fg(Color::Green).bold()
                    } else {
                        Style::default().fg(Color::Green)
                    },
                ),
                ".." => (
                    " ..",
                    if is_selected {
                        Style::default().fg(Color::Cyan).bold()
                    } else {
                        Style::default().fg(Color::DarkGray)
                    },
                ),
                s if s.starts_with('.') => (
                    &name[..],
                    if is_selected {
                        Style::default().fg(Color::DarkGray).bold()
                    } else {
                        Style::default().fg(Color::DarkGray)
                    },
                ),
                _ => (
                    &name[..],
                    if is_selected {
                        Style::default().fg(Color::Cyan).bold()
                    } else {
                        Style::default().fg(Color::White)
                    },
                ),
            };
            let display = if name == SELECT_ENTRY || name == ".." {
                format!(" {prefix}")
            } else {
                format!(" /{name}")
            };
            ListItem::new(display).style(style)
        })
        .collect();

    let list = List::new(items)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Directories ")
                .border_style(Style::default().fg(Color::DarkGray)),
        )
        .highlight_style(Style::default().bg(Color::DarkGray).fg(Color::White))
        .highlight_symbol(" > ");
    frame.render_stateful_widget(list, left_layout[2], &mut state.list_state);

    let help_text = if state.typing_path {
        "  type path    enter go    esc cancel"
    } else {
        "  enter open    space select    / type path    esc cancel"
    };
    let help = Paragraph::new(help_text).style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, left_layout[3]);

    let right_layout = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(6),
            Constraint::Length(3),
        ])
        .split(h_split[1]);

    let probe_header = if state.probe.has_any_data() {
        Paragraph::new("  FAVOR data detected")
            .style(Style::default().fg(Color::Green).bold())
            .block(Block::default().padding(Padding::top(1)))
    } else {
        Paragraph::new("  No FAVOR data at this path")
            .style(Style::default().fg(Color::DarkGray))
            .block(Block::default().padding(Padding::top(1)))
    };
    frame.render_widget(probe_header, right_layout[0]);

    let summary = state.probe.summary_lines();
    let probe_lines: Vec<Line> = summary
        .iter()
        .map(|(text, status)| {
            let (icon, color) = match status {
                ProbeStatus::Good => ("  +  ", Color::Green),
                ProbeStatus::Partial => ("  ~  ", Color::Yellow),
                ProbeStatus::Missing => ("  -  ", Color::DarkGray),
                ProbeStatus::Info => ("     ", Color::DarkGray),
            };
            Line::from(vec![
                Span::styled(icon, Style::default().fg(color)),
                Span::styled(text, Style::default().fg(color)),
            ])
        })
        .collect();

    let probe_detail = Paragraph::new(probe_lines).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Data scan ")
            .border_style(if state.probe.has_any_data() {
                Style::default().fg(Color::Green)
            } else {
                Style::default().fg(Color::DarkGray)
            }),
    );
    frame.render_widget(probe_detail, right_layout[1]);

    let verdict = if state.probe.has_any_data() {
        let tier = state
            .probe
            .detected_tier()
            .map(|t| t.as_str())
            .unwrap_or("?");
        Paragraph::new(Line::from(vec![
            Span::styled("  Ready — ", Style::default().fg(Color::Green)),
            Span::styled(
                format!("{tier} tier detected"),
                Style::default().fg(Color::Green).bold(),
            ),
        ]))
    } else {
        Paragraph::new(Line::from(vec![Span::styled(
            "  Will download data after setup",
            Style::default().fg(Color::Yellow),
        )]))
    };
    frame.render_widget(verdict, right_layout[2]);
}

/// Interactive multi-select for add-on packs.
/// `installed` contains pack IDs already present on disk — shown as pre-checked with "(installed)".
/// Returns list of selected pack IDs, or None if cancelled.
pub fn select_packs(packs: &[&Pack], installed: &[String]) -> io::Result<Option<Vec<String>>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;

    let mut cursor: usize = 0;
    let mut checked: Vec<bool> = packs
        .iter()
        .map(|p| installed.iter().any(|id| id == p.id))
        .collect();

    loop {
        let draw_checked = &checked;
        terminal.draw(|frame| {
            draw_pack_selector(frame, packs, draw_checked, installed, cursor);
        })?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press {
                continue;
            }
            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    cursor = cursor.saturating_sub(1);
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if cursor + 1 < packs.len() {
                        cursor += 1;
                    }
                }
                KeyCode::Char(' ') => {
                    checked[cursor] = !checked[cursor];
                }
                KeyCode::Char('a') => {
                    let all_on = checked.iter().all(|&c| c);
                    for c in checked.iter_mut() {
                        *c = !all_on;
                    }
                }
                KeyCode::Enter => {
                    let selected: Vec<String> = packs
                        .iter()
                        .zip(&checked)
                        .filter(|(_, &on)| on)
                        .map(|(p, _)| p.id.to_string())
                        .collect();
                    return Ok(Some(selected));
                }
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_pack_selector(
    frame: &mut Frame,
    packs: &[&Pack],
    checked: &[bool],
    installed: &[String],
    cursor: usize,
) {
    let area = frame.area();

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
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    let always: Vec<String> = Pack::required()
        .iter()
        .map(|p| format!("{} ({})", p.name, p.size_human))
        .collect();
    let info = Paragraph::new(Line::from(vec![
        Span::styled("  Always installed: ", Style::default().fg(Color::DarkGray)),
        Span::styled(always.join(", "), Style::default().fg(Color::DarkGray)),
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
                Style::default().fg(Color::Green).bold()
            } else if is_cursor {
                Style::default().fg(Color::White)
            } else {
                Style::default().fg(Color::DarkGray)
            };

            let name_style = if is_cursor {
                Style::default().fg(Color::Cyan).bold()
            } else if checked[i] {
                Style::default().fg(Color::White)
            } else {
                Style::default().fg(Color::DarkGray)
            };

            let size_style = if is_cursor {
                Style::default().fg(Color::Yellow)
            } else {
                Style::default().fg(Color::DarkGray)
            };

            let status_tag = if is_installed { " installed " } else { "" };
            let status_style = Style::default().fg(Color::Green);

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
            .border_style(Style::default().fg(Color::Cyan)),
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
    let summary = Paragraph::new(summary_text).style(Style::default().fg(Color::Yellow));
    frame.render_widget(summary, layout[3]);

    let help = Paragraph::new(
        "  up/down navigate    space toggle    a toggle all    enter done    esc skip",
    )
    .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[4]);
}

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

/// Interactive environment selector. Returns None if user skipped/cancelled.
pub fn select_environment() -> io::Result<Option<Environment>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;
    let mut selected: usize = 0;

    loop {
        terminal.draw(|frame| draw_environment(frame, selected))?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press {
                continue;
            }
            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    selected = selected.saturating_sub(1);
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if selected < ENV_OPTIONS.len() - 1 {
                        selected += 1;
                    }
                }
                KeyCode::Enter => return Ok(Some(ENV_OPTIONS[selected].env)),
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_environment(frame: &mut Frame, selected: usize) {
    let area = frame.area();

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
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    let mut selector_lines = Vec::new();
    for (i, opt) in ENV_OPTIONS.iter().enumerate() {
        let marker = if i == selected { " > " } else { "   " };
        let style = if i == selected {
            Style::default().fg(Color::Cyan).bold()
        } else {
            Style::default().fg(Color::DarkGray)
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
            .border_style(Style::default().fg(Color::Cyan)),
    );
    frame.render_widget(selector, layout[1]);

    let desc = Paragraph::new(format!("  {}", ENV_OPTIONS[selected].description))
        .wrap(ratatui::widgets::Wrap { trim: true })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(format!(" {} ", ENV_OPTIONS[selected].label))
                .border_style(Style::default().fg(Color::DarkGray)),
        );
    frame.render_widget(desc, layout[2]);

    let help = Paragraph::new("  up/down navigate    enter select    esc skip")
        .style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[3]);
}

struct MemoryPreset {
    label: &'static str,
    value: &'static str,
    bytes: u64,
}

const MEMORY_PRESETS: &[MemoryPreset] = &[
    MemoryPreset {
        label: "8 GB",
        value: "8GB",
        bytes: 8 * 1024 * 1024 * 1024,
    },
    MemoryPreset {
        label: "16 GB",
        value: "16GB",
        bytes: 16 * 1024 * 1024 * 1024,
    },
    MemoryPreset {
        label: "32 GB",
        value: "32GB",
        bytes: 32 * 1024 * 1024 * 1024,
    },
    MemoryPreset {
        label: "64 GB",
        value: "64GB",
        bytes: 64 * 1024 * 1024 * 1024,
    },
    MemoryPreset {
        label: "128 GB",
        value: "128GB",
        bytes: 128 * 1024 * 1024 * 1024,
    },
    MemoryPreset {
        label: "256 GB",
        value: "256GB",
        bytes: 256 * 1024 * 1024 * 1024,
    },
];

/// Interactive memory budget selector. Shows detected memory as context.
/// Returns None if user skipped, or Some("16GB") etc.
pub fn select_memory_budget(resources: &Resources) -> io::Result<Option<String>> {
    let _guard = TermGuard::enter()?;

    let backend = CrosstermBackend::new(io::stdout());
    let mut terminal = Terminal::new(backend)?;

    let detected = resources.memory_bytes;
    let mut selected: usize = MEMORY_PRESETS
        .iter()
        .rposition(|p| p.bytes <= detected * 100 / 80)
        .unwrap_or(1);

    let custom_idx = MEMORY_PRESETS.len();
    let mut custom_mode = false;
    let mut custom_buf = String::new();

    loop {
        let draw_selected = selected;
        let draw_custom = custom_mode;
        let draw_buf = custom_buf.clone();
        terminal.draw(|frame| {
            draw_memory_budget(
                frame,
                resources,
                draw_selected,
                custom_idx,
                draw_custom,
                &draw_buf,
            );
        })?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press {
                continue;
            }

            if custom_mode {
                match key.code {
                    KeyCode::Enter => {
                        if ResourceConfig::parse_memory_bytes(&custom_buf).is_some() {
                            return Ok(Some(custom_buf));
                        }
                    }
                    KeyCode::Esc => {
                        custom_mode = false;
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
                continue;
            }

            match key.code {
                KeyCode::Up | KeyCode::Char('k') => {
                    selected = selected.saturating_sub(1);
                }
                KeyCode::Down | KeyCode::Char('j') | KeyCode::Tab => {
                    if selected < custom_idx {
                        selected += 1;
                    }
                }
                KeyCode::Enter => {
                    if selected == custom_idx {
                        custom_mode = true;
                        custom_buf.clear();
                    } else {
                        return Ok(Some(MEMORY_PRESETS[selected].value.to_string()));
                    }
                }
                KeyCode::Esc | KeyCode::Char('q') => return Ok(None),
                _ => {}
            }
        }
    }
}

fn draw_memory_budget(
    frame: &mut Frame,
    resources: &Resources,
    selected: usize,
    custom_idx: usize,
    custom_mode: bool,
    custom_buf: &str,
) {
    let area = frame.area();

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
        .style(Style::default().fg(Color::Cyan).bold())
        .block(Block::default().padding(Padding::top(1)));
    frame.render_widget(title, layout[0]);

    let detected_text = format!(
        "  Detected: {}    Environment: {}",
        resources.memory_human(),
        resources.environment(),
    );
    let info = Paragraph::new(Line::from(vec![Span::styled(
        &detected_text,
        Style::default().fg(Color::Yellow),
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
                Style::default().fg(Color::Cyan).bold()
            } else {
                Style::default().fg(Color::White)
            };
            let note = if preset.bytes <= resources.memory_bytes * 100 / 80 {
                Span::styled("  (within detected)", Style::default().fg(Color::Green))
            } else {
                Span::styled("  (exceeds detected)", Style::default().fg(Color::Yellow))
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
        let color = if valid { Color::Green } else { Color::Red };
        items.push(ListItem::new(Line::from(vec![
            Span::styled(" > Custom: ", Style::default().fg(Color::Cyan).bold()),
            Span::styled(custom_buf, Style::default().fg(color)),
            Span::styled("_", Style::default().fg(Color::Cyan)),
            if valid {
                Span::styled("  (valid)", Style::default().fg(Color::Green))
            } else {
                Span::styled(
                    "  (e.g. 48GB, 12288MB)",
                    Style::default().fg(Color::DarkGray),
                )
            },
        ])));
    } else {
        let style = if is_custom_sel {
            Style::default().fg(Color::Cyan).bold()
        } else {
            Style::default().fg(Color::White)
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
            .border_style(Style::default().fg(Color::Cyan)),
    );
    frame.render_widget(list, layout[2]);

    let help_text = if custom_mode {
        "  type amount (e.g. 48GB)    enter confirm    esc cancel"
    } else {
        "  up/down navigate    enter select    esc skip"
    };
    let help = Paragraph::new(help_text).style(Style::default().fg(Color::DarkGray));
    frame.render_widget(help, layout[3]);
}
