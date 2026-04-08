use ratatui::style::{Color, Modifier, Style};

pub const ACCENT: Color = Color::Cyan;
pub const MUTED: Color = Color::DarkGray;
pub const OK: Color = Color::Green;
pub const WARN: Color = Color::Yellow;
pub const BAD: Color = Color::Red;
pub const FG: Color = Color::White;

pub const FOCUS_GLYPH: &str = "▌";

#[derive(Clone, Copy)]
pub enum Tone {
    Focus,
    Normal,
    Muted,
    Warn,
}

impl Tone {
    pub fn style(self) -> Style {
        match self {
            Tone::Focus => Style::default().fg(ACCENT).add_modifier(Modifier::BOLD),
            Tone::Normal => Style::default().fg(FG),
            Tone::Muted => Style::default().fg(MUTED),
            Tone::Warn => Style::default().fg(WARN),
        }
    }
}

pub fn hint_bar_style() -> Style {
    Style::default().fg(MUTED)
}

pub fn error_slot_style() -> Style {
    Style::default().fg(BAD)
}

pub const GLYPH_VCF: &str = "v";
pub const GLYPH_PARQUET: &str = "p";
pub const GLYPH_PHENO: &str = "P";
pub const GLYPH_KINSHIP: &str = "K";
pub const GLYPH_INGESTED: &str = "I";
pub const GLYPH_ANNOTATED: &str = "A";
pub const GLYPH_GENO_STORE: &str = "G";
pub const GLYPH_STAAR: &str = "S";
pub const GLYPH_ANNOROOT: &str = "R";
