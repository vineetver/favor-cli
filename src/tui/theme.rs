use ratatui::style::{Color, Style};

pub const ACCENT: Color = Color::Cyan;
pub const MUTED: Color = Color::DarkGray;
pub const OK: Color = Color::Green;
pub const WARN: Color = Color::Yellow;
pub const BAD: Color = Color::Red;
pub const FG: Color = Color::White;

pub const FOCUS_GLYPH: &str = "▌";

pub fn hint_bar_style() -> Style {
    Style::default().fg(MUTED)
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
