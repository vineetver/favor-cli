//! Format detection framework.

use std::path::Path;

use arrow::datatypes::{DataType, Field, Schema};

use parquet::file::reader::FileReader;

use crate::error::CohortError;

use super::{sniff_delimiter, Delimiter, InputFormat};

/// Detection confidence: 0.0 (no match) to 1.0 (certain).
pub enum DetectResult {
    Yes(f32),
    No,
}

/// What the registry returns after detection.
pub struct Detected {
    pub format: InputFormat,
    #[allow(dead_code)] // available for callers that need delimiter after detection
    pub delimiter: Option<Delimiter>,
    pub handler_name: &'static str,
}

/// A handler for one input file format.
#[allow(dead_code)] // public API — schema() and register() used as extension points
pub trait FormatHandler: Send + Sync {
    fn name(&self) -> &'static str;

    /// Can this handler process the given file?
    /// `header` is the first 8KB of the file.
    fn detect(&self, path: &Path, header: &[u8]) -> DetectResult;

    /// What columns does this format produce?
    fn schema(&self, path: &Path) -> Result<Schema, CohortError>;

    /// The InputFormat enum variant this handler maps to.
    fn input_format(&self) -> InputFormat;

    /// For tabular formats, the detected delimiter (if any).
    fn delimiter(&self, path: &Path) -> Option<Delimiter>;
}

fn path_lower(path: &Path) -> String {
    path.file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .to_lowercase()
}

fn first_line(header: &[u8]) -> Option<&str> {
    let text = std::str::from_utf8(header).ok()?;
    text.lines().next()
}

pub struct VcfHandler;

impl FormatHandler for VcfHandler {
    fn name(&self) -> &'static str {
        "VCF"
    }

    fn detect(&self, path: &Path, header: &[u8]) -> DetectResult {
        let name = path_lower(path);
        if name.ends_with(".vcf")
            || name.ends_with(".vcf.gz")
            || name.ends_with(".vcf.bgz")
            || name.ends_with(".bcf")
        {
            return DetectResult::Yes(0.9);
        }
        // Content-based: check magic bytes
        if header.starts_with(b"##fileformat=VCF") {
            return DetectResult::Yes(1.0);
        }
        // BGZF magic: 1f 8b 08 04 (gzip with extra field)
        if header.len() >= 4 && header[0] == 0x1f && header[1] == 0x8b {
            return DetectResult::Yes(0.5);
        }
        DetectResult::No
    }

    fn schema(&self, _path: &Path) -> Result<Schema, CohortError> {
        Ok(Schema::new(vec![
            Field::new("chromosome", DataType::Utf8, false),
            Field::new("position", DataType::Int32, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
            Field::new("rsid", DataType::Utf8, true),
            Field::new("qual", DataType::Utf8, true),
            Field::new("filter", DataType::Utf8, true),
        ]))
    }

    fn input_format(&self) -> InputFormat {
        InputFormat::Vcf
    }
    fn delimiter(&self, _path: &Path) -> Option<Delimiter> {
        None
    }
}

pub struct ParquetHandler;

impl FormatHandler for ParquetHandler {
    fn name(&self) -> &'static str {
        "Parquet"
    }

    fn detect(&self, path: &Path, header: &[u8]) -> DetectResult {
        if header.starts_with(b"PAR1") {
            return DetectResult::Yes(1.0);
        }
        let name = path_lower(path);
        if name.ends_with(".parquet") {
            return DetectResult::Yes(0.9);
        }
        DetectResult::No
    }

    fn schema(&self, path: &Path) -> Result<Schema, CohortError> {
        let file = std::fs::File::open(path)
            .map_err(|e| CohortError::Input(format!("Cannot open '{}': {e}", path.display())))?;
        let reader = parquet::file::reader::SerializedFileReader::new(file)
            .map_err(|e| CohortError::Input(format!("Bad parquet '{}': {e}", path.display())))?;
        let parquet_schema = reader.metadata().file_metadata().schema_descr();
        let fields: Vec<Field> = parquet_schema
            .root_schema()
            .get_fields()
            .iter()
            .map(|f| Field::new(f.name(), DataType::Utf8, true))
            .collect();
        Ok(Schema::new(fields))
    }

    fn input_format(&self) -> InputFormat {
        InputFormat::Parquet
    }
    fn delimiter(&self, _path: &Path) -> Option<Delimiter> {
        None
    }
}

pub struct TsvHandler;

impl FormatHandler for TsvHandler {
    fn name(&self) -> &'static str {
        "TSV"
    }

    fn detect(&self, path: &Path, header: &[u8]) -> DetectResult {
        let name = path_lower(path);
        if name.ends_with(".tsv")
            || name.ends_with(".tsv.gz")
            || name.ends_with(".txt")
            || name.ends_with(".txt.gz")
        {
            return DetectResult::Yes(0.9);
        }
        if let Some(line) = first_line(header) {
            if line.matches('\t').count() >= 2 {
                return DetectResult::Yes(0.6);
            }
        }
        DetectResult::No
    }

    fn schema(&self, path: &Path) -> Result<Schema, CohortError> {
        let delim = sniff_delimiter(path)?;
        let headers = super::read_headers(path, delim)?;
        let fields: Vec<Field> = headers
            .iter()
            .map(|h| Field::new(h.as_str(), DataType::Utf8, true))
            .collect();
        Ok(Schema::new(fields))
    }

    fn input_format(&self) -> InputFormat {
        InputFormat::Tabular
    }

    fn delimiter(&self, path: &Path) -> Option<Delimiter> {
        sniff_delimiter(path).ok()
    }
}

pub struct CsvHandler;

impl FormatHandler for CsvHandler {
    fn name(&self) -> &'static str {
        "CSV"
    }

    fn detect(&self, path: &Path, header: &[u8]) -> DetectResult {
        let name = path_lower(path);
        if name.ends_with(".csv") || name.ends_with(".csv.gz") {
            return DetectResult::Yes(0.9);
        }
        if let Some(line) = first_line(header) {
            let commas = line.matches(',').count();
            let tabs = line.matches('\t').count();
            if commas >= 2 && commas > tabs {
                return DetectResult::Yes(0.5);
            }
        }
        DetectResult::No
    }

    fn schema(&self, path: &Path) -> Result<Schema, CohortError> {
        let headers = super::read_headers(path, Delimiter::Comma)?;
        let fields: Vec<Field> = headers
            .iter()
            .map(|h| Field::new(h.as_str(), DataType::Utf8, true))
            .collect();
        Ok(Schema::new(fields))
    }

    fn input_format(&self) -> InputFormat {
        InputFormat::Tabular
    }

    fn delimiter(&self, _path: &Path) -> Option<Delimiter> {
        Some(Delimiter::Comma)
    }
}

pub struct FormatRegistry {
    handlers: Vec<Box<dyn FormatHandler>>,
}

impl FormatRegistry {
    pub fn new() -> Self {
        Self {
            handlers: vec![
                Box::new(VcfHandler),
                Box::new(ParquetHandler),
                Box::new(TsvHandler),
                Box::new(CsvHandler),
            ],
        }
    }

    /// Detect the best handler for a file using extension + content sniffing.
    pub fn detect(&self, path: &Path) -> Result<Detected, CohortError> {
        let mut header = [0u8; 8192];
        let n = std::fs::File::open(path)
            .and_then(|mut f| std::io::Read::read(&mut f, &mut header))
            .map_err(|e| CohortError::Input(format!("Cannot read '{}': {e}", path.display())))?;

        let mut best: Option<(&dyn FormatHandler, f32)> = None;
        for handler in &self.handlers {
            if let DetectResult::Yes(conf) = handler.detect(path, &header[..n]) {
                if best.is_none() || conf > best.unwrap().1 {
                    best = Some((handler.as_ref(), conf));
                }
            }
        }

        match best {
            Some((h, _)) => Ok(Detected {
                format: h.input_format(),
                delimiter: h.delimiter(path),
                handler_name: h.name(),
            }),
            None => Err(CohortError::Input(format!(
                "Cannot detect format for '{}'. Supported formats: {}",
                path.display(),
                self.handlers
                    .iter()
                    .map(|h| h.name())
                    .collect::<Vec<_>>()
                    .join(", ")
            ))),
        }
    }

    #[allow(dead_code)]
    pub fn register(&mut self, handler: Box<dyn FormatHandler>) {
        self.handlers.push(handler);
    }
}
