use std::path::PathBuf;

use crate::error::CohortError;
use crate::output::Output;

pub fn run(
    _input: PathBuf,
    _tissue: Option<String>,
    _disease: Option<String>,
    _output_path: Option<PathBuf>,
    _output: &dyn Output,
) -> Result<(), CohortError> {
    Err(CohortError::Input(
        "interpret: not yet implemented".to_string(),
    ))
}
