use std::fmt;

#[derive(Debug)]
pub enum CohortError {
    Input(String),
    DataMissing(String),
    Resource(String),
    Analysis(String),
    Internal(anyhow::Error),
    Cancelled,
}

impl CohortError {
    pub fn exit_code(&self) -> i32 {
        match self {
            Self::Input(_) => 1,
            Self::DataMissing(_) => 2,
            Self::Resource(_) => 3,
            Self::Analysis(_) => 4,
            Self::Internal(_) => 5,
            Self::Cancelled => 130,
        }
    }

    pub fn code_name(&self) -> &'static str {
        match self {
            Self::Input(_) => "input_error",
            Self::DataMissing(_) => "data_missing",
            Self::Resource(_) => "resource_error",
            Self::Analysis(_) => "analysis_error",
            Self::Internal(_) => "internal_error",
            Self::Cancelled => "cancelled",
        }
    }

    /// Prepend a context label to the error message while preserving the
    /// variant (and therefore the exit code). Used at parallel-worker
    /// boundaries where we want to attach `(worker_id, file path)` to an
    /// otherwise opaque error without losing its classification.
    pub fn with_context(self, ctx: impl std::fmt::Display) -> Self {
        match self {
            Self::Input(m) => Self::Input(format!("{ctx}: {m}")),
            Self::DataMissing(m) => Self::DataMissing(format!("{ctx}: {m}")),
            Self::Resource(m) => Self::Resource(format!("{ctx}: {m}")),
            Self::Analysis(m) => Self::Analysis(format!("{ctx}: {m}")),
            Self::Internal(e) => Self::Internal(e.context(ctx.to_string())),
            Self::Cancelled => Self::Cancelled,
        }
    }
}

impl fmt::Display for CohortError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Input(msg)
            | Self::DataMissing(msg)
            | Self::Resource(msg)
            | Self::Analysis(msg) => write!(f, "{msg}"),
            Self::Internal(err) => write!(f, "{err}"),
            Self::Cancelled => write!(f, "cancelled by user"),
        }
    }
}

impl From<anyhow::Error> for CohortError {
    fn from(err: anyhow::Error) -> Self {
        Self::Internal(err)
    }
}

impl From<std::io::Error> for CohortError {
    fn from(err: std::io::Error) -> Self {
        Self::Internal(err.into())
    }
}

impl From<toml::de::Error> for CohortError {
    fn from(err: toml::de::Error) -> Self {
        Self::Internal(err.into())
    }
}

impl From<toml::ser::Error> for CohortError {
    fn from(err: toml::ser::Error) -> Self {
        Self::Internal(err.into())
    }
}
