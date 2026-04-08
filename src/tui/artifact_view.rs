use std::path::PathBuf;

use crate::column::Col;

#[derive(Debug, Clone)]
pub enum Predicate {
    Eq(String),
}

#[derive(Debug, Clone)]
pub struct ArtifactFilter {
    pub column: Col,
    pub predicate: Predicate,
}

impl ArtifactFilter {
    pub fn new(column: Col, predicate: Predicate) -> Self {
        Self { column, predicate }
    }
}

pub enum ArtifactOutcome {
    Continue,
    OpenArtifact {
        path: PathBuf,
        filter: Option<ArtifactFilter>,
    },
}

pub trait ArtifactView {
    fn header(&self) -> Option<String> {
        None
    }
    fn outcome(&mut self) -> ArtifactOutcome {
        ArtifactOutcome::Continue
    }
}
