pub mod types;

pub mod annotate;
pub mod ingest;
pub mod meta_staar;
pub mod staar;

pub use types::{
    ArtifactKind, FormError, FormSchema, FormValues, RunRequest, SessionCtx, StageId,
};

pub trait Stage: Send + Sync + 'static {
    fn id(&self) -> StageId;
    fn label(&self) -> &'static str;
    fn inputs(&self) -> &'static [ArtifactKind];
    fn outputs(&self) -> &'static [ArtifactKind];
    fn form_schema(&self, ctx: &SessionCtx) -> FormSchema;
    fn build_command(&self, values: &FormValues) -> Result<RunRequest, FormError>;
}

use self::annotate::AnnotateStage;
use self::ingest::IngestStage;
use self::meta_staar::MetaStaarStage;
use self::staar::StaarStage;

pub const STAGES: &[&dyn Stage] = &[
    &IngestStage,
    &AnnotateStage,
    &StaarStage,
    &MetaStaarStage,
];

pub fn next_for(kind: ArtifactKind) -> Option<&'static dyn Stage> {
    STAGES
        .iter()
        .copied()
        .find(|s| s.inputs().contains(&kind))
}
