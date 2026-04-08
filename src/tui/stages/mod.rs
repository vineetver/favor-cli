#![allow(dead_code, unused_imports)]

pub mod types;

pub mod annotate;
pub mod ingest;
pub mod meta_staar;
pub mod staar;

pub use types::{
    ArtifactKind, FieldData, FieldValue, FormError, FormField, FormSchema, FormValues, PathKind,
    RunRequest, SessionCtx, StageGroup, StageId,
};

pub trait Stage: Send + Sync + 'static {
    fn id(&self) -> StageId;
    fn label(&self) -> &'static str;
    fn group(&self) -> StageGroup;
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
