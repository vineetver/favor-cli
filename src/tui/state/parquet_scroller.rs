use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use arrow::array::RecordBatch;
use arrow::datatypes::SchemaRef;
use parquet::arrow::arrow_reader::{
    ArrowReaderMetadata, ArrowReaderOptions, ParquetRecordBatchReaderBuilder, RowFilter,
};
use parquet::arrow::ProjectionMask;

use crate::error::CohortError;

pub trait RowFilterFactory: Send + Sync {
    fn build(&self, meta: &ArrowReaderMetadata) -> Result<RowFilter, CohortError>;
    fn describe(&self) -> &str;
}

pub struct ParquetScroller {
    path: PathBuf,
    meta: ArrowReaderMetadata,
    schema: SchemaRef,
    projection: ProjectionMask,
    filter: Option<Arc<dyn RowFilterFactory>>,
    row_group_count: usize,
    rg_idx: usize,
    rg_row: usize,
    batch: Option<RecordBatch>,
}

impl ParquetScroller {
    pub fn open(path: &Path) -> Result<Self, CohortError> {
        let file = File::open(path).map_err(|e| {
            CohortError::Resource(format!("open {}: {e}", path.display()))
        })?;
        let meta = ArrowReaderMetadata::load(&file, ArrowReaderOptions::new())
            .map_err(|e| CohortError::Input(format!("parquet metadata {}: {e}", path.display())))?;
        let schema = meta.schema().clone();
        let row_group_count = meta.metadata().num_row_groups();
        let projection = ProjectionMask::all();

        let mut s = Self {
            path: path.to_path_buf(),
            meta,
            schema,
            projection,
            filter: None,
            row_group_count,
            rg_idx: 0,
            rg_row: 0,
            batch: None,
        };
        if row_group_count > 0 {
            s.load_current_row_group()?;
        }
        Ok(s)
    }

    pub fn schema(&self) -> &SchemaRef {
        &self.schema
    }

    pub fn metadata(&self) -> &ArrowReaderMetadata {
        &self.meta
    }

    pub fn row_group_count(&self) -> usize {
        self.row_group_count
    }

    pub fn current_row_group(&self) -> usize {
        self.rg_idx
    }

    pub fn current_batch_len(&self) -> usize {
        self.batch.as_ref().map(|b| b.num_rows()).unwrap_or(0)
    }

    pub fn focus_row(&self) -> Option<usize> {
        let n = self.current_batch_len();
        if n == 0 {
            None
        } else {
            Some(self.rg_row.min(n - 1))
        }
    }

    pub fn focused_record(&self) -> Option<(&RecordBatch, usize)> {
        let row = self.focus_row()?;
        self.batch.as_ref().map(|b| (b, row))
    }

    pub fn filter_text(&self) -> Option<&str> {
        self.filter.as_ref().map(|f| f.describe())
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn goto_row_group(&mut self, i: usize) -> Result<(), CohortError> {
        if self.row_group_count == 0 {
            return Ok(());
        }
        let i = i.min(self.row_group_count - 1);
        self.rg_idx = i;
        self.load_current_row_group()
    }

    pub fn next_row_group(&mut self) -> Result<(), CohortError> {
        if self.row_group_count == 0 {
            return Ok(());
        }
        if self.rg_idx + 1 < self.row_group_count {
            self.rg_idx += 1;
            self.load_current_row_group()
        } else {
            Ok(())
        }
    }

    pub fn prev_row_group(&mut self) -> Result<(), CohortError> {
        if self.rg_idx == 0 {
            return Ok(());
        }
        self.rg_idx -= 1;
        self.load_current_row_group()
    }

    pub fn scroll_up(&mut self) {
        self.rg_row = self.rg_row.saturating_sub(1);
    }

    pub fn scroll_down(&mut self) {
        let n = self.current_batch_len();
        if n > 0 && self.rg_row + 1 < n {
            self.rg_row += 1;
        }
    }

    pub fn set_projection(&mut self, mask: ProjectionMask) -> Result<(), CohortError> {
        self.projection = mask;
        if self.row_group_count > 0 {
            self.load_current_row_group()?;
        }
        Ok(())
    }

    pub fn set_filter(
        &mut self,
        filter: Option<Arc<dyn RowFilterFactory>>,
    ) -> Result<(), CohortError> {
        self.filter = filter;
        if self.row_group_count > 0 {
            self.load_current_row_group()?;
        }
        Ok(())
    }

    pub fn total_rows_estimate(&self) -> u64 {
        self.meta
            .metadata()
            .row_groups()
            .iter()
            .map(|rg| rg.num_rows() as u64)
            .sum()
    }

    fn load_current_row_group(&mut self) -> Result<(), CohortError> {
        let file = File::open(&self.path).map_err(|e| {
            CohortError::Resource(format!("open {}: {e}", self.path.display()))
        })?;
        let mut builder =
            ParquetRecordBatchReaderBuilder::new_with_metadata(file, self.meta.clone())
                .with_row_groups(vec![self.rg_idx])
                .with_projection(self.projection.clone());
        if let Some(factory) = self.filter.clone() {
            let row_filter = factory.build(&self.meta)?;
            builder = builder.with_row_filter(row_filter);
        }
        let mut reader = builder
            .build()
            .map_err(|e| CohortError::Analysis(format!("parquet build: {e}")))?;

        let mut combined: Option<RecordBatch> = None;
        for batch in reader.by_ref() {
            let b = batch.map_err(|e| CohortError::Analysis(format!("parquet read: {e}")))?;
            combined = Some(match combined.take() {
                None => b,
                Some(prev) => arrow::compute::concat_batches(&prev.schema(), &[prev, b])
                    .map_err(|e| CohortError::Analysis(format!("concat: {e}")))?,
            });
        }
        self.batch = combined;
        self.rg_row = 0;
        Ok(())
    }
}
