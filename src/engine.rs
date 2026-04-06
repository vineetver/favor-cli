//! DataFusion execution engine.

use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use arrow::array::{self, Array};
use arrow::record_batch::RecordBatch;

use datafusion::common::DataFusionError;
use datafusion::datasource::file_format::parquet::ParquetFormat;
use datafusion::datasource::listing::{
    ListingOptions, ListingTable, ListingTableConfig, ListingTableUrl,
};
use datafusion::execution::memory_pool::{MemoryConsumer, MemoryPool, MemoryReservation};
use datafusion::logical_expr::SortExpr;
use datafusion::prelude::*;

use crate::error::FavorError;
use crate::resource::Resources;

#[derive(Debug)]
pub struct FavorPool {
    limit: usize,
    used: AtomicUsize,
}

impl FavorPool {
    pub fn new(resources: &Resources) -> Self {
        Self {
            limit: resources.memory_bytes as usize,
            used: AtomicUsize::new(0),
        }
    }
}

impl MemoryPool for FavorPool {
    fn grow(&self, _reservation: &MemoryReservation, additional: usize) {
        self.used.fetch_add(additional, Ordering::Relaxed);
    }

    fn shrink(&self, _reservation: &MemoryReservation, shrink: usize) {
        self.used.fetch_sub(shrink, Ordering::Relaxed);
    }

    fn try_grow(
        &self,
        _reservation: &MemoryReservation,
        additional: usize,
    ) -> Result<(), DataFusionError> {
        let old = self.used.fetch_add(additional, Ordering::Relaxed);
        if old + additional > self.limit {
            self.used.fetch_sub(additional, Ordering::Relaxed);
            return Err(DataFusionError::ResourcesExhausted(format!(
                "FavorPool: requested {} bytes, used {}, limit {}",
                additional, old, self.limit
            )));
        }
        Ok(())
    }

    fn reserved(&self) -> usize {
        self.used.load(Ordering::Relaxed)
    }

    fn register(&self, _consumer: &MemoryConsumer) {}
    fn unregister(&self, _consumer: &MemoryConsumer) {}
}

pub struct DfEngine {
    ctx: SessionContext,
    rt: tokio::runtime::Runtime,
}

impl DfEngine {
    pub fn new(resources: &Resources) -> Result<Self, FavorError> {
        let pool = Arc::new(FavorPool::new(resources));

        let config = SessionConfig::new()
            .with_target_partitions(resources.threads)
            .with_batch_size(8192)
            .with_prefer_existing_sort(true);

        let rt_env = datafusion::execution::runtime_env::RuntimeEnvBuilder::new()
            .with_memory_pool(pool)
            .with_temp_file_path(resources.temp_dir.clone())
            .build_arc()
            .map_err(|e| FavorError::Resource(format!("DataFusion runtime init: {e}")))?;

        let ctx = SessionContext::new_with_config_rt(config, rt_env);

        let rt = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(resources.threads.min(4))
            .enable_all()
            .build()
            .map_err(|e| FavorError::Resource(format!("Tokio runtime init: {e}")))?;

        Ok(Self { ctx, rt })
    }

    pub fn register_csv(&self, name: &str, path: &Path, delimiter: u8) -> Result<(), FavorError> {
        let ext = path.extension().and_then(|e| e.to_str()).unwrap_or("csv");
        self.rt.block_on(async {
            let opts = CsvReadOptions::new()
                .has_header(true)
                .delimiter(delimiter)
                .file_extension(ext);
            self.ctx
                .register_csv(name, &path.to_string_lossy(), opts)
                .await
                .map_err(|e| FavorError::Resource(format!("Register CSV '{name}': {e}")))
        })
    }

    pub fn register_parquet_file(&self, name: &str, path: &Path) -> Result<(), FavorError> {
        self.rt.block_on(async {
            self.ctx
                .register_parquet(name, &path.to_string_lossy(), ParquetReadOptions::default())
                .await
                .map_err(|e| FavorError::Resource(format!("Register parquet '{name}': {e}")))
        })
    }

    pub fn register_parquet_dir(&self, name: &str, path: &Path) -> Result<(), FavorError> {
        let url = format!("file://{}", path.display());
        self.rt.block_on(async {
            let table_path = ListingTableUrl::parse(&url)
                .map_err(|e| FavorError::Resource(format!("Invalid path '{url}': {e}")))?;
            let opts = ListingOptions::new(Arc::new(ParquetFormat::default()))
                .with_file_extension(".parquet");
            let schema = opts
                .infer_schema(&self.ctx.state(), &table_path)
                .await
                .map_err(|e| {
                    FavorError::Resource(format!(
                        "Cannot read parquet schema at {}: {e}",
                        path.display()
                    ))
                })?;
            let listing = ListingTable::try_new(
                ListingTableConfig::new(table_path)
                    .with_listing_options(opts)
                    .with_schema(schema),
            )
            .map_err(|e| FavorError::Resource(format!("Listing table '{name}': {e}")))?;
            self.ctx
                .register_table(name, Arc::new(listing))
                .map_err(|e| FavorError::Resource(format!("Register '{name}': {e}")))?;
            Ok(())
        })
    }

    /// Register a parquet file or directory with a declared sort order.
    /// Enables DataFusion to pick SortMergeJoinExec over HashJoinExec
    /// when both sides of a join are sorted on the same key.
    pub fn register_parquet_sorted(
        &self,
        name: &str,
        path: &Path,
        sort_columns: &[&str],
    ) -> Result<(), FavorError> {
        let is_dir = path.is_dir();
        let url = format!("file://{}", path.display());
        self.rt.block_on(async {
            let table_path = ListingTableUrl::parse(&url)
                .map_err(|e| FavorError::Resource(format!("Invalid path '{url}': {e}")))?;

            let sort_order: Vec<Vec<SortExpr>> = vec![sort_columns
                .iter()
                .map(|c| col(*c).sort(true, false))
                .collect()];

            let mut opts = ListingOptions::new(Arc::new(ParquetFormat::default()))
                .with_file_sort_order(sort_order);
            if is_dir {
                opts = opts.with_file_extension(".parquet");
            }

            let schema = opts
                .infer_schema(&self.ctx.state(), &table_path)
                .await
                .map_err(|e| {
                    FavorError::Resource(format!(
                        "Cannot read parquet schema at {}: {e}",
                        path.display()
                    ))
                })?;
            let listing = ListingTable::try_new(
                ListingTableConfig::new(table_path)
                    .with_listing_options(opts)
                    .with_schema(schema),
            )
            .map_err(|e| FavorError::Resource(format!("Listing table '{name}': {e}")))?;
            self.ctx
                .register_table(name, Arc::new(listing))
                .map_err(|e| FavorError::Resource(format!("Register '{name}': {e}")))?;
            Ok(())
        })
    }

    /// Column names from a registered table's Arrow schema — no query needed.
    pub fn table_columns(&self, table_name: &str) -> Result<Vec<String>, FavorError> {
        self.rt.block_on(async {
            let provider = self
                .ctx
                .table_provider(table_name)
                .await
                .map_err(|e| FavorError::Resource(format!("Table '{table_name}': {e}")))?;
            Ok(provider
                .schema()
                .fields()
                .iter()
                .map(|f| f.name().clone())
                .collect())
        })
    }

    pub fn collect(&self, sql: &str) -> Result<Vec<RecordBatch>, FavorError> {
        let sql_preview: String = if sql.len() > 200 {
            format!("{}...", &sql[..200])
        } else {
            sql.to_string()
        };
        self.rt.block_on(async {
            let df = self.ctx.sql(sql).await.map_err(|e| {
                FavorError::Analysis(format!("SQL parse failed: {e}\n  SQL: {sql_preview}"))
            })?;
            df.collect().await.map_err(|e| {
                FavorError::Analysis(format!("Query execution failed: {e}\n  SQL: {sql_preview}"))
            })
        })
    }

    pub fn execute(&self, sql: &str) -> Result<(), FavorError> {
        self.collect(sql)?;
        Ok(())
    }

    /// Return the physical plan as a string for debugging.
    pub fn explain(&self, sql: &str) -> Result<String, FavorError> {
        let batches = self.collect(&format!("EXPLAIN {sql}"))?;
        let mut lines = Vec::new();
        for batch in &batches {
            if let Some(a) = batch
                .column(1)
                .as_any()
                .downcast_ref::<array::StringArray>()
            {
                for i in 0..a.len() {
                    if !a.is_null(i) {
                        lines.push(a.value(i).to_string());
                    }
                }
            }
        }
        Ok(lines.join("\n"))
    }

    /// Returns 0 if query returns no rows.
    pub fn query_scalar(&self, sql: &str) -> Result<i64, FavorError> {
        let batches = self.collect(sql)?;
        for batch in &batches {
            if batch.num_rows() == 0 {
                continue;
            }
            let col = batch.column(0);
            if let Some(a) = col.as_any().downcast_ref::<array::Int64Array>() {
                return Ok(a.value(0));
            }
            if let Some(a) = col.as_any().downcast_ref::<array::UInt64Array>() {
                return Ok(a.value(0) as i64);
            }
            if let Some(a) = col.as_any().downcast_ref::<array::Int32Array>() {
                return Ok(a.value(0) as i64);
            }
        }
        Ok(0)
    }

    pub fn query_strings(&self, sql: &str) -> Result<Vec<String>, FavorError> {
        let batches = self.collect(sql)?;
        let mut result = Vec::new();
        for batch in &batches {
            let col = batch.column(0);
            // DataFusion 53 uses Utf8View for VARCHAR — handle both.
            if let Some(a) = col.as_any().downcast_ref::<array::StringArray>() {
                for i in 0..a.len() {
                    if !a.is_null(i) {
                        result.push(a.value(i).to_string());
                    }
                }
            } else if let Some(a) = col.as_any().downcast_ref::<array::StringViewArray>() {
                for i in 0..a.len() {
                    if !a.is_null(i) {
                        result.push(a.value(i).to_string());
                    }
                }
            }
        }
        Ok(result)
    }
}
