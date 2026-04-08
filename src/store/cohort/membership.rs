use std::collections::HashMap;
use std::fs::File;
use std::path::Path;

use crate::error::CohortError;

pub(super) fn load(chrom_dir: &Path) -> Result<HashMap<String, Vec<u32>>, CohortError> {
    use arrow::array::{StringArray, UInt32Array};

    let pq_path = chrom_dir.join("membership.parquet");
    let file = File::open(&pq_path)
        .map_err(|e| CohortError::Resource(format!("Open {}: {e}", pq_path.display())))?;
    let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
        .map_err(|e| CohortError::Resource(format!("Parquet open: {e}")))?
        .build()
        .map_err(|e| CohortError::Resource(format!("Parquet reader: {e}")))?;

    let mut gene_variants: HashMap<String, Vec<u32>> = HashMap::new();

    for batch_result in reader {
        let batch = batch_result.map_err(|e| CohortError::Resource(format!("Read: {e}")))?;
        let vvcf_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        let gene_arr = batch
            .column(1)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();

        for i in 0..batch.num_rows() {
            let gene = gene_arr.value(i).to_string();
            let vvcf = vvcf_arr.value(i);
            gene_variants.entry(gene).or_default().push(vvcf);
        }
    }

    for variants in gene_variants.values_mut() {
        variants.sort_unstable();
        variants.dedup();
    }

    Ok(gene_variants)
}
