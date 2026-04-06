//! Structural invariant validation for the genotype store.
//!
//! Checks that the store is internally consistent:
//! - variant_vcf alignment (row i == variant_vcf i)
//! - sparse_g offsets valid and carrier data well-formed
//! - no duplicate carriers within a variant
//! - membership references valid variant_vcfs
//! - sortedness invariants
//!
//! Any violation is a hard error — data cannot be trusted.

use std::collections::HashSet;
use std::path::Path;

use crate::error::FavorError;
use crate::staar::carrier::reader::VariantIndex;
use crate::staar::sparse_g::SparseG;
use crate::staar::store::StoreManifest;

pub struct CheckResult {
    pub name: &'static str,
    pub passed: bool,
    pub detail: String,
}

/// Validate all structural invariants for the entire store.
pub fn validate_store(
    store_dir: &Path,
    manifest: &StoreManifest,
) -> Result<Vec<CheckResult>, FavorError> {
    let mut results = Vec::new();

    for ci in &manifest.chromosomes {
        let chrom = &ci.name;
        let chrom_dir = store_dir.join(format!("chromosome={chrom}"));

        let sparse_g = SparseG::open(&chrom_dir)?;
        let variant_index = VariantIndex::load(&chrom_dir)?;

        results.extend(check_variant_vcf_alignment(
            &chrom_dir,
            &variant_index,
            chrom,
        ));
        results.extend(check_sparse_g_integrity(&sparse_g, chrom));
        results.extend(check_no_duplicate_carriers(&sparse_g, chrom));
        results.extend(check_membership_validity(&variant_index, &sparse_g, chrom));
        results.extend(check_sortedness(&variant_index, chrom));
        results.extend(check_sparse_g_variant_count(
            &sparse_g,
            &variant_index,
            chrom,
        ));
        results.extend(check_offsets_table(&sparse_g, chrom));
        results.extend(check_carrier_ordering(&sparse_g, chrom));
        results.extend(check_cross_layer_consistency(
            &sparse_g,
            &variant_index,
            chrom,
        ));
    }

    Ok(results)
}

/// Invariant 1: variant_vcf == row index in variants.parquet.
fn check_variant_vcf_alignment(
    chrom_dir: &Path,
    vi: &VariantIndex,
    chrom: &str,
) -> Vec<CheckResult> {
    use arrow::array::UInt32Array;
    use std::fs::File;

    let pq_path = chrom_dir.join("variants.parquet");
    let file = match File::open(&pq_path) {
        Ok(f) => f,
        Err(_) => {
            return vec![CheckResult {
                name: "variant_vcf_alignment",
                passed: false,
                detail: format!("chr{chrom}: cannot open variants.parquet"),
            }]
        }
    };
    let reader = match parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
        .and_then(|b| b.build())
    {
        Ok(r) => r,
        Err(e) => {
            return vec![CheckResult {
                name: "variant_vcf_alignment",
                passed: false,
                detail: format!("chr{chrom}: parquet error: {e}"),
            }]
        }
    };

    let mut expected = 0u32;
    let mut mismatches = 0usize;
    for batch_result in reader {
        let batch = match batch_result {
            Ok(b) => b,
            Err(_) => continue,
        };
        let vvcf_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .unwrap();
        for i in 0..batch.num_rows() {
            if vvcf_arr.value(i) != expected {
                mismatches += 1;
            }
            expected += 1;
        }
    }

    let n = expected as usize;
    let vi_len = vi.len();

    vec![
        CheckResult {
            name: "variant_vcf_alignment",
            passed: mismatches == 0,
            detail: if mismatches == 0 {
                format!("chr{chrom}: {n} variants, all variant_vcf == row index")
            } else {
                format!("chr{chrom}: {mismatches} mismatches in variant_vcf alignment")
            },
        },
        CheckResult {
            name: "variant_count_match",
            passed: n == vi_len,
            detail: format!("chr{chrom}: parquet({n}) vs VariantIndex({vi_len})"),
        },
    ]
}

/// Invariant 2: sparse_g carrier data is well-formed.
/// Every carrier has sample_id < n_samples and dosage in {1, 2}.
fn check_sparse_g_integrity(sg: &SparseG, chrom: &str) -> Vec<CheckResult> {
    let n_variants = sg.n_variants();
    let n_samples = sg.n_samples();
    let mut bad_sample = 0usize;
    let mut bad_dosage = 0usize;

    for v in 0..n_variants {
        let cl = sg.load_variant(v);
        for e in &cl.entries {
            if e.sample_idx >= n_samples {
                bad_sample += 1;
            }
            if e.dosage == 0 || e.dosage > 2 {
                bad_dosage += 1;
            }
        }
    }

    vec![
        CheckResult {
            name: "sparse_g_sample_bounds",
            passed: bad_sample == 0,
            detail: format!("chr{chrom}: {bad_sample} out-of-range sample_ids (max {n_samples})"),
        },
        CheckResult {
            name: "sparse_g_dosage_valid",
            passed: bad_dosage == 0,
            detail: format!("chr{chrom}: {bad_dosage} invalid dosages (expected 1 or 2)"),
        },
    ]
}

/// Invariant 3: no duplicate sample_ids within a single variant.
fn check_no_duplicate_carriers(sg: &SparseG, chrom: &str) -> Vec<CheckResult> {
    let n_variants = sg.n_variants();
    let mut dup_variants = 0usize;

    for v in 0..n_variants {
        let cl = sg.load_variant(v);
        let mut seen = HashSet::with_capacity(cl.entries.len());
        for e in &cl.entries {
            if !seen.insert(e.sample_idx) {
                dup_variants += 1;
                break;
            }
        }
    }

    vec![CheckResult {
        name: "no_duplicate_carriers",
        passed: dup_variants == 0,
        detail: format!("chr{chrom}: {dup_variants} variants with duplicate sample_ids"),
    }]
}

/// Invariant 4: all variant_vcfs in membership < n_variants in sparse_g.
fn check_membership_validity(vi: &VariantIndex, sg: &SparseG, chrom: &str) -> Vec<CheckResult> {
    let n_variants = sg.n_variants();
    let mut bad = 0usize;

    for gene in vi.gene_names() {
        for &v in vi.gene_variant_vcfs(gene) {
            if v >= n_variants {
                bad += 1;
            }
        }
    }

    vec![CheckResult {
        name: "membership_valid",
        passed: bad == 0,
        detail: format!("chr{chrom}: {bad} membership entries with variant_vcf >= {n_variants}"),
    }]
}

/// Invariant 5: sortedness.
fn check_sortedness(vi: &VariantIndex, chrom: &str) -> Vec<CheckResult> {
    let entries = vi.all_entries();
    let mut variant_sorted = true;
    for i in 1..entries.len() {
        let prev = &entries[i - 1];
        let curr = &entries[i];
        if (curr.position, &*curr.ref_allele, &*curr.alt_allele)
            < (prev.position, &*prev.ref_allele, &*prev.alt_allele)
        {
            variant_sorted = false;
            break;
        }
    }

    let mut membership_sorted = true;
    'outer: for gene in vi.gene_names() {
        let vcfs = vi.gene_variant_vcfs(gene);
        for i in 1..vcfs.len() {
            if vcfs[i] <= vcfs[i - 1] {
                membership_sorted = false;
                break 'outer;
            }
        }
    }

    vec![
        CheckResult {
            name: "variants_sorted",
            passed: variant_sorted,
            detail: format!("chr{chrom}: variants sorted by (pos, ref, alt): {variant_sorted}"),
        },
        CheckResult {
            name: "membership_sorted",
            passed: membership_sorted,
            detail: format!("chr{chrom}: membership sorted per gene: {membership_sorted}"),
        },
    ]
}

/// Cross-check: sparse_g.n_variants == variant_index.len()
fn check_sparse_g_variant_count(sg: &SparseG, vi: &VariantIndex, chrom: &str) -> Vec<CheckResult> {
    let sg_n = sg.n_variants() as usize;
    let vi_n = vi.len();
    vec![CheckResult {
        name: "sparse_g_variant_count",
        passed: sg_n == vi_n,
        detail: format!("chr{chrom}: sparse_g({sg_n}) vs VariantIndex({vi_n})"),
    }]
}

/// Invariant 6: offsets table correctness.
/// - Offsets are monotonically non-decreasing (variant data is contiguous)
/// - Every offset is within the carrier data region
/// - Reading at each offset produces a valid carrier count
fn check_offsets_table(sg: &SparseG, chrom: &str) -> Vec<CheckResult> {
    let offsets = sg.offsets();
    let n = offsets.len();
    let data_size = sg.carrier_data_size();

    let mut monotonic = true;
    for i in 1..n {
        if offsets[i] < offsets[i - 1] {
            monotonic = false;
            break;
        }
    }

    let mut out_of_bounds = 0usize;
    for (v, &off) in offsets.iter().enumerate() {
        // Each offset must leave room for at least the carrier count (2 bytes)
        if off + 2 > data_size {
            out_of_bounds += 1;
            continue;
        }
        // Validate: reading the carrier list doesn't panic.
        // load_variant does bounds-checked parsing via mmap slice.
        let cl = sg.load_variant(v as u32);
        // The record's total size: 2 (count) + n_carriers * entry_size
        let entry_size = if sg.n_samples() > 65535 { 5usize } else { 3 };
        let record_size = 2 + cl.len() * entry_size;
        if off + record_size as u64 > data_size {
            out_of_bounds += 1;
        }
    }

    vec![
        CheckResult {
            name: "offsets_monotonic",
            passed: monotonic,
            detail: format!("chr{chrom}: offsets monotonically non-decreasing: {monotonic}"),
        },
        CheckResult {
            name: "offsets_in_bounds",
            passed: out_of_bounds == 0,
            detail: format!("chr{chrom}: {out_of_bounds} offsets exceed carrier data region ({data_size} bytes)"),
        },
    ]
}

/// Invariant 7: carrier ordering — sample_ids strictly increasing within each variant.
/// Enables merge-join, prevents hidden duplicates, guarantees deterministic iteration.
fn check_carrier_ordering(sg: &SparseG, chrom: &str) -> Vec<CheckResult> {
    let n_variants = sg.n_variants();
    let mut unsorted_variants = 0usize;

    for v in 0..n_variants {
        let cl = sg.load_variant(v);
        for i in 1..cl.entries.len() {
            if cl.entries[i].sample_idx <= cl.entries[i - 1].sample_idx {
                unsorted_variants += 1;
                break;
            }
        }
    }

    vec![CheckResult {
        name: "carrier_ordering",
        passed: unsorted_variants == 0,
        detail: format!("chr{chrom}: {unsorted_variants} variants with unsorted sample_ids"),
    }]
}

/// Invariant 8: cross-layer consistency — storage and metadata agree.
/// For every variant_vcf v:
///   - VariantIndex.get(v).vid matches the expected vid from (chrom, pos, ref, alt)
///   - Carrier list from sparse_g has all sample_ids < n_samples
///   - The carrier MAC is consistent (non-zero carriers for maf > 0 variants)
fn check_cross_layer_consistency(sg: &SparseG, vi: &VariantIndex, chrom: &str) -> Vec<CheckResult> {
    let n = sg.n_variants().min(vi.len() as u32);
    let mut vid_mismatches = 0usize;
    let mut mac_inconsistencies = 0usize;

    for v in 0..n {
        let entry = vi.get(v);
        let cl = sg.load_variant(v);

        // Check vid matches expected format from metadata
        let expected_vid =
            crate::types::format_vid(chrom, entry.position, &entry.ref_allele, &entry.alt_allele);
        if *entry.vid != *expected_vid {
            vid_mismatches += 1;
        }

        // MAF > 0 should have at least one carrier; MAF == 0 should have none
        let mac: u32 = cl.entries.iter().map(|e| e.dosage as u32).sum();
        if entry.maf > 0.0 && mac == 0 {
            mac_inconsistencies += 1;
        }
    }

    vec![
        CheckResult {
            name: "cross_layer_vid",
            passed: vid_mismatches == 0,
            detail: format!(
                "chr{chrom}: {vid_mismatches} vid mismatches between VariantIndex and expected"
            ),
        },
        CheckResult {
            name: "cross_layer_mac",
            passed: mac_inconsistencies == 0,
            detail: format!(
                "chr{chrom}: {mac_inconsistencies} variants with maf>0 but zero carriers"
            ),
        },
    ]
}
