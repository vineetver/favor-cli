//! Memory budget for the dense kinship path.
//!
//! The dense AI-REML loop materializes Σ, Σ⁻¹, Σ⁻¹X and a handful of
//! n×n scratch matrices. Peak working set is roughly `5 · n² · 8` bytes.
//! The budget itself lives on `Resources::kinship_budget_bytes` (resolved
//! once at process start in `resource.rs`); this module is pure math —
//! callers pass the budget explicitly so kinship code never reaches for env.
//!
//! Sparse path users (`SparseColMat<u32, f64>` kinships + faer `Llt`
//! factor) bypass this guard entirely — their working set is `O(nnz_L)`,
//! not `O(n²)`.

use crate::error::CohortError;

/// Default cap exposed for tests so they can pass an explicit budget
/// without round-tripping through `Resources`. Production code reads
/// `Resources::kinship_budget_bytes` instead.
#[cfg(test)]
pub const DEFAULT_KINSHIP_MEM_BYTES: u64 = 16 * (1u64 << 30);

/// Check that dense AI-REML for `n` samples will fit in the configured
/// memory budget. Returns a `Resource` error with a clear hint when the
/// required allocation would exceed the cap.
pub fn check_memory_budget(n: usize, budget_bytes: u64) -> Result<(), CohortError> {
    let needed = bytes_for_dense_reml(n);
    if needed > budget_bytes {
        let needed_gb = needed as f64 / (1u64 << 30) as f64;
        let budget_gb = budget_bytes as f64 / (1u64 << 30) as f64;
        return Err(CohortError::Resource(format!(
            "Kinship-aware REML for n = {n} samples needs ~{needed_gb:.1} GB of dense matrix \
             memory; budget is {budget_gb:.1} GB. Set COHORT_KINSHIP_MEM_GB to raise the cap, \
             pass a sparse pedigree kinship to take the sparse path, or reduce sample count."
        )));
    }
    Ok(())
}

/// Conservative byte estimate for the dense REML working set at `n` samples:
/// `5 · n² · 8` covers Σ, Σ⁻¹, Σ⁻¹X, the AI scratch and a residual buffer.
pub fn bytes_for_dense_reml(n: usize) -> u64 {
    5u64 * (n as u64) * (n as u64) * 8
}

/// Decide between the dense and the sparse REML path. Sparse wins when
/// every kinship is already sparse-stored AND the dense working set would
/// not fit in the budget. A study with all-dense kinships always takes
/// the dense path regardless of size — we'd have nothing to gain by
/// pretending a fully populated GRM is sparse.
pub fn dense_path_fits(n: usize, budget_bytes: u64) -> bool {
    bytes_for_dense_reml(n) <= budget_bytes
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_huge_n() {
        let err = check_memory_budget(200_000, DEFAULT_KINSHIP_MEM_BYTES).unwrap_err();
        match err {
            CohortError::Resource(msg) => {
                assert!(msg.contains("COHORT_KINSHIP_MEM_GB"));
                assert!(msg.contains("200000"));
            }
            other => panic!("expected Resource, got {other:?}"),
        }
    }

    #[test]
    fn accepts_small_n() {
        check_memory_budget(500, DEFAULT_KINSHIP_MEM_BYTES).unwrap();
        check_memory_budget(2000, DEFAULT_KINSHIP_MEM_BYTES).unwrap();
    }
}
