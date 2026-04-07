#![allow(clippy::needless_range_loop)]
//! Exact selected inversion of Σ via the Takahashi recursion.
//!
//! Replaces the Hutchinson stochastic trace estimator with an exact
//! computation: the Erisman-Tinney recursion walks the Cholesky factor
//! of Σ to fill in the entries of Σ⁻¹ at the L sparsity pattern. Since
//! the L pattern is a superset of the union of all kinship patterns,
//! these entries are exactly what `tr(Σ⁻¹ K_l) = sum(Σ⁻¹ ⊙ K_l)` needs.
//!
//! This is the 1:1 replacement for upstream `R/glmmkin.R::R_fitglmm_ai`'s
//! `sum(Sigma_i * kins[[i]])` Frobenius product. Upstream computes that
//! sum exactly, but only after materializing the *full dense* `Sigma_i`
//! via `chol2inv(chol(Sigma))`. We compute the same sum exactly while
//! storing only the entries of `Sigma_inv` at the union sparsity pattern
//! (orders of magnitude less memory for pedigree kinships).
//!
//! ## The recursion
//!
//! Given the lower Cholesky factor `L` of `Σ` (so `Σ = L Lᵀ`), the entries
//! of `Z = Σ⁻¹ = L⁻ᵀ L⁻¹` satisfy
//!
//! ```text
//! Z[j,j] = 1/L[j,j]² − (1/L[j,j]) Σ_{k>j, (k,j)∈P_L} L[k,j] · Z[k,j]
//! Z[i,j] = -(1/L[j,j]) Σ_{k>j, (k,j)∈P_L} L[k,j] · Z[max(k,i), min(k,i)]
//!          for j < i with (i,j) ∈ P_L
//! ```
//!
//! Iterating `j` from `n-1` down to `0`. At iteration `j`, the lookups
//! `Z[k, *]` for `k > j` have already been filled in by earlier
//! iterations, so the recursion is a single sweep with no fixpoint.
//!
//! ## faer 0.22 plumbing
//!
//! The high-level `Llt<I, T>` wrapper does not expose the L factor's CSC
//! values, so we drop down to the simplicial API:
//!
//! - `prefactorize_symbolic_cholesky` → elimination tree + column counts
//! - `factorize_simplicial_symbolic_cholesky` → `SymbolicSimplicialCholesky`
//!   (exposes `col_ptr`, `row_idx`, `len_val`)
//! - `factorize_simplicial_numeric_llt(&mut L_values, ...)` — we own the
//!   values buffer
//! - `SimplicialLltRef::new(&symbolic, &L_values).solve_in_place_with_conj(...)`
//!   for triangular solves on the same factor
//!
//! Scratch memory is managed via `MemBuffer` / `MemStack` per the
//! `dyn_stack` API. The symbolic factor is built once per `fit_reml`
//! call and reused across every AI iteration; only the numeric factor
//! and the selected inverse are recomputed each iteration.

use faer::dyn_stack::{MemBuffer, MemStack, StackReq};
use faer::linalg::cholesky::llt::factor::LltRegularization;
use faer::reborrow::ReborrowMut;
use faer::sparse::linalg::cholesky::simplicial::{
    factorize_simplicial_numeric_llt, factorize_simplicial_numeric_llt_scratch,
    factorize_simplicial_symbolic_cholesky, factorize_simplicial_symbolic_cholesky_scratch,
    prefactorize_symbolic_cholesky, prefactorize_symbolic_cholesky_scratch,
    SimplicialLltRef, SymbolicSimplicialCholesky,
};
use faer::sparse::{SparseColMat, SparseColMatRef, SymbolicSparseColMatRef};
use faer::{Conj, MatMut, Par};

use crate::error::CohortError;

/// Cached symbolic Cholesky for the Σ pattern. Built once per `fit_reml`
/// call by [`SymbolicTakahashi::new`]; reused across every AI iteration.
pub struct SymbolicTakahashi {
    pub n: usize,
    pub symbolic: SymbolicSimplicialCholesky<u32>,
}

impl SymbolicTakahashi {
    /// Build the symbolic factor from the union pattern of Σ. The
    /// pattern must be square, symmetric (full storage), and contain
    /// the diagonal at every column.
    pub fn new(pattern: SymbolicSparseColMatRef<'_, u32>) -> Result<Self, CohortError> {
        let n = pattern.nrows();
        let nnz = pattern.row_idx().len();

        let mut etree_storage = vec![0_i32; n];
        let mut col_counts = vec![0_u32; n];

        // Pre-factorize: elimination tree + column counts. Then symbolic
        // factorize. Both share one scratch arena.
        let pre_req = StackReq::any_of(&[
            prefactorize_symbolic_cholesky_scratch::<u32>(n, nnz),
            factorize_simplicial_symbolic_cholesky_scratch::<u32>(n),
        ]);
        let mut pre_buf = MemBuffer::try_new(pre_req).map_err(|_| {
            CohortError::Resource("Takahashi symbolic scratch alloc failed".into())
        })?;
        let stack = MemStack::new(&mut pre_buf);

        let etree = prefactorize_symbolic_cholesky::<u32>(
            &mut etree_storage,
            &mut col_counts,
            pattern,
            stack,
        );

        let symbolic = factorize_simplicial_symbolic_cholesky::<u32>(
            pattern,
            etree,
            &col_counts,
            stack,
        )
        .map_err(|e| {
            CohortError::Analysis(format!(
                "Takahashi symbolic Cholesky failed (Σ pattern structurally singular?): {e:?}"
            ))
        })?;

        Ok(Self { n, symbolic })
    }
}

/// Per-iteration numeric state for the Takahashi solver. Owns the L
/// values for the current Σ factorization plus the selected inverse
/// `Z` at the L sparsity pattern.
pub struct TakahashiNumeric {
    pub symbolic: SymbolicSimplicialCholesky<u32>,
    pub l_values: Vec<f64>,
    pub selected: SelectedInverse,
}

impl TakahashiNumeric {
    /// Refactor Σ numerically at the current values, then run the
    /// Takahashi recursion to fill in the selected inverse.
    pub fn factor(
        symbolic: &SymbolicTakahashi,
        sigma: SparseColMatRef<'_, u32, f64>,
    ) -> Result<Self, CohortError> {
        let n = symbolic.n;
        let nnz_l = symbolic.symbolic.len_val();
        let mut l_values = vec![0.0_f64; nnz_l];

        let req = factorize_simplicial_numeric_llt_scratch::<u32, f64>(n);
        let mut buf = MemBuffer::try_new(req).map_err(|_| {
            CohortError::Resource("Takahashi numeric scratch alloc failed".into())
        })?;
        let stack = MemStack::new(&mut buf);

        factorize_simplicial_numeric_llt::<u32, f64>(
            &mut l_values,
            sigma,
            LltRegularization::default(),
            &symbolic.symbolic,
            stack,
        )
        .map_err(|e| {
            CohortError::Analysis(format!(
                "sparse Σ Cholesky failed at current τ (matrix not positive definite): {e:?}"
            ))
        })?;

        let selected = compute_selected_inverse(&symbolic.symbolic, &l_values)?;

        Ok(Self {
            symbolic: symbolic.symbolic.clone(),
            l_values,
            selected,
        })
    }

    /// Apply Σ⁻¹ to `rhs` in place via the cached Cholesky factor.
    pub fn solve_in_place(&self, mut rhs: MatMut<'_, f64>) {
        let llt = SimplicialLltRef::<'_, u32, f64>::new(&self.symbolic, &self.l_values);
        let req = self
            .symbolic
            .solve_in_place_scratch::<f64>(rhs.ncols());
        let mut buf = MemBuffer::try_new(req)
            .expect("solve_in_place scratch alloc failed (always EMPTY for simplicial)");
        let stack = MemStack::new(&mut buf);
        llt.solve_in_place_with_conj(Conj::No, rhs.rb_mut(), Par::Seq, stack);
    }
}

/// Selected inverse `Z` of Σ at the lower-triangular L sparsity pattern.
/// Stored as a CSC matrix where column `j` holds `Z[i, j]` for `i ≥ j`
/// in `(i, j) ∈ P_L`. Symmetry gives the upper triangle.
pub struct SelectedInverse {
    pub n: usize,
    pub col_ptr: Vec<u32>,
    pub row_idx: Vec<u32>,
    pub values: Vec<f64>,
}

impl SelectedInverse {
    /// Read `Z[i, j]`. Returns `0.0` if `(max(i,j), min(i,j))` is not in
    /// the L pattern. Caller is responsible for only querying entries
    /// inside the pattern (the kinship matrices we use are subsets of L
    /// by construction).
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        let (r, c) = if i >= j { (i, j) } else { (j, i) };
        let cs = self.col_ptr[c] as usize;
        let ce = self.col_ptr[c + 1] as usize;
        let slice = &self.row_idx[cs..ce];
        match slice.binary_search(&(r as u32)) {
            Ok(p) => self.values[cs + p],
            Err(_) => 0.0,
        }
    }

    /// Diagonal of Z, length n. Z[j, j] is the first entry of column j
    /// in the CSC layout (sorted ascending by row index).
    pub fn diag(&self) -> Vec<f64> {
        let n = self.n;
        let mut d = vec![0.0_f64; n];
        for j in 0..n {
            let s = self.col_ptr[j] as usize;
            debug_assert_eq!(self.row_idx[s] as usize, j);
            d[j] = self.values[s];
        }
        d
    }

    /// `tr(Z · K) = Σ_{(i,j) ∈ P(K)} Z[i,j] · K[i,j]`. Walks K's
    /// nonzero pattern and reads Z entries by binary search.
    pub fn trace_with(&self, k: &SparseColMat<u32, f64>) -> f64 {
        let n = self.n;
        debug_assert_eq!(k.nrows(), n);
        debug_assert_eq!(k.ncols(), n);
        let mut s = 0.0_f64;
        let cp = k.symbolic().col_ptr();
        let ri = k.symbolic().row_idx();
        let val = k.val();
        for j in 0..n {
            let cs = cp[j] as usize;
            let ce = cp[j + 1] as usize;
            for kk in cs..ce {
                let i = ri[kk] as usize;
                let v = val[kk];
                s += v * self.get(i, j);
            }
        }
        s
    }
}

/// Run the Takahashi recursion on a Cholesky factor `L` (stored in CSC
/// format via `symbolic` and `l_values`) to compute `Z = (L Lᵀ)⁻¹` at
/// the L sparsity pattern.
///
/// Iterates `j` from `n-1` down to `0`. For each column `j`, computes
/// the off-diagonal entries `Z[i, j]` for `i > j` first, then the
/// diagonal `Z[j, j]`. The off-diagonals depend on `Z[k, *]` for `k > j`
/// which have already been computed in earlier iterations.
fn compute_selected_inverse(
    symbolic: &SymbolicSimplicialCholesky<u32>,
    l_values: &[f64],
) -> Result<SelectedInverse, CohortError> {
    let n = symbolic.nrows();
    let col_ptr: Vec<u32> = symbolic.col_ptr().to_vec();
    let row_idx: Vec<u32> = symbolic.row_idx().to_vec();
    let nnz = row_idx.len();
    debug_assert_eq!(l_values.len(), nnz);

    // Sanity check: every column starts with the diagonal.
    for j in 0..n {
        let s = col_ptr[j] as usize;
        if s >= nnz || row_idx[s] as usize != j {
            return Err(CohortError::Internal(anyhow::anyhow!(
                "Cholesky factor column {j} does not start with the diagonal entry"
            )));
        }
    }

    let mut values = vec![0.0_f64; nnz];

    for j in (0..n).rev() {
        let s = col_ptr[j] as usize;
        let e = col_ptr[j + 1] as usize;
        let l_jj = l_values[s];
        if l_jj == 0.0 || !l_jj.is_finite() {
            return Err(CohortError::Internal(anyhow::anyhow!(
                "Cholesky diagonal at column {j} is zero or non-finite: {l_jj}"
            )));
        }
        let inv_l_jj = 1.0 / l_jj;

        // Off-diagonal entries first. For each (i, j) ∈ P_L with i > j,
        // compute Z[i, j] = -(1/L[j,j]) Σ_{kp > j, (kp,j) ∈ P_L} L[kp,j] · Z[max(kp,i), min(kp,i)].
        for i_pos in (s + 1)..e {
            let i_row = row_idx[i_pos] as usize;
            debug_assert!(i_row > j);

            let mut sum = 0.0_f64;
            for kp_pos in (s + 1)..e {
                let kp = row_idx[kp_pos] as usize;
                let l_kp_j = l_values[kp_pos];
                let (lr, lc) = if kp >= i_row { (kp, i_row) } else { (i_row, kp) };
                // Z[lr, lc] is in column lc, which is > j and therefore
                // already computed in an earlier iteration of the outer
                // j-loop. The lookup must succeed because (lr, lc) is in
                // L's pattern by the elimination tree fill rule.
                let cs_lc = col_ptr[lc] as usize;
                let ce_lc = col_ptr[lc + 1] as usize;
                let slice_lc = &row_idx[cs_lc..ce_lc];
                match slice_lc.binary_search(&(lr as u32)) {
                    Ok(p) => sum += l_kp_j * values[cs_lc + p],
                    Err(_) => {
                        return Err(CohortError::Internal(anyhow::anyhow!(
                            "Takahashi precondition violated: Z[{lr}, {lc}] not in L pattern \
                             (kp={kp}, i={i_row}, j={j})"
                        )))
                    }
                }
            }
            values[i_pos] = -inv_l_jj * sum;
        }

        // Diagonal: Z[j, j] = 1/L[j,j]² − (1/L[j,j]) Σ_{k > j, (k,j) ∈ P_L} L[k,j] · Z[k, j].
        // The Z[k, j] needed here are exactly the off-diagonal entries
        // we just stored at positions s+1..e in the values vector.
        let mut sum = 0.0_f64;
        for k_pos in (s + 1)..e {
            sum += l_values[k_pos] * values[k_pos];
        }
        values[s] = inv_l_jj * inv_l_jj - inv_l_jj * sum;
    }

    Ok(SelectedInverse {
        n,
        col_ptr,
        row_idx,
        values,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use faer::prelude::Solve;
    use faer::sparse::Triplet;

    /// 4x4 SPD test matrix, upper-triangle-only storage. This is what
    /// the simplicial Cholesky API expects (see `Assembler::new` for the
    /// reasoning).
    fn small_spd_upper() -> SparseColMat<u32, f64> {
        let entries = vec![
            Triplet::new(0_u32, 0_u32, 4.0),
            Triplet::new(1, 1, 5.0),
            Triplet::new(2, 2, 6.0),
            Triplet::new(3, 3, 7.0),
            Triplet::new(0, 1, 1.0),
            Triplet::new(0, 2, 0.5),
            Triplet::new(1, 3, 0.7),
        ];
        SparseColMat::<u32, f64>::try_new_from_triplets(4, 4, &entries).unwrap()
    }

    /// Same matrix as `small_spd_upper`, in full symmetric storage.
    /// Used by `trace_with` (which expects full storage like the
    /// production kinship matrices) and as the source for the dense
    /// inverse the Takahashi entries are checked against.
    fn small_spd_full() -> SparseColMat<u32, f64> {
        let entries = vec![
            Triplet::new(0_u32, 0_u32, 4.0),
            Triplet::new(1, 1, 5.0),
            Triplet::new(2, 2, 6.0),
            Triplet::new(3, 3, 7.0),
            Triplet::new(0, 1, 1.0),
            Triplet::new(1, 0, 1.0),
            Triplet::new(0, 2, 0.5),
            Triplet::new(2, 0, 0.5),
            Triplet::new(1, 3, 0.7),
            Triplet::new(3, 1, 0.7),
        ];
        SparseColMat::<u32, f64>::try_new_from_triplets(4, 4, &entries).unwrap()
    }

    fn dense_inverse_of(sigma: &SparseColMat<u32, f64>) -> Vec<Vec<f64>> {
        use faer::Mat;
        let n = sigma.nrows();
        let mut dense = Mat::<f64>::zeros(n, n);
        let cp = sigma.symbolic().col_ptr();
        let ri = sigma.symbolic().row_idx();
        let val = sigma.val();
        for j in 0..n {
            let s = cp[j] as usize;
            let e = cp[j + 1] as usize;
            for k in s..e {
                let i = ri[k] as usize;
                dense[(i, j)] = val[k];
            }
        }
        let inv = dense.col_piv_qr().solve(Mat::<f64>::identity(n, n));
        (0..n)
            .map(|i| (0..n).map(|j| inv[(i, j)]).collect())
            .collect()
    }

    #[test]
    fn takahashi_matches_dense_inverse_on_small_spd() {
        let sigma_upper = small_spd_upper();
        let sigma_full = small_spd_full();
        let n = sigma_upper.nrows();
        let symbolic = SymbolicTakahashi::new(sigma_upper.symbolic()).expect("symbolic");
        let numeric = TakahashiNumeric::factor(&symbolic, sigma_upper.as_ref()).expect("numeric");

        let dense_inv = dense_inverse_of(&sigma_full);

        // Every entry stored by the selected inverse must match the
        // exact dense inverse to a tight tolerance.
        for j in 0..n {
            let s = numeric.selected.col_ptr[j] as usize;
            let e = numeric.selected.col_ptr[j + 1] as usize;
            for k in s..e {
                let i = numeric.selected.row_idx[k] as usize;
                let z = numeric.selected.values[k];
                let exact = dense_inv[i][j];
                let err = (z - exact).abs();
                assert!(
                    err < 1e-10,
                    "Z[{i},{j}] = {z}, dense = {exact}, |err| = {err}"
                );
            }
        }
    }

    #[test]
    fn takahashi_diag_matches_dense_diag() {
        let sigma_upper = small_spd_upper();
        let sigma_full = small_spd_full();
        let symbolic = SymbolicTakahashi::new(sigma_upper.symbolic()).expect("symbolic");
        let numeric = TakahashiNumeric::factor(&symbolic, sigma_upper.as_ref()).expect("numeric");
        let diag = numeric.selected.diag();
        let dense_inv = dense_inverse_of(&sigma_full);
        for i in 0..sigma_upper.nrows() {
            let err = (diag[i] - dense_inv[i][i]).abs();
            assert!(err < 1e-10, "diag[{i}] {} vs {}", diag[i], dense_inv[i][i]);
        }
    }

    #[test]
    fn takahashi_trace_with_self_matches_dense() {
        let sigma_upper = small_spd_upper();
        let sigma_full = small_spd_full();
        let symbolic = SymbolicTakahashi::new(sigma_upper.symbolic()).expect("symbolic");
        let numeric = TakahashiNumeric::factor(&symbolic, sigma_upper.as_ref()).expect("numeric");

        // tr(Σ⁻¹ · Σ) = n. `trace_with` expects full symmetric storage
        // (like the production kinship matrices), so we pass the full
        // version here.
        let tr = numeric.selected.trace_with(&sigma_full);
        let n = sigma_upper.nrows() as f64;
        let err = (tr - n).abs();
        assert!(err < 1e-10, "tr(Σ⁻¹ Σ) = {tr}, expected {n}");
    }
}
