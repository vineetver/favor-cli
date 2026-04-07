//! Cached Σ pattern + symbolic Cholesky factor for the sparse path.
//!
//! The sparsity pattern of `Σ = diag(τ_g/W) + Σ_l τ_l · K_l` is constant
//! across AI iterations: it depends only on the union of the K_l patterns
//! plus the diagonal. So we compute it once, build the symbolic Cholesky
//! once, and reuse both each iteration. Per-iteration work is reduced to:
//!
//! 1. assemble the value vector at the current `τ`
//! 2. call `Llt::try_new_with_symbolic` with the cached symbolic factor
//!
//! Both steps are `O(nnz)` and avoid repeating the symbolic analysis,
//! which can dominate the AI step for large pedigrees.

use std::collections::BTreeSet;

use faer::sparse::linalg::solvers::{Llt, SymbolicLlt};
use faer::sparse::{SparseColMatRef, SymbolicSparseColMat};
use faer::Side;

use crate::error::FavorError;
use crate::staar::kinship::types::{GroupPartition, KinshipMatrix, VarianceComponents};

/// Per-entry contribution to Σ at one nonzero position.
///
/// Each entry sums kinship contributions from one or more `K_l[i,j]`
/// values, plus (for diagonal entries) the residual `τ_g(i) / W[i]` term.
#[derive(Clone)]
struct PatternEntry {
    /// `(kinship_idx, K_l[i,j])` pairs that sum into this entry.
    kinship_vals: Vec<(u32, f64)>,
    /// `Some(i)` if this entry is the diagonal at sample `i` and gets the
    /// residual term layered on top of the kinship sum.
    diagonal_sample: Option<u32>,
}

/// Cached pattern + symbolic factor for the sparse Σ rebuild loop.
///
/// Built once per [`fit_reml_sparse`](super::fit_reml_sparse) call. Both
/// the union pattern and the symbolic Cholesky factor are reused across
/// every AI iteration and every boundary refit.
pub struct Assembler {
    n: usize,
    /// CSC column pointer for the union pattern, length `n + 1`.
    pub(super) col_ptr: Vec<u32>,
    /// CSC row indices for the union pattern, length `nnz`.
    pub(super) row_idx: Vec<u32>,
    /// Per-nonzero contribution table, length `nnz`.
    contributions: Vec<PatternEntry>,
    /// Symbolic Cholesky factor of the union pattern. Cloning it bumps
    /// faer's internal `Arc`.
    symbolic: SymbolicLlt<u32>,
}

impl Assembler {
    /// Build the union pattern and the symbolic Cholesky factor. All
    /// kinships must be the sparse variant. Errors if the union pattern
    /// is structurally singular.
    pub fn new(n: usize, kinships: &[KinshipMatrix]) -> Result<Self, FavorError> {
        // Union of K_l patterns + the diagonal.
        let mut col_sets: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); n];
        for j in 0..n {
            col_sets[j].insert(j as u32);
        }
        for k in kinships {
            let sp = k.as_sparse().ok_or_else(|| {
                FavorError::Internal(anyhow::anyhow!(
                    "Assembler::new called with non-sparse kinship"
                ))
            })?;
            let cp = sp.symbolic().col_ptr();
            let ri = sp.symbolic().row_idx();
            for j in 0..n {
                let s = cp[j] as usize;
                let e = cp[j + 1] as usize;
                for kk in s..e {
                    col_sets[j].insert(ri[kk]);
                }
            }
        }

        // Flatten to CSC.
        let mut col_ptr: Vec<u32> = Vec::with_capacity(n + 1);
        let mut row_idx: Vec<u32> = Vec::new();
        col_ptr.push(0);
        for j in 0..n {
            for &r in &col_sets[j] {
                row_idx.push(r);
            }
            col_ptr.push(row_idx.len() as u32);
        }
        let nnz = row_idx.len();

        let mut contributions: Vec<PatternEntry> = (0..nnz)
            .map(|_| PatternEntry {
                kinship_vals: Vec::new(),
                diagonal_sample: None,
            })
            .collect();

        // Diagonal positions get the residual term.
        for j in 0..n {
            let s = col_ptr[j] as usize;
            let e = col_ptr[j + 1] as usize;
            let mut found = false;
            for k in s..e {
                if row_idx[k] as usize == j {
                    contributions[k].diagonal_sample = Some(j as u32);
                    found = true;
                    break;
                }
            }
            if !found {
                return Err(FavorError::Internal(anyhow::anyhow!(
                    "diagonal entry missing in union pattern at column {j}"
                )));
            }
        }

        // Per-kinship contributions, indexed back into the union pattern
        // by binary search (both slices are sorted).
        for (l, k) in kinships.iter().enumerate() {
            let sp = k.as_sparse().unwrap();
            let cp = sp.symbolic().col_ptr();
            let ri = sp.symbolic().row_idx();
            let val = sp.val();
            for j in 0..n {
                let s = cp[j] as usize;
                let e = cp[j + 1] as usize;
                let union_s = col_ptr[j] as usize;
                let union_e = col_ptr[j + 1] as usize;
                for kk in s..e {
                    let row = ri[kk];
                    let v = val[kk];
                    let union_slice = &row_idx[union_s..union_e];
                    let pos = union_slice
                        .binary_search(&row)
                        .map_err(|_| FavorError::Internal(anyhow::anyhow!(
                            "K_l pattern not contained in union pattern"
                        )))?;
                    contributions[union_s + pos]
                        .kinship_vals
                        .push((l as u32, v));
                }
            }
        }

        // Symbolic Cholesky of the union pattern. The symbolic factor is
        // reference-counted by faer; cloning it later is cheap.
        let symbolic_pattern = SymbolicSparseColMat::<u32>::new_checked(
            n,
            n,
            col_ptr.clone(),
            None,
            row_idx.clone(),
        );
        let symbolic = SymbolicLlt::<u32>::try_new(symbolic_pattern.as_ref(), Side::Lower)
            .map_err(|e| {
                FavorError::Analysis(format!(
                    "sparse symbolic Cholesky failed (Σ pattern is structurally singular?): {e:?}"
                ))
            })?;

        Ok(Self {
            n,
            col_ptr,
            row_idx,
            contributions,
            symbolic,
        })
    }

    /// Compute the value vector for Σ at the current τ. Length matches
    /// the union pattern's nnz; positions correspond to `col_ptr` /
    /// `row_idx`.
    pub fn assemble_values(
        &self,
        tau: &VarianceComponents,
        groups: &GroupPartition,
        weights: &[f64],
        row_to_group: &[usize],
    ) -> Vec<f64> {
        let mut values = vec![0.0_f64; self.row_idx.len()];
        for (k, contrib) in self.contributions.iter().enumerate() {
            let mut s = 0.0;
            for &(l, v) in &contrib.kinship_vals {
                s += tau.kinship(l as usize) * v;
            }
            if let Some(i) = contrib.diagonal_sample {
                let i = i as usize;
                let g = row_to_group[i];
                s += tau.group(g) / weights[i];
            }
            values[k] = s;
        }
        let _ = groups;
        values
    }

    /// Factor Σ at the current τ via the cached symbolic factor. Returns
    /// a fresh `Llt`.
    pub fn factor(
        &self,
        tau: &VarianceComponents,
        groups: &GroupPartition,
        weights: &[f64],
        row_to_group: &[usize],
    ) -> Result<Llt<u32, f64>, FavorError> {
        let values = self.assemble_values(tau, groups, weights, row_to_group);
        let symbolic_pattern = SymbolicSparseColMat::<u32>::new_checked(
            self.n,
            self.n,
            self.col_ptr.clone(),
            None,
            self.row_idx.clone(),
        );
        let pattern_ref = symbolic_pattern.as_ref();
        let mat_ref = SparseColMatRef::<'_, u32, f64>::new(pattern_ref, &values);
        Llt::<u32, f64>::try_new_with_symbolic(self.symbolic.clone(), mat_ref, Side::Lower)
            .map_err(|e| {
                FavorError::Analysis(format!(
                    "sparse Σ Cholesky failed at current τ (matrix not positive definite): {e:?}"
                ))
            })
    }

    pub fn n(&self) -> usize {
        self.n
    }
}
