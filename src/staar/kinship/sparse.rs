//! Sparse path for kinship-aware AI-REML.
//!
//! When all loaded kinships are sparse and dense Σ⁻¹ would blow the memory
//! budget, we route through this module instead of the dense `dense.rs`
//! kernel. The internals:
//!
//! - **Σ assembly**: build a sparse `SparseColMat<u32, f64>` from the kinship
//!   sparsity pattern + the residual diagonal in one pass per AI iteration.
//!   The pattern is invariant in τ, so the symbolic Cholesky is computed
//!   once outside the loop and cached.
//!
//! - **Σ⁻¹·v**: faer `Llt::solve_in_place`. Cost is `O(nnz_L)` per RHS.
//!
//! - **Trace term `tr(Σ⁻¹ K_l)`**: two strategies, picked by `nnz_L` and
//!   the user-configurable probe count.
//!   - [`hutchinson_trace`]: stochastic, deterministic seed, used as the
//!     general fallback. ~3% relative error with M=30 probes.
//!   - Takahashi selected inversion (lands as a follow-up): exact,
//!     bit-identical to GMMAT, preferred when `nnz_L` is moderate.
//!
//! - **Score path**: `KinshipState` carries a `SparseFactor` which the
//!   sparse-aware score kernel uses to compute Σ⁻¹·G_j on demand per
//!   variant column.
//!
//! This file is the entry point for everything sparse. Dense path lives in
//! `dense.rs`; the public dispatcher lives in `mod.rs`.

use std::collections::BTreeSet;

use faer::linalg::solvers::Solve;
use faer::sparse::linalg::solvers::{Llt, SymbolicLlt};
use faer::sparse::{SparseColMat, SparseColMatRef, SymbolicSparseColMat};
use faer::{Mat, Side};

use crate::error::FavorError;
use crate::staar::kinship::dense::{
    weights_or_one, BOUNDARY_FACTOR, REML_MAX_ITER, REML_MAX_OUTER_REFITS, REML_STEP_HALVE_MAX,
    REML_TOL,
};
use crate::staar::kinship::types::{
    GroupPartition, KinshipInverse, KinshipMatrix, KinshipState, VarianceComponents,
};

/// Default number of Hutchinson probes for the trace estimator. M=30 gives
/// ~3% relative error and reproducible results when seeded.
pub const DEFAULT_HUTCHINSON_PROBES: usize = 30;
/// Fixed seed for the Rademacher PRNG. Same seed → same trace estimates →
/// run-to-run reproducibility.
pub const HUTCHINSON_SEED: u64 = 0x5F3759DF_5F3759DF;

/// Owned wrapper around the faer sparse Cholesky factor of Σ. Cloned into
/// `KinshipState::inverse` (the `Sparse` arm) when fit_reml takes the
/// sparse path. Once held, every `Σ⁻¹·v` operation in scoring becomes
/// `factor.solve_in_place(v)`.
#[derive(Clone)]
pub struct SparseFactor {
    /// Faer's high-level sparse `LLᵀ` factor.
    pub llt: Llt<u32, f64>,
}

// ─── Σ pattern + assembler ──────────────────────────────────────────────────

/// Per-entry contribution to Σ at one nonzero position. Most entries pull
/// from a single kinship; a few (the diagonal) also pull from the residual
/// τ_g/W[i] term.
#[derive(Clone)]
struct EntryContrib {
    /// `(kinship_idx, K_l[i,j])` pairs that sum into this entry.
    kinship_vals: Vec<(u32, f64)>,
    /// If `Some(i)`, this entry is the diagonal at sample `i` and gets the
    /// residual τ_{g(i)}/W[i] contribution layered on top of the kinship sum.
    diagonal_sample: Option<u32>,
}

/// Cached pattern + symbolic factor for the sparse Σ rebuild loop.
///
/// The sparsity pattern of `Σ = τ_e·I + Σ_l τ_l · K_l` is constant across
/// AI iterations because it depends only on the union of the K_l patterns
/// plus the diagonal. We compute it once, build the symbolic Cholesky once,
/// and reuse both across every iteration. Per-iteration work: produce a
/// `Vec<f64>` of values and call `Llt::try_new_with_symbolic`.
pub struct SparseSigmaAssembler {
    n: usize,
    /// CSC column pointer for the union pattern, length `n + 1`.
    col_ptr: Vec<u32>,
    /// CSC row indices for the union pattern, length `nnz`.
    row_idx: Vec<u32>,
    /// Per-nonzero contribution table, length `nnz`.
    contributions: Vec<EntryContrib>,
    /// Symbolic Cholesky factor of the union pattern. Cloning it is an
    /// `Arc` bump (faer reference-counts the symbolic factor).
    symbolic: SymbolicLlt<u32>,
}

impl SparseSigmaAssembler {
    /// Build the union pattern and the symbolic Cholesky factor. All
    /// kinships must be the sparse variant.
    pub fn new(n: usize, kinships: &[KinshipMatrix]) -> Result<Self, FavorError> {
        // Union of K_l patterns + diagonal.
        let mut col_sets: Vec<BTreeSet<u32>> = vec![BTreeSet::new(); n];
        for j in 0..n {
            col_sets[j].insert(j as u32);
        }
        for k in kinships {
            let sp = k.as_sparse().ok_or_else(|| {
                FavorError::Internal(anyhow::anyhow!(
                    "SparseSigmaAssembler::new called with non-sparse kinship"
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

        // Build per-entry contribution table.
        let mut contributions: Vec<EntryContrib> = (0..nnz)
            .map(|_| EntryContrib {
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

        // Kinship contributions.
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
                    // Locate this row inside the union pattern of column j.
                    // Both slices are sorted ascending, so binary search.
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

        // Symbolic Cholesky factor of the union pattern. Build a temporary
        // owned `SymbolicSparseColMat` from clones of the pattern arrays so
        // the symbolic factor can be reused; faer Arc-clones internally.
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

    /// Compute the value vector for Σ at the given τ. Length matches the
    /// union pattern's nnz; positions correspond to `col_ptr` / `row_idx`.
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
        // groups passed for shape parity; the row_to_group lookup is the
        // bit we actually need.
        let _ = groups;
        values
    }

    /// Factor Σ at the given τ via the cached symbolic. Returns a fresh
    /// `Llt` factor.
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
        // Borrow the freshly built owned matrix as a Ref.
        let pattern_ref = symbolic_pattern.as_ref();
        let mat_ref = SparseColMatRef::<'_, u32, f64>::new(pattern_ref, &values);
        Llt::<u32, f64>::try_new_with_symbolic(self.symbolic.clone(), mat_ref, Side::Lower)
            .map_err(|e| {
                FavorError::Analysis(format!(
                    "sparse Σ Cholesky failed at current τ (matrix not positive definite): {e:?}"
                ))
            })
    }
}

// ─── Sparse matvec helper ───────────────────────────────────────────────────

/// Compute `K · v` for sparse symmetric K (only lower or full storage). The
/// kinship matrices we hold are full symmetric storage so we walk every
/// nonzero column slot.
fn sparse_matvec(k: &SparseColMat<u32, f64>, v: &Mat<f64>) -> Mat<f64> {
    debug_assert_eq!(v.ncols(), 1);
    let n = k.nrows();
    let mut out = Mat::<f64>::zeros(n, 1);
    let cp = k.symbolic().col_ptr();
    let ri = k.symbolic().row_idx();
    let val = k.val();
    for j in 0..n {
        let v_j = v[(j, 0)];
        let s = cp[j] as usize;
        let e = cp[j + 1] as usize;
        for kk in s..e {
            let row = ri[kk] as usize;
            out[(row, 0)] += val[kk] * v_j;
        }
    }
    out
}

// ─── Hutchinson stochastic trace estimator ──────────────────────────────────

/// Deterministic xorshift64 PRNG. Same `state` → same sequence. Used to
/// seed the Rademacher probes so two runs with identical inputs produce
/// identical Hutchinson estimates and identical p-values.
#[inline]
fn xorshift64_next(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

/// Stochastic estimate of `tr(Σ⁻¹ K_l)`:
///
/// ```text
/// tr(Σ⁻¹ K_l) ≈ (1/M) Σ_{m=1..M} z_m' · solve(Σ, K_l · z_m)
/// ```
///
/// where `z_m ∈ {-1, +1}^n` are independent Rademacher vectors. Each probe
/// costs one sparse matvec + one sparse Cholesky solve. With M=30 the
/// estimator has ~3% relative error and is fully reproducible given a
/// fixed seed.
pub fn hutchinson_trace(
    factor: &Llt<u32, f64>,
    k_l: &SparseColMat<u32, f64>,
    n_probes: usize,
    seed: u64,
) -> f64 {
    let n = k_l.nrows();
    let mut rng = seed.wrapping_mul(0x9E3779B97F4A7C15);
    if rng == 0 {
        rng = 1;
    }
    let mut accum = 0.0_f64;
    for _ in 0..n_probes {
        let mut z = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            let bit = xorshift64_next(&mut rng) & 1;
            z[(i, 0)] = if bit == 0 { -1.0 } else { 1.0 };
        }
        let kz = sparse_matvec(k_l, &z);
        let mut sigma_inv_kz = kz;
        factor.solve_in_place(&mut sigma_inv_kz);
        let mut probe = 0.0;
        for i in 0..n {
            probe += z[(i, 0)] * sigma_inv_kz[(i, 0)];
        }
        accum += probe;
    }
    accum / n_probes as f64
}

/// Stochastic estimate of the diagonal of Σ⁻¹ via Hutchinson. Used by the
/// score gradient for the group trace term `Σ_{i∈g} Σ⁻¹[i,i] / W[i]`.
///
/// ```text
/// diag(Σ⁻¹)[i] ≈ (1/M) Σ_m z_m[i] · solve(Σ, z_m)[i]
/// ```
///
/// Returns a length-`n` vector. Same seed → same output.
pub fn hutchinson_diag_inverse(
    factor: &Llt<u32, f64>,
    n: usize,
    n_probes: usize,
    seed: u64,
) -> Vec<f64> {
    let mut rng = seed.wrapping_mul(0xBF58476D1CE4E5B9);
    if rng == 0 {
        rng = 1;
    }
    let mut diag = vec![0.0_f64; n];
    for _ in 0..n_probes {
        let mut z = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            let bit = xorshift64_next(&mut rng) & 1;
            z[(i, 0)] = if bit == 0 { -1.0 } else { 1.0 };
        }
        let mut sigma_inv_z = z.clone();
        factor.solve_in_place(&mut sigma_inv_z);
        for i in 0..n {
            diag[i] += z[(i, 0)] * sigma_inv_z[(i, 0)];
        }
    }
    let inv_m = 1.0 / n_probes as f64;
    for v in &mut diag {
        *v *= inv_m;
    }
    diag
}

// ─── Sparse AI step ─────────────────────────────────────────────────────────

/// Per-iteration scratch from `ai_step_sparse`. Mirrors the dense `AiStep`.
struct SparseAiStep {
    d_tau: Vec<f64>,
    alpha: Mat<f64>,
    p_y: Mat<f64>,
    sigma_inv_x: Mat<f64>,
    cov: Mat<f64>,
    factor: Llt<u32, f64>,
    /// Diagonal of P at the current τ. Needed for the score's group trace term.
    diag_p: Vec<f64>,
}

/// Apply `P` to a column vector via the sparse factor:
/// `P v = Σ⁻¹ v − Σ⁻¹ X (X' Σ⁻¹ X)⁻¹ X' Σ⁻¹ v`.
fn apply_p(
    factor: &Llt<u32, f64>,
    sigma_inv_x: &Mat<f64>,
    cov: &Mat<f64>,
    x: &Mat<f64>,
    v: &Mat<f64>,
) -> Mat<f64> {
    let mut sinv_v = v.clone();
    factor.solve_in_place(&mut sinv_v);
    let xt_sinv_v = x.transpose() * &sinv_v;
    let cov_xt = cov * &xt_sinv_v;
    let proj = sigma_inv_x * &cov_xt;
    let n = v.nrows();
    let mut out = Mat::<f64>::zeros(n, 1);
    for i in 0..n {
        out[(i, 0)] = sinv_v[(i, 0)] - proj[(i, 0)];
    }
    out
}

/// One AI-REML iteration on the sparse path.
///
/// Mirrors `dense::ai_step` algorithmically — same score and AI matrix
/// formulae from `R/glmmkin.R::R_fitglmm_ai_var` — but never materializes
/// `Σ⁻¹`. The trace term `tr(Σ⁻¹ K_l)` and the group `diag(Σ⁻¹)[i]` term
/// come from the Hutchinson estimator (selected inversion is a follow-up).
fn ai_step_sparse(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    row_to_group: &[usize],
    tau: &VarianceComponents,
    fixtau: &[bool],
    assembler: &SparseSigmaAssembler,
    n_probes: usize,
) -> Result<SparseAiStep, FavorError> {
    let n = y.nrows();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    let k = x.ncols();
    debug_assert_eq!(fixtau.len(), n_comp);

    // Factor Σ at the current τ.
    let factor = assembler.factor(tau, groups, weights, row_to_group)?;

    // Σ⁻¹ X via column-wise solves on a fresh copy of X.
    let mut sigma_inv_x = x.clone();
    factor.solve_in_place(&mut sigma_inv_x);

    // (X' Σ⁻¹ X)⁻¹.
    let xt_sinv_x = x.transpose() * &sigma_inv_x;
    let cov = xt_sinv_x.col_piv_qr().solve(Mat::<f64>::identity(k, k));

    // Fixed effects: α̂ = cov · X' · Σ⁻¹ · y.
    let mut sinv_y = y.clone();
    factor.solve_in_place(&mut sinv_y);
    let xt_sinv_y = x.transpose() * &sinv_y;
    let alpha = &cov * &xt_sinv_y;

    // Projected response PY = Σ⁻¹ (Y - Xα).
    let eta = x * &alpha;
    let mut residual = Mat::<f64>::zeros(n, 1);
    for i in 0..n {
        residual[(i, 0)] = y[(i, 0)] - eta[(i, 0)];
    }
    let mut p_y = residual.clone();
    factor.solve_in_place(&mut p_y);

    // diag(P) — Hutchinson estimate for diag(Σ⁻¹) minus the X-projection
    // contribution which is computable analytically from sigma_inv_x and cov.
    let diag_sigma_inv = hutchinson_diag_inverse(&factor, n, n_probes, HUTCHINSON_SEED);
    let mut diag_p = vec![0.0_f64; n];
    for i in 0..n {
        let mut hi = 0.0;
        for a in 0..k {
            for b in 0..k {
                hi += sigma_inv_x[(i, a)] * cov[(a, b)] * sigma_inv_x[(i, b)];
            }
        }
        diag_p[i] = diag_sigma_inv[i] - hi;
    }

    // Score:
    //   kinship-l: U_l = (PY)' K_l (PY) − tr(P K_l)
    //   group-g:   U_g = Σ_{i∈g} (PY[i])²/W[i] − Σ_{i∈g} diag_p[i]/W[i]
    let mut score = vec![0.0_f64; n_comp];
    let mut k_py: Vec<Mat<f64>> = Vec::with_capacity(l);
    for li in 0..l {
        let k_sp = kinships[li].as_sparse().unwrap();
        let v = sparse_matvec(k_sp, &p_y);

        // tr(P K_l) = tr(Σ⁻¹ K_l) − tr(Σ⁻¹X cov X'Σ⁻¹ K_l).
        // First term: Hutchinson with the same fixed seed (different
        // mixing constant so it decorrelates from the diag estimator).
        let tr_sinv_k = hutchinson_trace(
            &factor,
            k_sp,
            n_probes,
            HUTCHINSON_SEED ^ ((li as u64).wrapping_mul(0xC2B2AE3D27D4EB4F)),
        );
        // Second term: tr(cov · X' Σ⁻¹ K_l Σ⁻¹ X). Define B = Σ⁻¹X (n×k).
        // X' Σ⁻¹ K_l Σ⁻¹ X = B' K_l B is k×k. Compute explicitly.
        let k_b = {
            let mut out = Mat::<f64>::zeros(n, k);
            for c in 0..k {
                let col = sigma_inv_x.col(c).to_owned();
                let mut col_mat = Mat::<f64>::zeros(n, 1);
                for i in 0..n {
                    col_mat[(i, 0)] = col[i];
                }
                let kb_col = sparse_matvec(k_sp, &col_mat);
                for i in 0..n {
                    out[(i, c)] = kb_col[(i, 0)];
                }
            }
            out
        };
        let bt_k_b = sigma_inv_x.transpose() * &k_b;
        let mut tr_correction = 0.0;
        for a in 0..k {
            for b in 0..k {
                tr_correction += cov[(a, b)] * bt_k_b[(b, a)];
            }
        }
        let tr_p_k = tr_sinv_k - tr_correction;

        let mut py_k_py = 0.0;
        for i in 0..n {
            py_k_py += p_y[(i, 0)] * v[(i, 0)];
        }
        score[li] = py_k_py - tr_p_k;
        k_py.push(v);
    }

    for gi in 0..g {
        let mut s_data = 0.0;
        let mut s_trace = 0.0;
        for &row in groups.group(gi) {
            let i = row as usize;
            let w_i = weights[i];
            s_data += p_y[(i, 0)] * p_y[(i, 0)] / w_i;
            s_trace += diag_p[i] / w_i;
        }
        score[l + gi] = s_data - s_trace;
    }

    // AI matrix: AI[i,j] = (PY)' V_i P V_j (PY). Same form as dense, only
    // P · v changes from a dense matvec to `apply_p`.
    let mut ai = Mat::<f64>::zeros(n_comp, n_comp);

    let mut p_k_py: Vec<Mat<f64>> = Vec::with_capacity(l);
    for li in 0..l {
        p_k_py.push(apply_p(&factor, &sigma_inv_x, &cov, x, &k_py[li]));
    }

    let mut v_py: Vec<Mat<f64>> = Vec::with_capacity(g);
    let mut p_v_py: Vec<Mat<f64>> = Vec::with_capacity(g);
    for gi in 0..g {
        let mut vpy = Mat::<f64>::zeros(n, 1);
        for &row in groups.group(gi) {
            let i = row as usize;
            vpy[(i, 0)] = p_y[(i, 0)] / weights[i];
        }
        let p_vpy = apply_p(&factor, &sigma_inv_x, &cov, x, &vpy);
        v_py.push(vpy);
        p_v_py.push(p_vpy);
    }

    // (l, l').
    for li in 0..l {
        for lj in 0..l {
            let mut s = 0.0;
            for i in 0..n {
                s += k_py[li][(i, 0)] * p_k_py[lj][(i, 0)];
            }
            ai[(li, lj)] = s;
        }
    }
    // (l, g).
    for li in 0..l {
        for gj in 0..g {
            let mut s = 0.0;
            for i in 0..n {
                s += k_py[li][(i, 0)] * p_v_py[gj][(i, 0)];
            }
            ai[(li, l + gj)] = s;
            ai[(l + gj, li)] = s;
        }
    }
    // (g, g').
    for gi in 0..g {
        for gj in 0..g {
            let mut s = 0.0;
            for &row in groups.group(gi) {
                let i = row as usize;
                s += v_py[gi][(i, 0)] * p_v_py[gj][(i, 0)];
            }
            ai[(l + gi, l + gj)] = s;
        }
    }

    // Solve AI Δτ = score over the free components only.
    let mut d_tau = vec![0.0_f64; n_comp];
    let n_free = fixtau.iter().filter(|&&f| !f).count();
    if n_free > 0 {
        let mut ai_free = Mat::<f64>::zeros(n_free, n_free);
        let mut score_free = Mat::<f64>::zeros(n_free, 1);
        let idx_map: Vec<usize> = (0..n_comp).filter(|&i| !fixtau[i]).collect();
        for (a, &fi) in idx_map.iter().enumerate() {
            score_free[(a, 0)] = score[fi];
            for (b, &fj) in idx_map.iter().enumerate() {
                ai_free[(a, b)] = ai[(fi, fj)];
            }
        }
        let dt_free = ai_free.col_piv_qr().solve(score_free);
        for (a, &fi) in idx_map.iter().enumerate() {
            d_tau[fi] = dt_free[(a, 0)];
        }
    }

    Ok(SparseAiStep {
        d_tau,
        alpha,
        p_y,
        sigma_inv_x,
        cov,
        factor,
        diag_p,
    })
}

// ─── Sparse fit_reml outer loop ─────────────────────────────────────────────

fn fit_reml_inner_sparse(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: &[f64],
    row_to_group: &[usize],
    init_tau: VarianceComponents,
    fixtau: &[bool],
    assembler: &SparseSigmaAssembler,
    n_probes: usize,
) -> Result<KinshipState, FavorError> {
    let n = y.nrows();
    let k = x.ncols();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;

    let mut tau = init_tau;
    let mut alpha_prev = Mat::<f64>::zeros(k, 1);
    let mut iter_used = 0;
    let mut last: Option<SparseAiStep> = None;

    for iter in 0..REML_MAX_ITER {
        iter_used = iter + 1;
        let tau0 = tau.clone();
        let step_out = ai_step_sparse(
            y,
            x,
            kinships,
            groups,
            weights,
            row_to_group,
            &tau,
            fixtau,
            assembler,
            n_probes,
        )?;

        // Newton step + step-halving on negative τ.
        let mut step = step_out.d_tau.clone();
        let mut try_count = 0;
        loop {
            for i in 0..n_comp {
                if !fixtau[i] {
                    tau.as_slice_mut()[i] = tau0.as_slice()[i] + step[i];
                }
            }
            for i in 0..n_comp {
                let v = tau.as_slice()[i];
                if v < REML_TOL && tau0.as_slice()[i] < REML_TOL {
                    tau.as_slice_mut()[i] = 0.0;
                }
            }
            if tau.as_slice().iter().all(|&t| t >= 0.0) {
                break;
            }
            for s in step.iter_mut() {
                *s /= 2.0;
            }
            try_count += 1;
            if try_count > REML_STEP_HALVE_MAX {
                return Err(FavorError::Analysis(format!(
                    "sparse AI-REML step-halving failed in {REML_STEP_HALVE_MAX} attempts"
                )));
            }
        }
        for v in tau.as_slice_mut().iter_mut() {
            if *v < REML_TOL {
                *v = 0.0;
            }
        }

        let mut rel_alpha = 0.0_f64;
        for i in 0..k {
            let denom = step_out.alpha[(i, 0)].abs() + alpha_prev[(i, 0)].abs() + REML_TOL;
            let r = (step_out.alpha[(i, 0)] - alpha_prev[(i, 0)]).abs() / denom;
            if r > rel_alpha {
                rel_alpha = r;
            }
        }
        let mut rel_tau = 0.0_f64;
        for i in 0..n_comp {
            let denom = tau.as_slice()[i].abs() + tau0.as_slice()[i].abs() + REML_TOL;
            let r = (tau.as_slice()[i] - tau0.as_slice()[i]).abs() / denom;
            if r > rel_tau {
                rel_tau = r;
            }
        }

        alpha_prev = step_out.alpha.clone();
        last = Some(step_out);

        if 2.0 * rel_alpha.max(rel_tau) < REML_TOL {
            break;
        }

        if tau.as_slice().iter().any(|&t| t.abs() > REML_TOL.powi(-2)) {
            return Err(FavorError::Analysis(format!(
                "sparse AI-REML diverged: |τ_max| > {}",
                REML_TOL.powi(-2)
            )));
        }
    }

    let last = last
        .ok_or_else(|| FavorError::Analysis("sparse AI-REML produced no iterations".into()))?;

    let total_var: f64 = tau.as_slice().iter().sum();
    let h2: Vec<f64> = if total_var > 0.0 {
        (0..l).map(|i| tau.kinship(i) / total_var).collect()
    } else {
        vec![0.0; l]
    };

    let _ = n;
    let _ = last.diag_p;

    Ok(KinshipState {
        tau,
        inverse: KinshipInverse::Sparse(SparseFactor { llt: last.factor }),
        sigma_inv_x: last.sigma_inv_x,
        cov: last.cov,
        p_y: last.p_y,
        h2,
        n_iter: iter_used,
        outer_refits: 0,
    })
}

/// Outer boundary-refit loop on the sparse path. Caller is the public
/// `fit_reml` dispatcher in `mod.rs`.
pub fn fit_reml_sparse(
    y: &Mat<f64>,
    x: &Mat<f64>,
    kinships: &[KinshipMatrix],
    groups: &GroupPartition,
    weights: Option<&[f64]>,
    init_tau: VarianceComponents,
    n_probes: usize,
) -> Result<KinshipState, FavorError> {
    let n = y.nrows();
    let l = kinships.len();
    let g = groups.n_groups();
    let n_comp = l + g;
    let w_vec = weights_or_one(n, weights);

    // row → group lookup, used to assemble the diagonal residual term.
    let mut row_to_group = vec![0_usize; n];
    for gi in 0..groups.n_groups() {
        for &row in groups.group(gi) {
            row_to_group[row as usize] = gi;
        }
    }

    let assembler = SparseSigmaAssembler::new(n, kinships)?;

    let mut fixtau = vec![false; n_comp];
    let mut state = fit_reml_inner_sparse(
        y,
        x,
        kinships,
        groups,
        &w_vec,
        &row_to_group,
        init_tau,
        &fixtau,
        &assembler,
        n_probes,
    )?;

    for refit_iter in 0..REML_MAX_OUTER_REFITS {
        let mut changed = false;
        for i in 0..n_comp {
            if !fixtau[i] && state.tau.as_slice()[i] < BOUNDARY_FACTOR * REML_TOL {
                fixtau[i] = true;
                changed = true;
            }
        }
        if !changed {
            state.outer_refits = refit_iter;
            return Ok(state);
        }
        let mut tau_refit = state.tau.clone();
        for i in 0..n_comp {
            if fixtau[i] {
                tau_refit.as_slice_mut()[i] = 0.0;
            }
        }
        state = fit_reml_inner_sparse(
            y,
            x,
            kinships,
            groups,
            &w_vec,
            &row_to_group,
            tau_refit,
            &fixtau,
            &assembler,
            n_probes,
        )?;
        state.outer_refits = refit_iter + 1;
    }

    Err(FavorError::Analysis(format!(
        "sparse AI-REML boundary refit did not stabilize within {REML_MAX_OUTER_REFITS} iterations"
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::staar::kinship::dense::fit_reml_dense;
    use faer::sparse::Triplet;

    fn pedigree_kinship(n_fam: usize, sibs: usize) -> KinshipMatrix {
        let n = n_fam * sibs;
        let mut entries: Vec<Triplet<u32, u32, f64>> = Vec::new();
        for f in 0..n_fam {
            let base = f * sibs;
            for i in 0..sibs {
                for j in 0..sibs {
                    let v = if i == j { 1.0 } else { 0.5 };
                    entries.push(Triplet::new((base + i) as u32, (base + j) as u32, v));
                }
            }
        }
        KinshipMatrix::from_triplets(n, entries, "ped".into()).unwrap()
    }

    fn dense_pedigree_kinship(n_fam: usize, sibs: usize) -> KinshipMatrix {
        let n = n_fam * sibs;
        let mut k = Mat::<f64>::zeros(n, n);
        for f in 0..n_fam {
            let base = f * sibs;
            for i in 0..sibs {
                for j in 0..sibs {
                    k[(base + i, base + j)] = if i == j { 1.0 } else { 0.5 };
                }
            }
        }
        // KinshipMatrix::new will route this back to sparse since density is
        // low; force the dense variant for the comparison test.
        KinshipMatrix::Dense {
            matrix: k,
            label: "ped-dense".into(),
        }
    }

    fn synthetic_pheno(n: usize) -> Mat<f64> {
        let mut y = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            let fam = i / 5;
            let g = ((fam as f64 * 0.91).sin()) * 0.6;
            let e = (i as f64 * 1.7).sin() * 0.4;
            y[(i, 0)] = g + e;
        }
        y
    }

    fn intercept(n: usize) -> Mat<f64> {
        let mut x = Mat::<f64>::zeros(n, 1);
        for i in 0..n {
            x[(i, 0)] = 1.0;
        }
        x
    }

    #[test]
    fn sparse_assembler_pattern_includes_diagonal() {
        let kin = pedigree_kinship(20, 5);
        let n = kin.n();
        let assembler = SparseSigmaAssembler::new(n, std::slice::from_ref(&kin)).unwrap();
        // Diagonal must be present in every column.
        for j in 0..n {
            let s = assembler.col_ptr[j] as usize;
            let e = assembler.col_ptr[j + 1] as usize;
            let mut has_diag = false;
            for k in s..e {
                if assembler.row_idx[k] as usize == j {
                    has_diag = true;
                    break;
                }
            }
            assert!(has_diag, "missing diagonal at column {j}");
        }
    }

    #[test]
    fn sparse_factor_succeeds_at_warm_start() {
        let kin = pedigree_kinship(20, 5);
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let assembler = SparseSigmaAssembler::new(n, std::slice::from_ref(&kin)).unwrap();
        let mut tau = VarianceComponents::zeros(1, 1);
        tau.set_kinship(0, 0.5);
        tau.set_group(0, 0.5);
        let weights = vec![1.0; n];
        let row_to_group = vec![0_usize; n];
        let factor = assembler.factor(&tau, &groups, &weights, &row_to_group);
        assert!(factor.is_ok(), "Σ at warm-start τ should be SPD");
    }

    #[test]
    fn hutchinson_trace_reproducible() {
        let kin = pedigree_kinship(40, 5);
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let weights = vec![1.0; n];
        let row_to_group = vec![0_usize; n];
        let assembler = SparseSigmaAssembler::new(n, std::slice::from_ref(&kin)).unwrap();
        let mut tau = VarianceComponents::zeros(1, 1);
        tau.set_kinship(0, 0.4);
        tau.set_group(0, 0.6);
        let factor = assembler
            .factor(&tau, &groups, &weights, &row_to_group)
            .unwrap();
        let k_sp = kin.as_sparse().unwrap();
        let t1 = hutchinson_trace(&factor, k_sp, 30, 12345);
        let t2 = hutchinson_trace(&factor, k_sp, 30, 12345);
        assert_eq!(t1, t2, "same seed must produce identical estimate");
    }

    #[test]
    fn sparse_fit_reml_pedigree_matches_dense_within_tolerance() {
        // 20 families × 5 sibs = 100 samples. Sparse path with M=200 probes
        // should land within ~5% of the dense REML τ values.
        let kin_sparse = pedigree_kinship(20, 5);
        let kin_dense = dense_pedigree_kinship(20, 5);
        let n = kin_sparse.n();
        let groups = GroupPartition::single(n);
        let y = synthetic_pheno(n);
        let x = intercept(n);

        let weights = vec![1.0_f64; n];

        // Warm start identical to mod.rs::fit_reml.
        let y_mean = (0..n).map(|i| y[(i, 0)]).sum::<f64>() / n as f64;
        let y_var = (0..n)
            .map(|i| (y[(i, 0)] - y_mean).powi(2))
            .sum::<f64>()
            / (n as f64 - 1.0).max(1.0);
        let mut warm_init = VarianceComponents::zeros(1, 1);
        let warm = y_var / 2.0;
        warm_init.set_kinship(0, warm / kin_sparse.mean_diagonal());
        warm_init.set_group(0, warm);

        let dense_state = fit_reml_dense(
            &y,
            &x,
            std::slice::from_ref(&kin_dense),
            &groups,
            &weights,
            warm_init.clone(),
        )
        .expect("dense fit");

        // High probe count keeps Hutchinson noise small enough that the
        // tolerance below holds reliably.
        let sparse_state = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin_sparse),
            &groups,
            None,
            warm_init,
            200,
        )
        .expect("sparse fit");

        let dense_tau_k = dense_state.tau.kinship(0);
        let sparse_tau_k = sparse_state.tau.kinship(0);
        let dense_tau_g = dense_state.tau.group(0);
        let sparse_tau_g = sparse_state.tau.group(0);

        let scale = (dense_tau_k.abs() + dense_tau_g.abs()).max(1e-3);
        let err_k = (dense_tau_k - sparse_tau_k).abs() / scale;
        let err_g = (dense_tau_g - sparse_tau_g).abs() / scale;

        // 15% relative tolerance is generous for M=200 probes on n=100;
        // tightening it requires more probes or selected inversion.
        assert!(
            err_k < 0.15,
            "τ_kinship: dense={dense_tau_k:.4} sparse={sparse_tau_k:.4} rel_err={err_k:.3}"
        );
        assert!(
            err_g < 0.15,
            "τ_group:   dense={dense_tau_g:.4} sparse={sparse_tau_g:.4} rel_err={err_g:.3}"
        );

        // Sparse path produced a sparse factor (not a dense Σ⁻¹).
        assert!(matches!(sparse_state.inverse, KinshipInverse::Sparse(_)));
    }

    #[test]
    fn sparse_fit_reml_reproducible() {
        let kin = pedigree_kinship(15, 4);
        let n = kin.n();
        let groups = GroupPartition::single(n);
        let y = synthetic_pheno(n);
        let x = intercept(n);

        let mut warm_init = VarianceComponents::zeros(1, 1);
        warm_init.set_kinship(0, 0.5);
        warm_init.set_group(0, 0.5);

        let s1 = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin),
            &groups,
            None,
            warm_init.clone(),
            64,
        )
        .expect("first sparse fit");
        let s2 = fit_reml_sparse(
            &y,
            &x,
            std::slice::from_ref(&kin),
            &groups,
            None,
            warm_init,
            64,
        )
        .expect("second sparse fit");

        assert_eq!(s1.tau.kinship(0), s2.tau.kinship(0));
        assert_eq!(s1.tau.group(0), s2.tau.group(0));
        assert_eq!(s1.n_iter, s2.n_iter);
    }
}
