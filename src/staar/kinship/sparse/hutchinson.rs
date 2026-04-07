//! Stochastic trace and diagonal estimators for the sparse path.
//!
//! Both estimators are Hutchinson-style: pick `M` Rademacher probes, apply
//! Σ⁻¹ via the cached Cholesky factor, and average. Reproducible given a
//! fixed seed — same seed in, same float out, run-to-run.
//!
//! These are *approximations*. Standard error scales as `1/√M`; with
//! `M = 30` the relative error is around 3%. The default sparse path uses
//! the exact Takahashi recursion in [`super::takahashi`] instead. This
//! module is kept as a fallback under `SparseSolverKind::Hutchinson` and
//! is used by the regression test that confirms the stochastic path still
//! produces valid-within-tolerance estimates against the dense path.

use faer::sparse::linalg::solvers::Llt;
use faer::sparse::SparseColMat;
use faer::Mat;

/// Default Hutchinson probe count. `M = 30` gives ~3% relative error and
/// converges quickly enough that the AI loop still hits its outer
/// convergence criterion in a reasonable number of iterations.
pub const DEFAULT_HUTCHINSON_PROBES: usize = 30;

/// Fixed seed for the Rademacher PRNG. Same seed → same probes → same
/// estimates → run-to-run reproducibility for both p-values and
/// debugging.
pub const HUTCHINSON_SEED: u64 = 0x5F3759DF_5F3759DF;

/// Deterministic xorshift64. Equivalent stream to a small Mersenne or
/// PCG would work; xorshift is the cheapest deterministic choice.
#[inline]
pub fn xorshift64_next(state: &mut u64) -> u64 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *state = x;
    x
}

/// Sparse `K · v` for full-storage symmetric kinship. Walks every
/// nonzero column slot — the kinship matrices we hold are full
/// symmetric storage, not lower-triangular.
pub fn sparse_matvec(k: &SparseColMat<u32, f64>, v: &Mat<f64>) -> Mat<f64> {
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

/// Stochastic estimate of `tr(Σ⁻¹ K_l)`:
///
/// ```text
/// tr(Σ⁻¹ K_l) ≈ (1/M) Σ_{m=1..M} z_m' · solve(Σ, K_l · z_m)
/// ```
///
/// where each `z_m ∈ {-1, +1}^n` is a Rademacher vector. Cost per probe
/// is one sparse matvec plus one sparse Cholesky solve. With `M = 30`
/// the relative error is ~3%.
///
/// **Approximation, not exact.** The default sparse path uses
/// [`super::takahashi`] instead.
pub fn trace_k_estimate(
    factor: &Llt<u32, f64>,
    k_l: &SparseColMat<u32, f64>,
    n_probes: usize,
    seed: u64,
) -> f64 {
    use faer::linalg::solvers::Solve;

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

/// Stochastic estimate of `diag(Σ⁻¹)`. Used by the score's group trace
/// term `Σ_{i∈g} Σ⁻¹[i,i] / W[i]`.
///
/// ```text
/// diag(Σ⁻¹)[i] ≈ (1/M) Σ_m z_m[i] · solve(Σ, z_m)[i]
/// ```
///
/// Returns a length-`n` vector. Same seed → same output. The default
/// sparse path uses the Takahashi selected inverse instead.
pub fn diag_inverse_estimate(
    factor: &Llt<u32, f64>,
    n: usize,
    n_probes: usize,
    seed: u64,
) -> Vec<f64> {
    use faer::linalg::solvers::Solve;

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
