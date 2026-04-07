# upstream pinning

This file pins the upstream R sources we port from, and records the rename
map so contributors can read our code side-by-side with the R. Update the
SHAs and re-walk the audit when bumping upstream versions.

## pinned commits

| repo | branch | sha | date |
|---|---|---|---|
| hanchenphd/GMMAT | master | `473b34268ade7f8b920f67657c6f4759880d35d3` | 2025-11-21 |
| UW-GAC/GENESIS   | devel  | `dd480380c7705d36d1dbeb09c73b08f60623701a` | 2026-03-04 |

GMMAT is the primary source — our dense and sparse AI-REML kernels are
ports of `R/glmmkin.R` and `src/fitglmm.cpp`. GENESIS is a cross-check;
its formulas in `R/runAIREMLgaussian.R` agree with GMMAT, and its
`R/iterateAIREMLworkingY.R` validates the two-loop structure we use for
the binary PQL path (outer loop updates the working response, inner loop
runs AI-REML).

## function map

The R-to-rust mapping. Citations in code headers should use this table.

| upstream function | upstream file:lines | our file | our function |
|---|---|---|---|
| `glmmkin` (top-level entry) | `R/glmmkin.R:1-120` | `mod.rs` | `fit_reml` (dispatcher), `fit_pql_glmm` |
| `glmmkin.fit` (family + boundary refit loop) | `R/glmmkin.R:122-279` | `dense.rs::fit_reml_dense` outer loop, `sparse.rs::fit_reml_sparse` outer loop | both — they each have a copy of the boundary refit |
| `glmmkin.ai` (AI-REML iteration loop, interleaves working-response update) | `R/glmmkin.R:281-422` | `dense.rs::fit_reml_inner`, `sparse.rs::fit_reml_inner_sparse` | iteration loop only — we factor the working-response update out into `pql.rs` for binary, and Gaussian skips it |
| `R_fitglmm_ai` (one AI step, sparse-aware) | `R/glmmkin.R:662-710` | `sparse.rs::ai_step_sparse` | one AI step |
| `R_fitglmm_ai_dense` (one AI step, pure-R dense) | `R/glmmkin.R:712-769` | `dense.rs::ai_step` | one AI step |
| `fitglmm_ai` (one AI step, C++ dense, used when n < 2^15.5) | `src/fitglmm.cpp:541-` | `dense.rs::ai_step` | same math as above; C++ is just a faster path for the same formulas |
| `glmmkin.fit` PQL outer loop (Gaussian no-op, binary working-response update interleaved with AI step) | `R/glmmkin.R:122-279` + `R/glmmkin.R:401-405` | `pql.rs::fit_pql_glmm` | two-loop factoring; see PQL note below |
| GENESIS `.iterateAIREMLworkingY` (separate PQL outer loop) | `R/iterateAIREMLworkingY.R:17-63` | `pql.rs::fit_pql_glmm` | structural cross-check — GENESIS factors PQL the same way we do |
| GENESIS `.runAIREMLgaussian` (AI-REML inner loop) | `R/runAIREMLgaussian.R:1-143` | `dense.rs::fit_reml_inner` | cross-check on score, AI matrix, convergence |

## constant provenance

The four constants in `dense.rs`:

| constant | value | upstream? | notes |
|---|---|---|---|
| `REML_TOL` | `1e-5` | yes — `glmmkin.R:122` `tol = 1e-5` | identical |
| `REML_MAX_ITER` | `500` | yes — `glmmkin.R:122` `maxiter = 500` | identical |
| `BOUNDARY_FACTOR` | `1.01` | yes — appears throughout `glmmkin.R` as the literal `1.01 * tol` (lines 125, 187, 191, 203, 207, 245, 253, 260, 267, 274, 369, 375, 385, 486, 492, 502) | upstream uses it as the threshold `tau < 1.01 * tol` to decide a component has hit the boundary |
| `REML_MAX_OUTER_REFITS` | `10` | **no** — upstream's boundary refit `while` loop at `glmmkin.R:194-210` is unbounded | our cap is a safety net; on hit we should return an error pointing the user at the diagnostic |
| `REML_STEP_HALVE_MAX` | `50` | **no** — upstream's step-halving loop at `glmmkin.R:362-366` is unbounded | same — our cap; error on hit |

For PR A we keep the safety caps but document them as ours, not upstream's.
For PR C (test gaps) we add error-path tests for both caps being hit.

## convergence criterion

GMMAT (`glmmkin.R:405`):

```r
2 * max(abs(alpha - alpha0)/(abs(alpha) + abs(alpha0) + tol),
        abs(tau   - tau0)  /(abs(tau)   + abs(tau0)   + tol)) < tol
```

This is what `dense.rs::fit_reml_inner` and `sparse.rs::fit_reml_inner_sparse`
implement. Identical formulation. Cite this line at the convergence check.

GENESIS uses a different criterion (`abs(sigma2.kplus1 - sigma2.k) < tol*var(Y)` —
absolute tolerance scaled by trait variance, no alpha term). We follow GMMAT.
Both are valid; the choice is upstream-author preference.

## warm-start

GMMAT (`glmmkin.R:309-316`):

```r
tau[idxtau] <- rep(var(Y) / (q + ng), q2)        # split variance evenly
tau[i+ng]   <- tau[i+ng] / mean(diag(kins[[i]])) # scale kinship slots
```

`mod.rs:104-119` matches this. Cite at the warm-start.

## sparse path is not bit-identical to upstream

This is the audit's most important finding. Upstream's sparse path
(`R_fitglmm_ai`, `glmmkin.R:662-710`) does the trace term `tr(Sigma_i K_l)`
exactly via the Frobenius inner product `sum(Sigma_i * kins)`. To do that
it computes `Sigma_i = chol2inv(chol(Sigma))` — which materializes the full
inverse as a dense `n x n` matrix even when the Cholesky factor is sparse,
because the inverse of a sparse SPD matrix is generally dense.

So upstream's "sparse" path is really "sparse Cholesky factorization, dense
inverse, exploit kinship sparsity in the Frobenius product". It costs `O(n^2)`
memory.

Our sparse path (`sparse.rs::ai_step_sparse`) does NOT materialize `Sigma_i`.
It uses `factor.solve_in_place()` for `Sigma_i v` and a Hutchinson stochastic
estimator for `tr(Sigma_i K_l)`. This costs `O(nnz_L * m)` memory where m is
the probe count (default 30). It is an approximation with ~3% relative error
on the trace term, which propagates into tau via the AI step.

The Takahashi recursion in Phase 4 (#26 #27) replaces Hutchinson with an
exact computation of `Sigma_i` entries restricted to the union sparsity pattern
of all kinship matrices. Once that lands the sparse path is mathematically
1:1 with upstream `R_fitglmm_ai` while still using strictly less memory than
upstream (we never materialize the full `Sigma_i`, only the entries the
Frobenius product needs).

Until Phase 4 lands, the sparse path is "good approximation, ~3% trace
error". After Phase 4 it is "1:1 with upstream, exact at the entry pattern".
Document this clearly in PR A.

## PQL note

GMMAT does NOT have a separate PQL outer loop. Working-response update is
interleaved with the AI step inside `glmmkin.ai` (`glmmkin.R:401-405`):

```r
mu      <- family$linkinv(eta)
mu.eta  <- family$mu.eta(eta)
Y       <- eta - offset + (y - mu) / mu.eta
sqrtW   <- mu.eta / sqrt(...)
```

These four lines run between AI iterations within the same convergence loop.
For Gaussian family, `linkinv` is identity and `mu.eta` is 1, so the update
is a no-op and `Y = y` throughout.

GENESIS factors this differently (`R/iterateAIREMLworkingY.R:17-63`): outer
loop updates the working response, inner loop runs the full AI-REML to
convergence, repeat until the working response stops changing. This is a
nested-loop structure.

Our `pql.rs::fit_pql_glmm` follows GENESIS, not GMMAT. The two-loop and
single-loop structures converge to the same point but follow different paths.
We cite both upstreams: GENESIS for the structure, GMMAT for the formulas.

## variable rename map

When reading our dense and sparse AI step against the upstream R, use this
table to translate names. Add it as a header comment to `reml/ai.rs` once
Phase 1 lands.

| upstream R    | our rust              | meaning |
|---|---|---|
| `Y`           | `y`                   | response (working response Y for binary, raw y for Gaussian) |
| `X`           | `x`                   | covariate matrix |
| `tau` / `theta` | `tau` (`VarianceComponents`) | variance components, layout `[group_1..group_g, kinship_1..kinship_q]` upstream — **note we use `[kinship_1..kinship_q, group_1..group_g]`, the reverse order** |
| `q`           | `kinships.len()`      | number of kinship matrices |
| `ng`          | `groups.len()`        | number of group components |
| `group.idx`   | `GroupPartition`      | sample-to-group assignment |
| `W`           | `weights` (or `1.0` for Gaussian) | IRLS weights, `mu.eta^2 / variance(mu)` |
| `kins[[i]]`   | `kinships[i]`         | i-th kinship matrix |
| `Sigma`       | local `sigma`         | covariance matrix Σ = diag(τ_g/W) + Σ τ_l K_l |
| `Sigma_i`     | `sigma_inv`           | Σ⁻¹ |
| `Sigma_iX`    | `sigma_inv_x`         | Σ⁻¹ X |
| `XSigma_iX`   | local `xt_sigma_inv_x`| X' Σ⁻¹ X |
| `cov`         | `cov`                 | (X' Σ⁻¹ X)⁻¹ |
| `Sigma_iXcov` | `sigma_inv_x_cov`     | Σ⁻¹ X (X' Σ⁻¹ X)⁻¹ |
| `P`           | (not materialized in sparse path) | projection P = Σ⁻¹ - Σ⁻¹ X (X'Σ⁻¹X)⁻¹ X' Σ⁻¹ |
| `PY`          | `p_y`                 | P y |
| `wPY`         | `w_p_y`               | PY / W (component-wise) |
| `diagP`       | `diag_p`              | diag(P) / W |
| `APY`         | local `a_py`          | K_l PY |
| `PAPY`        | local `p_a_py`        | P K_l PY |
| `score`       | `score`               | score vector U for free components |
| `AI`          | `ai`                  | average information matrix |
| `Dtau`        | `delta_tau`           | Newton step |
| `fixtau`      | `fixed_mask`          | which components are fixed at the boundary |
| `idxtau`      | `free_idx`            | indices of free components |
| `q2`          | `n_free`              | number of free components |

The `tau` ordering is the only thing that confuses readers. Upstream R puts
groups first (`tau[1..ng]` are groups, `tau[ng+1..]` are kinships). Our
`VarianceComponents` puts kinships first. We chose this because the kinship
slots are the conceptually primary parameters. When porting any formula
with `idxtau[i] <= ng` (group component check), translate to `idx >= L`
where L is `kinships.len()`.

## what to do with this file

- Phase 1 (restructure): each new file in `reml/`, `dense/`, `sparse/` gets
  a header citing the relevant entry from "function map" above. Use the
  exact line ranges.
- Phase 2 (audit): every math function gets a `// upstream:` line at the
  function header pointing to the upstream lines. The variable rename map
  goes inline as a comment block above `ai_step` so reviewers can read the
  rust against the R.
- Phase 4 (Takahashi): updates the "sparse path is not bit-identical"
  section to read "1:1 with upstream after Takahashi landed in PR B".
- Bumping upstream: refresh the SHAs at the top, re-walk the line ranges,
  re-confirm constants haven't moved.
