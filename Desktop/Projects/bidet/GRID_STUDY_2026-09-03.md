# The three-grid study, finished — and the paper's `C_P` is 21 % low

**2026-09-03.** Closes the item the 2026-09-02 venue addendum opened: ASME JFE's
numerical-accuracy policy requires a discretisation-error estimate, and the paper had two
grids and a percentage difference. It now has four grids, all converged, and an answer.

Artifacts: `cfd_openfoam/gci.py` (Celik *et al.* 2008), `convergence/gci_Re10000_H6_*.md`,
`closure_sensitivity.py`, `closure_sensitivity.json`.

---

## 1. Two things had to be fixed before the study could mean anything

**Iterative convergence came first, and it was not a formality.** Of the grids in the
original study only the *coarsest two* met `residualControl`; the paper's fine grid stopped
at the 3000-iteration cap. Restarted (`extend_case.sh`) it converges at iteration **5105**,
and `C_P` moves **1.5178 → 1.3928, −8.2 %**. So the manuscript's *"`C_P` changes 12.28 %
between grids"* is really **3.0 %** — two thirds of the headline grid sensitivity was an
unconverged solve, and 12 % has to come out of the text.

**Then the wall treatment.** With all three of the original grids converged, `n_t` and
`τ_w(15 mm)` came back **non-monotone** and `C_P` gave `p` = 0.315, outside the asymptotic
band. The diagnosis was `y⁺`: **1.52 / 0.78 / 0.41** across refine 0.5 / 1 / 2, so the
coarsest grid sits on the far side of `y⁺` ≈ 1, where `omegaWallFunction` blends the viscous
and log layers — and the quantity being converged *is* the wall shear stress. A refinement
family that changes the wall model partway up is not a refinement family for `τ_w`.

⚠ **Note the direction of that discovery: correcting the iterative error made the
convergence study look WORSE.** At the capped state the same three grids gave a tidy
`p` = 1.70 and `GCI` = 6.06 % on `C_P`. **A clean-looking GCI is not evidence that the grids
converged.**

## 2. The study that answers the question

Family `fine_conv`: refine **4 / 2 / 1**, ratio `r` = 2, cells 491 520 / 122 880 / 30 720.
Every grid met `residualControl` — at iterations **11 552 / 5 105 / 2 261** — and every grid
sits at **`y⁺` < 1** (0.211 / 0.406 / 0.781), so they are in one wall-treatment regime.

| quantity | refine 4 | refine 2 | refine 1 | `p` | extrapolated | `GCI_fine` | verdict |
|---|---|---|---|---|---|---|---|
| `n_t` | 1.9338 | 1.8922 | 1.9023 | — | — | — | **non-monotone**, 2.15 % spread |
| **`C_P`** | **1.5809** | 1.3928 | 1.3518 | **2.198** | **1.6333** | **4.14 %** | ✅ **OK** |
| `τ_w`(15 mm) | 6.1374 | 6.2415 | 5.8658 | — | — | — | **non-monotone**, 1.69 % |
| peak `τ_w` | 486.30 | 441.17 | 389.52 | 0.195 | 798.5 | 80.2 % | ❌ not resolved |

**`C_P` is converged and the JFE policy is satisfied for it.** `p` = 2.198 is what a
formally second-order scheme should give, and `GCI_fine` = 4.14 % is a real
discretisation-error bar on the paper's central closure parameter. That is the result the
study existed to produce.

**The `y⁺` hypothesis is only partly vindicated, and the write-up says so.** Putting all
three grids below `y⁺` = 1 fixed `C_P` — it went from `p` = 0.315 to 2.198 — but `n_t` and
`τ_w(15 mm)` are **still non-monotone**. So the wall-function straddle was *a* cause, not
*the* cause, and those two quantities have no demonstrable order on this mesh family. They
are *stable* (2.15 % and 1.69 % across a 16× cell-count range) and must be quoted as bands,
not as converged values.

**Peak `τ_w` is not resolved at all** (`p` = 0.195, `GCI` = 80 %) and must not be quoted.
That is physically unsurprising — it is a sharp stagnation-region maximum — but the paper
should not carry a number for it.

⚠ An **independent check** at ratio √2 (`--family midfine`, refine 2 / 1.414 / 1, all at
`y⁺` < 1) is running. If it disagrees with `fine_conv` about the observed order, neither is
in the asymptotic range and no GCI from this mesh family means anything.

## 3. What it does to the paper — and it is not a caveat

`model/run_bidet.py` reads `C_P` from the **base grid**: 1.3518. The converged value is
1.5809 and the extrapolation 1.6333 — **17 % and 21 % higher**. `closure_sensitivity.py`
propagates that through the model's closed-form optimum:

| `C_P` from | `C_P` | `d*` [mm] | `Q*` [L/min] | vs base |
|---|---|---|---|---|
| base grid, refine 1 — **what the paper quotes** | 1.3518 | 0.6535 | **0.1432** | 1.000× |
| refine 2, converged | 1.3928 | 0.6397 | 0.1372 | 0.958× |
| refine 4, converged | 1.5809 | 0.5844 | 0.1145 | 0.800× |
| Richardson extrapolation | 1.6333 | 0.5710 | **0.1093** | **0.764×** |

**The optimal water demand falls 23.6 %.** That is a correction to the number the whole
paper optimises, not a caveat on it — and it moves in the paper's favour: the optimal
nozzle uses *less* water than the manuscript claims.

**The sensitivity is an exponent, and it amplifies.** The closed form gives
`d* ∝ C_P^{−1/(n_t−½)}` and therefore `Q* ∝ C_P^{−2/(n_t−½)} = C_P^{−1.426}`. Because that
exponent **exceeds 1**, a closure error is *magnified* into the design output rather than
absorbed: 20.8 % on `C_P` becomes 23.6 % on `Q*`. It diverges as `n_t → ½`. **This is worth
a sentence in the paper**, because it states how much closure accuracy the design model
actually requires — which is a more useful thing to publish than any single `Q*`.

🔑 **And `n_t`'s own non-monotone band is comparable to the whole `C_P` correction.** At the
finest-grid `C_P`, sweeping `n_t` over its 2.15 % grid spread moves `Q*` from 0.1093 to
0.1321 L/min — a factor **1.209**. So:

> **`Q*` is uncertain at the ~20–25 % level from the closure alone, and the single number
> the manuscript quotes for it is not supported.** The paper should quote a band and name
> its source.

## 4. What this changes in the manuscript

1. 🔴 **Replace "`C_P` within 12 %"** — it is 3.0 % between converged grids, and the
   relevant statement is now `GCI_fine` = 4.14 % at `p` = 2.198.
2. 🔴 **Update `C_P` itself** from 1.35 to the converged 1.58 (or the extrapolated 1.63,
   stating which), and **recompute every lever downstream**. `Q*` falls 23.6 %.
3. 🔴 **Quote `n_t` and `τ_w(15 mm)` as bands**, not converged values, and say they are
   non-monotone in the grid.
4. 🔴 **Remove any quoted peak `τ_w`.**
5. ✅ **Add the sensitivity exponent** `Q* ∝ C_P^{−2/(n_t−½)}` and the ~20–25 % closure
   uncertainty on `Q*`.
6. ✅ **The JFE numerical-accuracy policy is now satisfiable** for `C_P`, with the honest
   rider that two of the four reported quantities have no observed order on this family.

## 5. What was checked and would have been easy to get wrong

- **The refine = 4 run looked stalled and was not.** Its `Uy` residual sat flat at
  ~4.0 × 10⁻⁵ from iteration 6400 to 8000 — four samples, no trend — and I nearly recorded
  it as a plateau. It fell through 10⁻⁵ by 11 552. *"It has plateaued"* needed 3500 more
  iterations of evidence than the samples that suggested it.
- **`gci.py` refuses rather than reports** where Celik *et al.* exclude the case: it prints
  `NOT CONVERGENT: non-monotone` instead of a GCI, and flags `p` outside [0.5, 2.5] instead
  of quoting it silently.
- **Cell counts are read from `blockMesh`'s own log**, not recomputed from `gen_case.py`,
  and the code asserts they scale as `refine²` before applying the procedure — so a grading
  change could not silently break the family.
