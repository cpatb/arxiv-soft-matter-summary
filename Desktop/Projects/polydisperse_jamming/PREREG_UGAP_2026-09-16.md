# PREREG — the U-GAP block: measuring a(u), b(u) where the K-ladder reads them

**Written 2026-09-16, committed BEFORE any `phase12_ug*` file exists.** Blind predictions
are printed by `predict_ugap.py` → `out/predict_ugap.txt`, which **refuses to run if any
U-GAP data file is present**. Grade with `grade_ugap.py` → `out/grade_ugap.txt`.

Tag `ug` (N = 1500) and `ug3k` (N = 3000). These match **none** of `mixr*`, `extr*`,
`tern*`, `kl*`, so M0–M6, F1–F3, X0–X5, T0–T6, NA1–NA4 and G0–G8 cannot be re-decided by
anything here; `grade_ugap.py` asserts that structurally.

---

## 1. What this block is for, and what already decided its design

`CLAIMS_KLADDER_2026-09-15.md` C1a(ii): the tabulated closure has **no measured knot above
u = 0.667**, and **4 of K = 5's 10 pairs — its most numerous pair type — sit at
u = 4^(−1/4) = 0.70711**, inside that gap. The rung that decides "and then SATURATES" is
evaluated where the closure has no data.

Two diagnostics, both run and committed before this file was written, turned that caveat
into a design.

**`out/diag_ugap.txt` — the region MATTERS; it is not a degenerate parameter.** Sweeping
a(u\*) and b(u\*) over the monotone-admissible range `[a(0.667), 1] × [b(0.667), 1]`:

| quantity | measured |
|---|---|
| worst K = 5 joint φ excursion | **9.59 scale** (4.7× TAB's own in-sample rms 2.02) |
| worst K = 2 / K = 3 excursion (control) | **0.00 scale**, exactly |
| worst K = 4 excursion | 3.49 scale — through the *interpolant's shape* at u = 0.630, not through a pair at u\* |
| rows where the sweep carries R(5) − R(4) across zero | **12 of 12** |
| candidate separation (excursion / median seR 0.48) | **20.09×**, against the standing 1.6× bar |

This is the mirror of the ternary block's withdrawn `a(1/4)` inversion, and it comes out
the **opposite** way: there φ was bit-identical over the whole sweep, so the parameter
could not testify; here it moves 4.7× the in-sample rms, so it decides.

**`out/diag_uident.txt` — why r = 1.25 came back UNDETERMINED, and what does *not* fix it.**

* The profile bar is Δχ² = 1 **per point**, so the interval's width is set by the
  **curvature**, not the sample size. Demonstrated rather than argued: the same r = 1.5
  point set at n = 8, 16, 32 returns w(a) = **0.0400** three times, identically.
  **More seeds and more compositions at the same sensitivity buy nothing under this rule**,
  and seeds buy less than that — the bar is in units of `PHI_SCALE`, not of the seed error.
* The closure is a **MIN over branches** with an escape at a = 1 or b = 1 that collapses it
  to the **constant φ = β**. At r = 1.25 that escape scores χ² = **7.09 against a bar of 8**,
  so every (a, b) that kills the binding branch fits. The joint 95 % region there covers
  a ∈ [0.0001, 0.9999] **and** b ∈ [0.0001, 0.9999] and **contains the constant model**, so
  **u = 0.8 constrains neither parameter nor any combination of them** — the fallback of
  "constrain the combination the K-ladder uses" is not available at that ratio.
* `a` enters only the coarse branch and `b` only the fine one, and they swap at the
  crossover x\* ≈ 0.356. The block that ran sampled x ∈ [0.05, 0.60] — **6 points below,
  2 above**. A grid that straddles x\* is what the identifiability ceiling prefers; a grid
  *above* it (x ∈ [0.30, 0.75]) is far worse and fails `a` at **every** ratio tested.

**Therefore the ratio is 1.41421, not 1.25.** K = 5 reads the closure at u = 0.70711. The
noise-free ceiling is already identifiable there on the *old* grid and fails only at
r = 1.25, where the filling signal max(φ − β) is **1.4 scale** against 7.1 at r = 1.5. The
block that would rescue u = 0.8 is not the block the K-ladder needs.

---

## 2. The design — which axis each cell moves, and what is held fixed

**Moved:** size ratio `r` (three levels) and fines volume fraction `x` (nine levels).
**Held fixed and identical to every binary cell on disk:** `p_hat` = 1e-3, `mu_p` = 0.3,
seeds (1, 2, 3), `MIN_PHASE_COUNT` = 50, the protocol in `packing.build_and_compress`, and
the runner — `phase12_mixture.py` **unchanged**, driven by its existing
`--ratio` / `--xgrid` / `--tag` options.

**N = 1500 is not a choice.** `mixture_closure.binary_points` drops any record whose
`n_requested` is not 1500, so a knot measured at 3000 could never enter the closure. The
N = 3000 arm is a separate, labelled control (§5, U5) and is **never pooled**.

| arm | r | u | why |
|---|---|---|---|
| **A** | 4^(1/4) = 1.41421 | **0.70711** | the knot the K-ladder actually reads |
| **B** | 1.25 | 0.80000 | the far end of the gap, and the direct retest of the ratio that failed |
| **C** | 1.5 | 0.66667 | **the control** — an already-identifiable knot, re-measured on the *new* grid |

Arm C is what makes arms A and B attributable. Without it, a disagreement between the new
a(0.70711) and the shipped PCHIP could be the ratio **or** the composition grid, and the
block would not be able to say which.

**Composition grid** (`ugap.X_GRID`), 9 points, four below x\* and five above:

```
x = 0.05, 0.12, 0.20, 0.28, 0.36, 0.45, 0.55, 0.70, 0.88
```

`x` = 0.05 and 0.20 are kept deliberately: r = 1.25 and r = 1.5 already have **those exact
cells** on disk under tag `mixr`, at the same seeds, N, `p_hat` and `mu_p`, and
`phase12_mixture.cell_key` carries no tag — so those four cells must come back
**bit-identically** (gate U0). `x` = 0.88 is the largest composition affordable at
`MIN_PHASE_COUNT` on **all three** ratios at N = 1500 (the binding one is r = 1.5, limit
0.896).

**27 cells, 27 admitted, 0 refused, 81 runs at N = 1500**; arm A repeated at N = 3000 on
the **identical cell list** (9 cells, 27 runs).

### Why this block does NOT buy more seeds

Because the registered rule cannot see them. §1 measured that. Spending the budget on the
composition axis instead is the whole point of running `diag_uident.py` first.

---

## 3. What "identifiable" means — quoted verbatim, not restated

From `grade_mixture_ext.tabulated_closure`, unchanged, and reused in `ugap.identifiable`:

> Identifiable = both profile intervals narrower than 0.25 AND containing their own best
> fit.

`ugap.IDENT_WIDTH = 0.25`. **Not movable.** A block designed to meet a bar must not be
graded against a bar it moved. `grade_ugap.py` asserts the constant equals the one
`grade_mixture_ext` uses.

---

## 4. The identifiability ceiling of THIS grid, computed before the run

Noise-free φ generated from the shipped closure at each arm's own u, passed through the
rule above (`out/predict_ugap.txt`). This is an **upper bound on the block, not a promise
about it**.

| arm | r | u | x\* | pts below / above | w(a) | w(b) | ident? |
|---|---|---|---|---|---|---|---|
| A | 1.41421 | 0.70711 | 0.3565 | 4 / 5 | 0.055 | 0.055 | YES |
| B | 1.25000 | 0.80000 | 0.3617 | 5 / 4 | 0.115 | 0.085 | YES |
| C | 1.50000 | 0.66667 | 0.3539 | 4 / 5 | 0.050 | 0.050 | YES |

**Robustness — how wrong the shipped a(u\*) can be before arm A fails:** arm A stays
identifiable up to a true **a(u\*) = 0.9800** (against the shipped 0.9398). It fails at
0.99. The old 8-point grid's limit is 0.960, which is why the grid was widened.

---

## 5. The registered gates

Bars are **TAB's own in-sample binary residual statistics** (rms **2.02**, max **5.74**)
and the verbatim identifiability rule. **No new bar is invented anywhere in this file.**

**U0 — wiring, RUN not asserted.** The four (r, x) cells the `mixr` block already holds —
(1.25, 0.05), (1.25, 0.20), (1.5, 0.05), (1.5, 0.20) — must reproduce **bit-identically**,
`diff == 0.0` exactly, seed by seed. *Predicted 12 of 12 runs.* A failure means the new
block is not on the old protocol and **nothing else in the block is graded**.
(The last prereg asserted a reduction "by inspection" and it raised `IndexError`; this one
runs.)

**U1 — is a(u) measurable at u = 0.70711?** Apply the verbatim rule to arm A.
*Predicted **IDENTIFIABLE**.* This is the question the block was commissioned to answer and
it can fail outright.

**U1b — the attribution control.** Arm C must be identifiable **and** its (a, b) must agree
with the existing `mixr` values a = 0.9249, b = 0.8985 to within the **union of the two
profile intervals**. *Predicted PASS.* If U1b fails, U1 and U1c are reported as
**NOT INFORMATIVE about the ratio**, because the composition grid moved the answer on a
ratio where nothing else changed.

**U1c — the retest of the failed ratio, and a registered reading of its null.** Apply the
rule to arm B. *Predicted **IDENTIFIABLE**, with lower confidence than U1.* Registered in
advance: **if arm B comes back UNDETERMINED on a grid whose ceiling passes at a true
a = 0.98, that is evidence a(0.8) > 0.98**, i.e. the filling at u = 0.8 is weaker than any
value the closure currently carries. A null here is a bound, not a shrug.

**U2 — the DIRECTION, and it is the prediction most likely to fail.** r = 1.25 is held out
of the closure (it failed identifiability, so `tabulated_closure` never put a knot there),
so its φ data is a legal test of the extrapolation — and the shipped closure
**over-predicts filling there on 8 of 8 compositions**, mean **−1.03 scale**, rms 1.19,
against a fitted-knot control at r = 1.5 of mean +0.11, 3/8 negative.
*Predicted: the shipped closure's residual at arm A is **negative on ≥ 7 of 9**
compositions* (sign-test p ≤ 0.09), *and the fitted a(0.70711) comes out **above** 0.9398.*
Both halves are reported; the sign count is the gate.

**U3 — the K-ladder consequence, and this prereg predicts its own claim is at risk.**
Reinsert the measured knot and recompute R(5) − R(4) on the existing K-ladder cells.
*Predicted: the saturation verdict SURVIVES — R(5) − R(4) ≤ 0 on ≥ 2 of 3 `y_S` levels at
both N*, which is the pattern C1 measured. **The registered blind table says this fails if
a(u\*) ≳ 0.97**:

| N | y_S | R5−R4 now | at a = 0.9398 | at a = 0.97 | at a = 0.99 |
|---|---|---|---|---|---|
| 1500 | 0.03 | +0.20 | +0.45 | +1.94 | +2.84 |
| 1500 | 0.05 | −1.58 | −1.33 | +0.55 | +1.66 |
| 1500 | 0.08 | −1.77 | −1.52 | +0.19 | +1.16 |
| 3000 | 0.03 | +0.63 | +0.88 | +2.37 | +3.27 |
| 3000 | 0.05 | −0.80 | −0.55 | +1.33 | +2.45 |
| 3000 | 0.08 | −1.21 | −0.96 | +0.76 | +1.73 |

U2 predicts a > 0.9398 and U3 predicts saturation survives; **they are in tension by
construction and the block decides between them.** The a = 0.9398 column re-inserts the
closure's own interpolated value as an explicit knot, so the difference from "now"
(+0.25 on every row) is PCHIP reshaping and not physics; it is printed so the two cannot
be confused. **The K3→K4 step and the G7 sign result are untouched by this block** — K = 2,
3, 4 pairs all sit at u ≤ 0.630 — except through the interpolant leakage of caveat C-iii.

**U4 — monotonicity.** With the new knots inserted, a(u) and b(u) must remain monotone
non-decreasing in u. *Predicted PASS.* A failure is a statement about the two-parameter
tabulated form, not about the measurement.

**U5 — the N = 3000 arm, and what it is NOT.** Declared in advance: the K-ladder measured a
**block-wide +1.9 scale offset** when N doubles (G5), so a knot fitted at N = 3000 is **not
comparable** to a closure calibrated at N = 1500, and a difference in a between the two N
is **expected** and is **not** evidence about the closure. What U5 grades is only that arm
A's identifiability is not a finite-size accident: *predicted identifiable at N = 3000 too*.

**U6 — POST-HOC, declared here so it cannot be presented as registered.** Whether the new
knot collapses G8's 11.9× nested-rule convention spread (J 34.9, ADDtab 108.1, SEQd21 42.8,
SEQd32 48.0, SEQd10 186.3, SEQd43 509.3). Reported as a search, not a test.

---

## 6. What this block cannot decide — carried from `out/predict_ugap.txt`

These four caveats are printed by the blind prediction file and **must appear in the claims
file**; `grade_ugap.py` asserts they do. (On 2026-09-15 `predict_kladder.txt` flagged this
very gap before any cell ran and the claims shipped without it until a referee asked.)

* **C-i** One friction (`mu_p` = 0.3), one pressure rung (`p_hat` = 1e-3), one N (1500).
  Nothing here says a(u) is pressure- or friction-independent.
* **C-ii** Arm A's knot is at u = 0.70711 **exactly**. It does **not** fill the gap: the
  interval (0.70711, 1) still has no measured point, and PCHIP's shape there remains an
  extrapolation between a measured knot and a proved limit.
* **C-iii** Inserting a knot moves a(u), b(u) at **neighbouring** u through the
  interpolant's derivatives — measured at up to **3.49 scale** on the K = 4 cells, whose
  pairs all sit at u ≤ 0.630. Any K = 4 number recomputed after this block is **not**
  independent of it.
* **C-iv** The identifiability bar is in units of `PHI_SCALE`, not of the seed error. It
  does not tighten with seeds or compositions, so **a PASS is not a claim about precision.**

## 7. Rules this block is run under

* A failure is the model's range, **not** a refit. No bar moves.
* Read the **exit code**, never the artifact, when the artifact is what the run rewrites.
  Never run two copies of `grade_ugap.py` concurrently over one output file.
* Print N before and after every filter.
* A step in a residual is (data − model) and **both move**: any R quoted across the K axis
  comes with its decomposition.
* Before quoting a null, check it is not a restatement of an axis that was varied.
* Re-read each registered **sentence** beside the line of code that grades it.
