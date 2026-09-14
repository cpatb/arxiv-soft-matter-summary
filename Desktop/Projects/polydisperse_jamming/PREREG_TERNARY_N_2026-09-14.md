# PREREG_TERNARY_N_2026-09-14 — is the ternary under-prediction a FINITE-PHASE-SIZE artifact?

**Addendum to `PREREG_TERNARY_2026-09-14.md`, written after T1 failed and before any
N = 3000 cell existed.** It exists because grading T1 surfaced a confound that the
original design could not break, and the honest response to that is a run, not a
paragraph.

---

## 1. The problem

T1 failed: the additive closure under-predicts all six ternaries, mean **+7.53** scale,
rms **7.83** against a bar of 2.8. The registered reading is non-additivity.

But the six residuals are **perfectly rank-correlated with the number of COARSE grains**:

| (`y_L`,`y_M`,`y_S`) | `n_L` at N = 1500 | TAB residual (scale) |
|---|---|---|
| (0.55, 0.30, 0.15) | 66 | +10.1 |
| (0.60, 0.25, 0.15) | 74 | +9.9 |
| (0.65, 0.20, 0.15) | 82 | +8.5 |
| (0.55, 0.35, 0.10) | 84 | +6.7 |
| (0.60, 0.30, 0.10) | 96 | +5.7 |
| (0.70, 0.20, 0.10) | 121 | +4.3 |

**Spearman(residual, `n_L`) = −1.0000 (p = 0.0000), on 6 of 6.** Also
Spearman(residual, `y_S`) = +0.878.

All six clear `MIN_PHASE_COUNT` = 50, and `d_max/L` = 0.187–0.212 is below the 0.25
flag — so nothing was violated. But `n_L` and `y_S` are **confounded by construction**:
number fraction goes as `y_i/d_i³`, so raising the fines volume at fixed N necessarily
strips grains from the coarse class. Within this design the two cannot be separated, and
a perfect rank correlation with a sample-size variable is exactly what a finite-size
artifact looks like.

**So T1's verdict is currently ambiguous between:**

* **(1) NON-ADDITIVITY** — the residual is set by composition (`y_S`), and `n_L` merely
  tracks it.
* **(2) FINITE PHASE SIZE** — the residual is set by how few coarse grains there are, and
  composition merely tracks it. A coarse skeleton of 66 grains in a periodic cell may
  simply pack differently from the thermodynamic limit the closure describes.

## 2. The discriminator

**Hold the composition FIXED and double N.** That moves `n_L` by exactly 2× while leaving
`y_L`, `y_M`, `y_S`, CV, `d32`, every size ratio and every pair weight **identical**. It
is the one axis the original design holds fixed and the confound needs.

* N = **3000**, all six registered compositions plus the monodisperse control, seeds
  1/2/3 = **21 cells**, tag **`tern3k`** (matches neither `phase12_tern_shard*` nor
  `mixr*`/`extr*`, so T0–T6, X0–X5 and M0–M6 are all closed against it).
* `n_L` goes 66→132, 74→148, 82→164, 84→168, 96→192, 121→242.
* `d_max/L` **falls** as N^(−1/3), 0.187–0.212 → 0.148–0.168, further below the flag.
* `beta` is **re-measured at N = 3000** from that block's own monodisperse control and
  the N = 3000 residuals are computed against it. Using the N = 1500 `beta` would fold a
  finite-size shift in the reference into the very quantity under test.

## 3. Registered predictions and bars

Residual sign convention unchanged: **MEASURED − PREDICTED**, positive = under-predicted.

| Gate | Statement | Bar |
|---|---|---|
| **NA0** | The N = 3000 monodisperse control is a valid `beta` for its own block | reported, not barred; the N-shift in `beta` is printed |
| **NA1** | **THE DISCRIMINATOR.** Mean TAB residual over the six compositions at N = 3000 | **NON-ADDITIVE ⇔ ≥ +4.7** (within 2/3 of the N = 1500 value of +7.53). **FINITE-SIZE ⇔ ≤ +2.8** (the T1 bar: at that point the additive closure passes and T1's verdict is withdrawn). Between 2.8 and 4.7 ⇒ **UNDECIDED**, reported as such |
| **NA2** | Per-composition agreement between N = 1500 and N = 3000 | \|Δresidual\| ≤ 2× the standard error of the difference of the two three-seed means, reported per composition with the count that pass |
| **NA3** | **The correlation itself.** Spearman(residual, `n_L`) recomputed over the **12** points (6 compositions × 2 N values) | if the residual is set by `n_L`, the pooled correlation stays strongly negative; if it is set by composition, the N = 3000 points sit at the SAME residuals with DOUBLED `n_L`, which **destroys** the pooled correlation. Registered: NON-ADDITIVE ⇔ \|rho_pooled\| ≤ 0.6 |
| **NA4** | Spearman(residual, `y_S`) within the N = 3000 block alone | reported; `y_S` is the surviving candidate if `n_L` is excluded |

**NA1 is written so it can overturn my own headline.** If the mean residual falls to the
T1 bar, the ternary block does not show non-additivity, `BENCH_HYPOTHESES.md` §1 item D
and §2 are withdrawn, and the finding becomes "the additive closure is sound and the
earlier failure was a coarse-phase sampling effect" — which is a publishable negative and
a different next experiment.

**No bar in `PREREG_TERNARY_2026-09-14.md` is touched.** T1 keeps its verdict at
N = 1500 whatever this block says; NA1 decides what that verdict MEANS.

## 4. What this still cannot settle

Doubling N doubles every class's count, so it cannot distinguish "too few COARSE grains"
from "too few grains overall". It is not meant to: the coarse class is the one with 66,
the fine class has 1147, and no finite-size story runs through a class with a thousand
members. If NA1 comes out non-additive, a residual finite-size dependence at the 1–2
scale level is not excluded, and NA2 measures it.
