# CLAIMS — `mustarcap`: the mixed-cap repair. **ALL GATES PASS.**

Graded 2026-09-22 against `PREREG_MUSTARCAP_2026-09-22.md`, written before launch.
Artifact: `out/grade_mustarcap.txt`. **108 of 108 records delivered.**

---

## The one-line result

> **The `μ` ladder's monotonicity violation was a cap artefact, and the repair proves it
> three independent ways — including by turning two of three MISSED blind predictions
> into HITS without touching a tolerance.**

---

## 1. The gates

| gate | registered | delivered | verdict |
|---|---|---|---|
| **M1a IDENTITY** (primary) | previously-converged cells bit-identical at `1e6` | **55 compared, 0 mismatches** | ✅ PASS |
| **M2 RESCUE** | ≥ 15 of 17 previously-capped cells converge | **17 of 17** | ✅ PASS |
| **M3 MONOTONICITY** (the point) | no rise in `φ(μ)` above `PHI_SCALE` at any CV | **0 rises** | ✅ PASS |
| **M4 DIRECTION** | previously-capped cells shift **UP** | **+0.01319 mean = +6.28 × `PHI_SCALE`, 17 of 17 positive** | ✅ PASS |
| M5 `μ*`, `A` | reported, not gated | `A(0.5)/A(0)` = 0.886 vs [0.90, 1.20] | **UNDETERMINED**, as registered |

## 2. The ladder, before and after

Matched rung `p̂ = 3e-5`, the fit-free channel:

| CV | μ=0.0 | μ=0.3 | μ=0.6 | μ=1.0 | μ=3.0 |
|---|---|---|---|---|---|
| **before** (mixed cap) | 0.69084 | 0.63667 | **0.59853** (0/3) | 0.61486 | **0.57693** (0/3) |
| **after** (all 1e6) | 0.69084 | 0.63667 | **0.62181** (3/3) | 0.61486 | **0.61236** (3/3) |

**The `μ = 0.6` cell moved +0.0233 and the `μ = 3.0` cell +0.0354** — both upward, both
from 0/3 converged to 3/3. `μ = 0.3` and `μ = 1.0` were not re-run and did not move.

> 🔑 **The +0.01633 rise at `μ` 0.6 → 1.0 — the ladder's single monotonicity violation, at
> 7.8 × `PHI_SCALE` — is GONE.** It was entirely the under-compression of a capped cell.

## 3. Why this is a demonstration and not just a repair

Three things had to be true together, and each could have failed:

1. **M1a: the cells that had converged did not move at all** (55/55, exactly zero). A cap
   cannot affect a run that never reached it — so whatever moved, the cap is why.
2. **M4: every rescued cell moved UP**, 17 of 17, mean +6.28 × `PHI_SCALE`, range
   +0.00023 to +0.03946. This is the *directional* prediction: a capped cell stopped
   before reaching its pressure target, so it is **under-compressed and reads low**.
   Monotonicity alone (M3) is necessary but not sufficient; M4 is what distinguishes
   "the cap explained it" from "it happens to be monotone now".
3. 🔑 **F4: all three of the ORIGINAL blind point predictions now HIT.**

| blind prediction (written 2026-09-19) | mixed-cap | **repaired** |
|---|---|---|
| `φ(μ=0.6, CV=0.0)` = 0.56394 | +0.00019 HIT | +0.00019 **HIT** |
| `φ(μ=0.6, CV=0.5)` = 0.62042 | −0.02189 🔴 MISS | **+0.00139 HIT** |
| `φ(μ=3.0, CV=0.5)` = 0.61346 | −0.03653 🔴 MISS | **−0.00110 HIT** |

**The ±0.0030 tolerance was never widened.** Two predictions that missed by 7× and 12×
the tolerance now sit inside it. The predictions were right; the protocol was wrong.

## 4. What did NOT change, and must not be claimed

- 🔴 **`μ*` is still UNDETERMINED.** `A(0.5)/A(0)` = 0.886 remains outside [0.90, 1.20],
  so **`μ*` and `A` are not separated by this design** whatever `μ*` does — exactly the
  outcome `PREREG_MUSTAR` registered as *likely*. **Repairing a cap does not repair an
  identifiability problem**, and this campaign was never registered as able to.
  (`μ*(0.5)/μ*(0)` = 1.687 against [1.05, 1.45] is reported and is moot.)
- ⚠ **The algebraic form's CV = 0.5 fit is DEGENERATE** — `A = 12326`, `μ* = 0.0000` on
  4 points and 3 parameters. Its ratio would have printed as `1.6 × 10⁵`. The grader now
  checks each fitted parameter against its own admissible range and **refuses the ratio**
  rather than quoting a runaway fit.
- ⚠ **2 of 108 cells still cap at `1e6`**: `μ_p = 0.0, CV = 0.5`, seeds 1 and 3, bottom
  rung. **They are reported, not rescued** — the prereg says a third cap raise is not the
  response. They are frictionless *and* widest-CV, the stiffest corner, and they were
  already capped at `1e6` before this campaign. **So the CV = 0.5 ladder's top rung is
  still censored (1/3 converged)** and M5's CV = 0.5 fit correctly uses 4 points, not 5.
- ⚠ **`1e6` was barely enough.** The rescued cells used 209 501 – **981 426** iterations
  and **10 cells needed more than 500 000**. A campaign at this corner should not assume
  headroom.

## 5. The protocol conclusions this pair of campaigns establishes

- **`MAX_ITER = 200000` is too low wherever `CV > 0` and `μ_p ≥ 0.6`** — 0 of 3 at three
  separate cells, now all rescued.
- **Capped cells bias `φ` DOWNWARD**, measured twice independently: +0.0049 on
  `wcbidicap`'s headline difference, and **+0.01319 mean over 17 cells** here.
- **The identity argument for raising a cap is now tested three times and passed three
  times**: 58/58 (`wcbidicap` C1), 55/55 (M1a), and 52/52 on partial data.
- ⚠ **The campaign that failed was the one deliberately run at the default cap** — the
  comment above `mustar`'s block in `scripts/launch_capmustar.sh` reads *"Back to the
  DEFAULT cap: this campaign is not about caps."* It ran in the same launch as the
  campaign that had just proved the default inadequate.

## 6. Free by-product: a reproducibility measurement on a published number

`mustarcap` includes `μ_p = 0, CV = 0, lognormal` — the twelve cells
`\phiJFrictionless` is built from. **All twelve reproduce to 1.2 × 10⁻¹⁴**, with `iters`
and `n_contacts` exactly equal. `test_anchor_reproduction.py` pins that, and also pins
**non-takeover**: `load_fits` dedups with `idxmax`, which returns the first row at a tie,
so which campaign feeds the published macro is decided by filename sort order. It has not
moved and now cannot move silently.
