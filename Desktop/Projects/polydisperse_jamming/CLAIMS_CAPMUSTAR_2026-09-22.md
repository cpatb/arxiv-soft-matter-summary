# CLAIMS — `wcbidicap` and `mustar`, the two campaigns launched 2026-09-19

Graded 2026-09-22. Artifacts: `out/grade_wcbidicap.txt`, `out/grade_mustar.txt`.
Preregs: `PREREG_WCBIDI_2026-09-19.md` Amendment 3, `PREREG_MUSTAR_2026-09-19.md`.

Both campaigns completed: **72 of 72 records each.**

---

## 1. `wcbidicap` — the cap re-run. **ALL REGISTERED GATES PASS.**

| gate | registered | delivered | verdict |
|---|---|---|---|
| C1 IDENTITY (primary) | previously-converged cells bit-identical at `1e6` | **58 compared, 0 mismatches** | ✅ PASS |
| C2 RESCUE | the 14 previously-capped cells converge | **14 of 14** | ✅ PASS |
| C3 ARGMAX (the point) | argmax of `φ_0^bidi(CV)` in {0.2, 0.3, 0.4} on converged-only data | **CV = 0.4 on both channels**, gaps +0.00357 (1 s.e. 0.00158) and +0.00439 (1 s.e. 0.00091) | ✅ PASS |

**C1 is the gate that licenses the rest**, and it is the second time this argument has been
tested in this repository and passed: a cell that stopped because it *converged* cannot be
affected by raising a limit it never reached. The 14 rescued cells used **205 276 – 475 581**
iterations — every one under half the new cap, and every one above the old one.

### 🔑 CLAIM W1 — the provisional (T) verdict is now a RESULT

`PREREG_WCBIDI` §5 declared the first campaign's argmax **worthless** if more than 2 of 72
cells capped; 14 did, and the argmax rung had **0 of 3** fully-converged ladders. After the
re-run, **all 6 CV levels have 3 of 3 fully-converged ladders** and the argmax is unmoved at
CV = 0.4.

> **The frictionless branch inherits the bidisperse family's non-monotonicity**, as
> translation requires. The registered alternative — a monotone frictionless branch, or an
> argmax at CV = 0.5 that would have meant the turnover was a censoring artefact — is not
> what the data show.

### 🔴 CLAIM W2 — G4's blind band was passing on CENSORED data, and the miss is now visible

`Δφ_0` (the CV 0 → 0.5 shift, fit-free rung channel, the channel Amendment 2 registered as
the one to quote):

| data | `Δφ_0` | registered band [+0.0125, +0.0146] |
|---|---|---|
| first campaign, mixed convergence | +0.01406 | inside |
| **raised cap, all converged** | **+0.01895 ± 0.00100** | 🔴 **OUTSIDE** |

**The band is not moved.** C1 proves the 58 converged cells did not move at all, so the
entire **+35 %** shift comes from the 14 cells that had capped — 13 of them at the bottom
rung. **A capped cell stopped before reaching its pressure target, so it is
under-compressed and its `φ` reads LOW.** This is the direct, quantitative measurement of
that bias, and it is the evidence the `mustar` diagnosis below rests on.

⚠ C4 is a **diagnostic, not a gate**, in both the prereg and the grader — it is printed and
never appended to the failure list. The grader prints `ALL GATES PASSED` while C4 visibly
misses, which is correct behaviour but is worth stating plainly rather than leaving to be
noticed.

---

## 2. `mustar` — **NOT READABLE.** The verdict its own gates give is right but shallow.

| gate | registered | delivered | verdict |
|---|---|---|---|
| F1 WIRING | the re-run `μ = 3.0, CV = 0` cells reproduce an existing record | 3/3 | ✅ PASS (but see §4) |
| F2 HEALTH | 72 of 72 converge | **55 of 72** | 🔴 FAIL |
| F3a | `A(0.5)/A(0)` in [0.90, 1.20] | **1.510** | 🔴 FAIL |
| F3b | `μ*(0.5)/μ*(0)` in [1.05, 1.45] | **2.219** | 🔴 FAIL (moot: F3a failed) |
| F4 | three blind points within ±0.0030 | **1 hit, 2 miss** (−0.02189, −0.03653) | 🔴 FAIL |
| F6 | `μ*(CV)` monotone | 0.188, 0.216, 0.416 | ✅ PASS |

`F3a` returns **UNDETERMINED as the prereg registered as likely**: if `A` moves as much as
`μ*`, the two mechanisms are not separated by this design whatever `μ*` does.

### 🔑 CLAIM M1 — the binding defect is that the μ ladder is a **MIXED-CAP LADDER**

`grade_mustar` fits `φ(μ) = φ_∞ + A exp(−μ/μ*)` along five frictions. Rebuilding that ladder
through the grader's own loader, with the cap and convergence printed beside every point:

| CV | μ=0.0 | μ=0.3 | μ=0.6 | μ=1.0 | μ=3.0 |
|---|---|---|---|---|---|
| 0.0 | 0.63182 (3/3, 1e6) | 0.57535 (3/3, 2e5) | 0.56413 (3/3, **2e5**) | 0.56148 (3/3, 1e6) | 0.56054 (3/3, **2e5**) |
| 0.2 | 0.64202 (**0/3**, 2e5) | 0.58755 (3/3, 2e5) | 0.57824 (**1/3**, 2e5) | 0.57316 (3/3, 1e6) | 0.56795 (**0/3**, 2e5) |
| 0.5 | 0.69084 (**1/3**, 1e6) | 0.63667 (3/3, 1e6) | 0.59853 (**0/3**, 2e5) | 0.61486 (3/3, 1e6) | 0.57693 (**0/3**, 2e5) |

Two statements follow, and neither needs a threshold beyond the declared `PHI_SCALE`:

1. **The protocol changes along the very axis the fit measures a scale in.**
   `μ ∈ {0.6, 3.0}` came from the `mustar` tag at the default `MAX_ITER = 200000`;
   `μ ∈ {0, 1.0}` came from `capfix_lg*` at `1000000`.
2. **The ladder's ONE monotonicity violation sits exactly at the cap boundary.**
   `φ_J` cannot rise with friction. Across all fifteen points there is exactly one rise:
   **CV = 0.5, μ 0.6 → 1.0, +0.01633 = 7.8 × `PHI_SCALE`** — the step from a
   **0/3-converged cell at 2e5** to a **3/3-converged cell at 1e6**, in the direction
   under-compression predicts and W2 measured.

> **So no fitted `μ*` from this ladder is a measurement of anything.** At the matched rung,
> **0 of 3 cells converged** at (0.6, 0.5), (3.0, 0.2) and (3.0, 0.5) — three of the six
> `(μ, CV)` cells with `CV > 0` that the whole `μ*(CV)` question is about. Both of F4's
> misses are at CV = 0.5.

### CLAIM M2 — the non-convergence tracks BOTH axes under test

Capped cells per `(μ, CV)`, of 12 each: (0.6, 0.2) 2, (0.6, 0.5) 6, (3.0, 0.2) 3,
(3.0, 0.5) 6. **Zero at CV = 0 for either friction.** This is the mechanism already on the
card — friction and polydispersity both raise `φ_J`, so at a fixed pressure ladder the
packing is stiffer and the solver cannot reach the bottom rung inside the cap — with the
extra fact that here it does not merely narrow the ladder, it **fails to converge at all.**

### ⏳ REPAIR REGISTERED AND RUNNING

`PREREG_MUSTARCAP_2026-09-22.md`: the same grid at `PJ_MAX_ITER = 1e6`, 108 records,
graded by `grade_mustarcap.py` (written before the cells landed). **M1a is the primary
identity gate**; **M3, the monotonicity test, is not passable by construction**, because
only the `μ = 0.6` end of the violating step is being re-run.

---

## 3. What these two campaigns together establish about the PROTOCOL

- **`MAX_ITER = 200000` is too low wherever `CV > 0` and `μ_p ≥ 0.6`.** Not marginally:
  0 of 3 at three separate cells.
- **Capped cells bias `φ` DOWNWARD**, measured at +0.0049 on a headline difference (+35 %).
- **The identity argument for raising a cap is now tested twice and passed twice**
  (58/58 here; 52/52 so far on the partial `mustarcap` data).
- ⚠ **A campaign deliberately run "at the default cap because this campaign is not about
  caps"** (`scripts/launch_capmustar.sh`, the comment above the `mustar` block) was the one
  that failed. The two campaigns ran in the same launch, and the first one proved the
  default inadequate hours before the second one used it.

## 4. ⚠ A DEFECT IN `grade_mustar.py`'s F1, RECORDED AND NOT SILENTLY FIXED

F1 is named *"the re-run `μ_p = 3.0, CV = 0` cells reproduce an existing record **exactly**"*
and the prereg says **BIT-IDENTICAL**. It is implemented as `abs(v - x) < 1e-12`.

Checking the stated criterion instead: the re-run records and the pre-existing ones differ
in **`phi` only**, by **1.1e-16 to 1.2e-14** (≤ ~99 ULP), while `iters`, `n_contacts`,
`z_all`, `z_nr`, `rattler_frac` and `p_hat` are **exactly identical**. The trajectory is the
same; only the derived scalar moves.

> **F1's VERDICT stands — the wiring is right — and its WORDING does not.** The fact it was
> concealing is a real one and is now on record: **this repository's older records are not
> bit-reproducible by the current code at the ~1e-14 level**, while same-code-era runs are
> bit-identical. `grade_mustarcap.py` reports this as a **diagnostic, not a gate**, because
> it is not the registered criterion and failing a campaign on it would fail it for the
> wrong reason.
