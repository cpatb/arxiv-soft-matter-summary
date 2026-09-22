# PREREG — `mustarcap`: the μ ladder is a MIXED-CAP ladder, and that is why it bends

Written 2026-09-22, **before launch**, after grading `wcbidicap` and `mustar`
(`out/grade_wcbidicap.txt` route, `/tmp/gmu.log` route; both committed as
`CLAIMS_MUSTARCAP_2026-09-22.md` inputs).

## 1. The defect this campaign exists to repair

`grade_mustar.py` fits `φ(μ) = φ_∞ + A exp(−μ/μ*)` along a five-point friction ladder at the
matched rung `p̂ = 3e-5`. Rebuilding that ladder **through the grader's own loader** and
printing the iteration cap and convergence status beside every point gives:

| CV | μ=0.0 | μ=0.3 | μ=0.6 | μ=1.0 | μ=3.0 |
|---|---|---|---|---|---|
| 0.0 | 0.63182 (3/3, 1e6) | 0.57535 (3/3, 2e5) | 0.56413 (3/3, **2e5**) | 0.56148 (3/3, 1e6) | 0.56054 (3/3, **2e5**) |
| 0.2 | 0.64202 (**0/3**, 2e5) | 0.58755 (3/3, 2e5) | 0.57824 (**1/3**, 2e5) | 0.57316 (3/3, 1e6) | 0.56795 (**0/3**, 2e5) |
| 0.5 | 0.69084 (**1/3**, 1e6) | 0.63667 (3/3, 1e6) | 0.59853 (**0/3**, 2e5) | 0.61486 (3/3, 1e6) | 0.57693 (**0/3**, 2e5) |

Two facts follow, and neither needs a threshold:

1. **The ladder mixes two caps, and the cap correlates with the axis under test.**
   `μ ∈ {0.6, 3.0}` came from the `mustar` tag at the default `MAX_ITER = 200000`;
   `μ ∈ {0.0, 1.0}` came from `capfix_lg*` at `1000000`. The fitted quantity is `μ*`, a
   scale *along μ*, and the protocol changes along μ.
2. **The ladder's one monotonicity violation sits exactly at the cap boundary.**
   `φ_J` cannot rise with friction. Across all fifteen points there is exactly one rise:
   CV = 0.5, μ 0.6 → 1.0, **+0.01633 = 7.8 × `PHI_SCALE`**, i.e. resolved well above the
   protocol's own run-to-run scatter — and it is precisely the step from a **0/3-converged
   cell at 2e5** to a **3/3-converged cell at 1e6**.

**The sign is the one under-compression predicts.** A capped cell stopped *before* reaching
its pressure target, so it is under-compressed and its `φ` reads **low**. `wcbidicap`
measured that bias directly this same session: rescuing its 14 capped cells moved a headline
difference from **+0.01406 to +0.01895** (+35 %), capped cells reading low.

> Therefore the `mustar` verdict is not "UNDETERMINED because `A` moved". `F3a` is correct
> but shallow. The binding statement is that **the object being fitted is not a
> single-protocol ladder**, so no fitted `μ*` from it is a measurement of anything.

**Registered grid.** `lognormal`, `μ_p ∈ {0.0, 0.6, 3.0}` × `CV ∈ {0.0, 0.2, 0.5}` ×
seeds {1,2,3} × 4 rungs = **108 records**, tag **`mustarcap`**, `PJ_MAX_ITER = 1000000`.
The whole grid rather than only the 17 capped cells, for the reason Amendment 3 of
`PREREG_WCBIDI` gives: mixing rungs from two caps inside one ladder is the pooling the
`_it{cap}` key suffix exists to prevent. `μ = 0.3` and `μ = 1.0` are **not** re-run —
`μ = 1.0` is already at 1e6, and `μ = 0.3` converged 3/3, which gate M1 below is precisely
the argument for reusing.

## 2. Gates

### M1 — IDENTITY. **PRIMARY, and it gates everything else.**

The cells that **already converged** under the 200000 cap must return **bit-identical**
under 1000000. A cell that stopped because it converged cannot be affected by raising a
limit it never reached.

> **REGISTERED: M1 passes on every previously-converged `mustar` record** (55 of 72 at all
> rungs). If any moves, the cap is not the only thing that changed and **nothing from this
> re-run may be read** — the grader prints M2–M5 but marks them unreadable.

This is the same gate as `wcbidicap`'s C1, which passed 58/58, so it is a repeat of a test
that has succeeded once — that is a reason to expect a pass, not a reason to skip it.

### M2 — RESCUE.

> **REGISTERED: at least 15 of the 17 previously-capped `mustar` cells converge at 1e6.**

Basis, and it is an extrapolation being stated as one: `wcbidicap` rescued **14 of 14**, and
its rescued cells used **205 276 – 475 581** iterations, every one under half the new cap.
`mustarcap`'s capped cells are at *higher friction*, which is the slower direction, so a
miss here is a real outcome and is reported, not re-raised.

### M3 — MONOTONICITY. **The point of the exercise.** Fit-free, no threshold needed.

> **REGISTERED: on the rescued data, `φ(μ)` is non-increasing in `μ_p` at every CV**, in the
> sense that no rise exceeds `PHI_SCALE = 2.1e-3`.
>
> The registered alternative: **a rise above `PHI_SCALE` survives**. That would mean the
> bend is not a cap artefact, and the explanation offered in §1 is wrong — in which case
> the +0.01633 must be explained by something else before any `μ*` is quoted.

⚠ **This gate can fail, and it can fail in a way that refutes the reason for the campaign.**
The currently-violating step is μ 0.6 → 1.0 at CV = 0.5; only the **0.6** end is being
re-run, so the gate is not passable by construction.

### M4 — HOW FAR THE ANSWER MOVES, which is the quantity `wcbidicap` says to expect.

Report, for every re-run point, `φ_new − φ_old` split by whether the cell had converged.

> **REGISTERED, directional and therefore falsifiable: for the previously-capped cells the
> shift is POSITIVE** (under-compressed cells read low), and **zero for the converged ones
> by M1**. A negative mean shift on the capped cells refutes the under-compression reading.

### M5 — `μ*` and `A`, re-fitted. Reported, and **explicitly not the primary gate.**

`F3`'s registered bands (`μ*(0.5)/μ*(0) ∈ [1.05, 1.45]`, `A(0.5)/A(0) ∈ [0.90, 1.20]`) are
**not moved**. They were missed on the mixed-cap data at 2.219 and 1.510.

> ⚠ **`F3a` remains binding: if `A` again leaves [0.90, 1.20], `μ*` and `A` are not
> separated by this design and the verdict is UNDETERMINED whatever `μ*` does.** Repairing
> the cap does not repair an identifiability problem, and this campaign is not registered as
> able to.

`F5`'s form test and the algebraic alternative are re-reported unchanged. `F4`'s three blind
point predictions are re-scored at the **unchanged** ±0.0030 tolerance; two of three missed
on mixed-cap data (−0.02189, −0.03653, both at CV = 0.5). **They are not re-registered** —
a blind prediction is scored once against the protocol it was written for and once against
the repaired one, and both are printed.

## 3. What would make this campaign worthless

- **M1 fails** → the cap is not the only difference between the tags and the whole μ ladder
  needs rebuilding from one campaign, not patched.
- **Fewer than 15 of 17 rescued** → the 1e6 cap is also too low at high friction, and the
  right response is to report the ladder as unmeasurable at this design rather than to raise
  the cap a third time.
- **M3 passes but `A` again moves** → the cap was a real defect *and* the design still
  cannot separate the two mechanisms. This is the **most likely** outcome and is registered
  as such: the campaign would then have repaired a protocol defect without delivering `μ*`.

## 4. What is NOT being run

- No new seeds. `n = 3` is the design and adding seeds after seeing the answer would be
  choosing the design from the result.
- No `μ = 0.3` or `μ = 1.0` re-run (see §1).
- No re-registration of `F3`/`F4`'s bands.
- Nothing at `N = 3000`; this ladder is the `N = 1500` one throughout.
