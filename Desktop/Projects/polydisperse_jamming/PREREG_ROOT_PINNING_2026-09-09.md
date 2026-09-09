# Pre-registration — is `phi_c` a material root, or the bottom of the ladder?

**Written 2026-09-09, BEFORE the deciding cells were launched.** The launch command and
the grading script are named below; neither existed before this file.

---

## 1. What prompted it

`analyse_phase8._pressure_root` (estimator E2, the PRIMARY one) fits
`p_ss = c (phi - phi_c)^b` and is bounded above by `lo - 1e-6`, where `lo` is the lowest
rung it is given. Measured across **every** group of `phase8_prod2*` and
`phase8_p5uniform*` under all four gate sets — 46 resolved fits — the returned root sits

    0.0003 to 0.0124 below the lowest admitted rung, and NEVER more.

The ladder spacing is 0.020, so every root lands inside the first interval below the
ladder. Worse, the quality gates supply a natural experiment: at `cv` = 0 and 0.2 the
quasi-static gate removes exactly the bottom rung, moving `lo` up by one spacing. The
root moves with it:

| cell | `lo` shift | root shift |
|---|---|---|
| lognormal cv 0.00 seed 1 | +0.0200 | +0.0269 |
| lognormal cv 0.00 seed 2 | +0.0200 | +0.0250 |
| lognormal cv 0.00 seed 3 | +0.0200 | +0.0168 |
| lognormal cv 0.50 seed 2 | +0.0200 | +0.0150 |
| uniform  cv 0.50 seed 1  | +0.0200 | +0.0144 |
| uniform  cv 0.50 seed 2  | +0.0200 | +0.0144 |

Mean response **+0.0188 of root per +0.0200 of ladder**, i.e. a slope of 0.94. A root
that follows its own fitting window one-for-one is not obviously a property of the
material.

And the ladders are centred on the project's **own published isotropic `phi_J`**
(`phase8_criticalstate.PHI_J_ISO`), so if the root is pinned to the ladder then

    delta phi_c(CV)  ==  delta phi_J_iso(CV)   by construction.

Measured: `delta phi_J_iso` = 0.06232 (lognormal 0 -> 0.5) and 0.02278 (uniform), against
reported `delta phi_c` = 0.0623 and 0.0215. **100.0 % and 94.3 %.** So the claim that the
critical-state slope 0.1247 "matches the independently measured isotropic slope 0.1246"
may be two readings of one number, not corroboration — the second sharing the first's
assumption is the failure mode this project has already met twice.

⚠ This does NOT touch **P3**, the matched-`phi` pressure collapse (5.4x-43x), which uses
no fit, no root and no ladder placement.

## 2. The test

A ladder at `cv` = 0, `mu_p` = 0.3, seeds 1-3, offsets **-0.06, -0.04, -0.02, 0.00,
+0.02** from the same centre 0.574533 — i.e. rungs at phi = 0.5145, 0.5345, 0.5545,
0.5745, 0.5945. Same family, N, `gamma_max`, `shear_steps` and protocol as prod2.
Tag `downlad`, so it pools with nothing (`phase8_downlad*.jsonl`).

The up-ladder's answer for this cell is `phi_c` = **0.5695 +- 0.0030** (ungated E2,
n = 3 seeds).

## 3. The two hypotheses, and what each predicts

**H_pinned — the root is the bottom of the window.** E2 on the down-ladder returns a
root within one spacing below the new bottom rung:

    phi_c(down)  in  [0.4945, 0.5145]     i.e. a shift of about -0.055

**H_real — 0.5695 is a material critical state.** Then the two lowest rungs are BELOW
it and cannot sustain a steady shear stress, so, independent of any fit:

    p_ss_hat(phi = 0.5145) and p_ss_hat(phi = 0.5345)  < 1e-6   (prod2's jammed rungs
                                                                run 2e-4 to 2e-2)
    z_nr_ss(phi = 0.5145) < 4     (the frictional isostatic count)
    and E2 returns phi_c(down) in [0.565, 0.575], a shift of < 0.005

**Separation.** The two predictions for `phi_c(down)` are 0.055 apart: **18 seed
standard deviations** and 2.75 ladder spacings. The test is informative.

**Neither passes.** If the low rungs are unjammed AND the root still moves by more than
0.02, then `phi_c` is bracketed physically but E2 is not the estimator that measures it,
and the pressure-root family must be replaced (E4 is stress-free but returns NaN at
`cv` = 0.5, so it cannot span the contrast).

## 4. What each outcome does to the shipped claims

| outcome | `delta phi_c` = +0.044 to +0.062 | slope "matches" 0.1246 | P5 2.90x | P3 5.4-43x |
|---|---|---|---|---|
| H_real | stands | stands, and is corroboration | stands | stands |
| H_pinned | **withdrawn**: it is `delta phi_J_iso` restated | **withdrawn**: one number, twice | **withdrawn**: 2.90x is the isotropic 2.73x | **stands** |

## 5. Commands

    bash scripts/launch_0909b.sh            # via WSL; 15 cells
    python grade_root_pinning.py            # grades this file's predictions

---

## 6. CORRECTIONS TO §1, made after the run, by a referee

**The registered predictions in §3 and the outcome table in §4 are untouched.** What is
corrected here is §1's *motivating* arithmetic, which was computed before the +0.10 and
+0.12 rungs landed the same day and no longer reproduces. A referee caught it.

- **"0.0003 to 0.0124 below the lowest rung, across 46 resolved fits, and NEVER more"
  is wrong on the current data.** Recomputed by `grade_root_pinning.py`, which now
  prints the census every run so it cannot go stale again: **57 resolved fits, offsets
  0.0002 to 0.0323, median 0.0043; 7 of 57 exceed one ladder spacing and 0 exceed two.**
  Every offset beyond the old ceiling is on a 7-rung ladder under the `steady` or `both`
  gate, where the gate drops the bottom **two** rungs rather than one — so §1's other
  claim, "the quasi-static gate removes exactly and only the bottom rung", is true of
  the quasi-static gate on the 5-rung ladders and not of the other gates on the 7-rung
  ones.
- The census was never the load-bearing evidence and it should not have been written as
  if it were. The load-bearing evidence is the **response**: +1.15 of the ladder move
  when the bottom moves, +0.06 when the top moves, and the down-ladder result below.

## 7. OUTCOME, graded 2026-09-09 — `out/root_pinning_verdict.txt`

**SPLIT, the third branch §3 named.** Neither hypothesis's fitted band was hit.

| | |
|---|---|
| H_pinned, root in [0.4945, 0.5145] | **fail** |
| H_real, root in [0.565, 0.575] | **fail** |
| H_real, fit-free: the rungs below the claimed root cannot carry stress | **PASS** |

- **Fit-free (3 seeds):** at φ = 0.5145 / 0.5345 / 0.5545 / 0.5745 / 0.5945 the steady
  grain pressure is 1.16e-07 / 2.72e-07 / 7.65e-07 / 3.74e-06 / 4.73e-04 and `z_nr` is
  2.38 / 2.84 / 3.30 / 3.94 / 4.58, with the rattler fraction falling 0.796 → 0.066.
  The two rungs below the up-ladder's claimed root 0.5695 are unjammed on **both**
  channels. **The critical state is genuinely bracketed, between 0.5545 and 0.5945.**
- **E2 refuses on all three seeds**, with its exponent saturated at the bound `b` = 8.00.
  The roots it would have returned are 0.4511 / 0.4522 / 0.4594 — far *below* the
  ladder, not merely under its bottom rung. This is a **third failure mode**: a power law
  cannot represent a rung below its own root, where `p → 0` and `log p → −∞`, so the fit
  runs the exponent to its ceiling instead of moving the root. E2 therefore produces a
  number **only when the whole ladder sits above the root**, and in that case the number
  it produces tracks the bottom rung.
- **E4 gives 0.5763 ± 0.0013 on the up-ladder and 0.5763 ± 0.0013 on the down-ladder** —
  identical across a 0.060 move of the window. ⚠ That robustness is *structural*, not a
  surprise: E4 interpolates between the two rungs that bracket the crossing and both
  ladders contain them. It is the right property for the quantity, and it is what E2
  lacks; it is not an independent second measurement.

### What this does to the shipped claims

| claim | status |
|---|---|
| Δ`φ_c` = +0.044 … +0.061, sign and rough size | **stands** — but as a statement about a ladder-relative root, not a measured material constant |
| "the critical-state slope 0.1247 **matches the independently measured** isotropic 0.1246" | 🔴 **withdrawn as corroboration** — the ladders are centred on those very `φ_J` values and the root tracks the ladder |
| `φ_c` ≈ `φ_J` to within 0.005 at `μ_p` = 0.3 | 🔴 **withdrawn** — forced by the ladder placement |
| P5's 2.82× family separation | **stands as a measurement**, but see `p5_degeneracy.py`: the same separation is present in `φ_J`, with no shear in it |
| P3, the matched-`φ` pressure collapse 5.4×–43× | **untouched** — no fit, no root, no ladder placement |
| P4, the frictionless anchor | **untouched** — E4 there is a genuine bracketed interpolation and `Z` = 2d = 6 is exact |
