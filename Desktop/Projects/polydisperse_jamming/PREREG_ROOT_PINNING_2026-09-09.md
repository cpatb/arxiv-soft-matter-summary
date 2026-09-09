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
