# PREREG — `epssig`: is the apex-recoil signature a readout of the viscosity regularisation?

Written 2026-09-23 ~03:05 CDT, before any cell. Launcher `scripts/epssig_2026-09-23.sh`, grader
`grade_epssig.py` (both written before the data). Same binary `db413e87dd21`, cs_fit, L9, WC_CLOCK 1,
θ₀ 5°, TEND 6, CFL 0.10, MP_EXPO 2; x 0.85 (φ 0.527, φ₀ 0.62), U₀ 2, σ\* 300 Pa, δ 1.3169 — the
`kmerdrop` ladder's settings.

## Why (solver source read before registering)
`wc_drop_impact.c`: η = η_s (1 − φ/φ_J(f̂))^(−2), with ε = 1 − φ/φ_J floored at EPS_MIN (default
0.01, η_max = 10⁴ η_s). A JAMMED region has ε ≤ 0, so its viscosity IS the cap η_max. A fully
activated but UNjammed region sits at ε = 1 − φ/φ_m:
- w₂ 0.25: 1 − 0.527/0.53258 = **0.0105**, 1.05× the floor, i.e. 0.91 η_max;
- w₂ 0.75: at φ_m 0.50282 it is jammed, ε floored.

So at the default EPS_MIN the "jammed vs thickened" viscosity contrast is small and set by the
floor. The recoil signature S2 (`PREREG_JAMSIG`: J0 ≥ 0.125, J+ ≤ 0.078) could be a readout of
EPS_MIN — the knob-in-both-model-and-detector trap of lessons.md.

## Cells (4)
| label | w₂ | φ_m | JAM_FHAT | EPS_MIN | η_max / η_s |
|---|---|---|---|---|---|
| e003_w0.25 | 0.25 | 0.53258 | 1.0 | 0.003 | 1.1e5 |
| e003_w0.75 | 0.75 | 0.50282 | 0.7937 | 0.003 | 1.1e5 |
| e03_w0.25 | 0.25 | 0.53258 | 1.0 | 0.03 | 1.1e3 |
| e03_w0.75 | 0.75 | 0.50282 | 0.7937 | 0.03 | 1.1e3 |
Default-EPS_MIN references: `kmerdrop` w0.25 (S2 0.172, J0) and w0.75 (S2 0.047, J+).

## Predictions (two-sided: live below, dead above)
- **G1 (the signature is not a floor artefact below the default):** at EPS_MIN 0.003, where the cap is
  10× higher and the thickened w₂ 0.25 sits far below it (ε 0.0105 vs 0.003), S2(w0.75) < 0.1 ≤
  S2(w0.25).
- **G2 (it IS carried by the contrast):** at EPS_MIN 0.03, BOTH arms' high-activation regions are
  clamped at the same η_max (w₂ 0.25's ε 0.0105 < 0.03), so the contrast is removed. Prediction:
  the S2 gap |S2(w0.25) − S2(w0.75)| falls below HALF its default value (0.125/2 = 0.0625).
- **Reading of the four outcomes, fixed now:** G1 ✓ G2 ✓ → S2 is a signature of the jammed-vs-
  thickened VISCOSITY CONTRAST, real in the model and set by how stiff a jammed region is. G1 ✗ →
  S2 is not robust to the cap and is withdrawn as a candidate. G1 ✓ G2 ✗ → S2 does not depend on the
  viscosity contrast; the mechanism is something else and is unexplained.
- **Health:** rc 0 and t/τ ≥ 5.9, and the banner's `eps_min=` must equal the cell's EPS_MIN; else
  EXCLUDED. A dead cell's S2 is a partial value and decides nothing. Deaths are expected at 0.003.
- J is reported per cell (its definition, f̂ > f_c, does not involve the cap, but the dynamics do).
