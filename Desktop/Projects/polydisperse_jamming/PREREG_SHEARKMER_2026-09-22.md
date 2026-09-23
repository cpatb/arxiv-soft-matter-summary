# PREREG — `shearkmer`: the doublet onset and sharpness in SHEARED packs

Written 2026-09-22 (fourth block, late), before launch. Why: `kmer`'s H_σ\* PASS (σ\*(doublets)/
σ\*(singlets) = 0.783) and H_δ result (δ unchanged) come from ISOTROPIC packs; Wyart–Cates
activation is about SHEARED contacts — the exact caveat the CV lead carried until `actshape`.

## Grid
Constant-pressure shear, `actshape`'s protocol exactly (rate-matched CV 0: 102 000 steps, gain
0.15, γ_max 3, φ_start 0.55, p̂ 3e-4, N = 1000 GRAINS), dumps on, tag `shearkmer`:
w₂ ∈ {0, 0.5, 1} (tangent doublets of one primary, `--doublet-w`) × μ_p ∈ {0, 0.3} × seeds {1,2,3}
= **18 cells**.

## Gates
- **SK0 identity:** the w₂ = 0 cells ARE `actshape`'s CV-0, p̂ 3e-4 cells (same key); `p_hat_series`,
  `z_series` and `force` must be bit-identical to `phase9_actshape_c00a` records.
- **SK1 health:** steady, servo_held, quasistatic (actshape's W4).
- **SK2 overlap (K4 of PREREG_KMER):** > 0.1 % of contacts beyond u/r 0.05 → excluded.

## Predictions (from the isotropic `kmer` result, before any sheared cell)
- **S1** σ\*(w₂ = 1)/σ\*(0) at μ_p 0.3 ∈ **[0.55, 0.95]** (R\* = 1/4 for every contact, σ = p).
- **S2** |δ₁₀/₉₀(w₂ = 1) − δ₁₀/₉₀(0)| < **0.05** at μ_p 0.3 (sharpness unchanged, as isotropic).
- **S3** σ\*(w₂ = 0.5) lies between the two endpoints.
Seed-paired, t_crit(df 2) = 4.303 where a contrast is tested.
⚠ Z from sheared clump packs is member-level (shear.py's contact graph is not grain-aware) and is
not read.

## Launch 1 INVALID (recorded, not deleted)
The first launch built SINGLETS for every w₂ > 0 cell: `shear._prepare` constructs its own
`PackingState` for insertion and did not carry `doublet_w` into it. Caught because σ\*/κ and δ
came back IDENTICAL across w₂ = 0 / 0.5 / 1 to every printed digit — byte-identity between cells
that should differ (lessons.md 2026-08-31). The 12 w₂ > 0 records and dumps are kept in
`out/invalid_shearkmer_v1/`. Fixed in `shear.py`; the grader gains **SK3 WIRING** (realised doublet
count must match w₂), which would have failed launch 1 on every cell. The w₂ = 0 cells are valid
(SK0 6/6 bit-identical) and are not re-run. Predictions S1–S3 unchanged.

## RESULT (graded 2026-09-22 ~21:00 CDT, `out/grade_shearkmer.txt`, launch 3, 18/18)

Gates: SK0 6/6 bit-identical, SK3 wiring 12/12, SK1 health 18/18, SK2 overlap 0 excluded
(sheared μ 0.3 doublet packs do NOT interpenetrate, unlike isotropic μ 3 clump packs).

- **S1 PASS:** σ\*(pure doublets)/σ\*(singlets) = **0.786** at μ_p 0.3 — the isotropic 0.783 reproduced
  under shear. (μ_p 0: 0.853, reported.)
- **S2 FAIL:** under shear, frictional doublets DO soften the activation: δ₁₀/₉₀ 1.254 → **1.063**
  (−0.191 ± 0.008, seed-paired). The isotropic "δ unchanged" does not transfer to sheared packs.
  (μ_p 0: −0.002 ± 0.013, unchanged.)
- **S3 PASS:** σ\*(w₂ 0.5) lies between the endpoints (0.171 in [0.153, 0.194]).

**Consequence, computed without a Weibull assumption** (`shearkmer_sjam.py`; reported, not
registered): at fixed chemistry and primary, pure doublets change the stress to jam by
×0.67 / 0.74 / 0.79 / 0.83 / 0.94 / **1.05** at f_c = 0.10 / 0.25 / 0.37 / 0.50 / 0.75 / 0.90 (μ_p 0.3).
The δ-channel offset the referee flagged is REAL and appears only at the window edge, where it
slightly reverses the sign. At fixed φ doublets also LOWER f_c (the window opens), moving the
operating point away from that edge — so the net effect at fixed weighed φ is toward jamming.
Frictionless: a uniform ×0.85, as a pure onset-scale shift would give.
