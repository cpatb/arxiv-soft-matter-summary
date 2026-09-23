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
