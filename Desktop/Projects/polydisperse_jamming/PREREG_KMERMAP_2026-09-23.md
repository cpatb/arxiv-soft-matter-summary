# PREREG — `kmermap` (A3): the drop-impact regime map in (x, doublet fraction, U₀), all channels

Written 2026-09-22 ~21:50 CDT, before any cell. Launcher `scripts/kmermap_2026-09-23.sh` (reuses the
G0-gated jamwin binary `db413e87dd21`). Solver `wc_drop_impact.c`, cs_fit base, L9, WC_CLOCK 1,
θ₀ 5°, TEND 6, MP_EXPO 2, DIL_COEFF 0, PHI0 0.62.

## Inputs, all from this project's DEM (no bench data)
| | singlets (w₂ 0) | pure doublets (w₂ 1) | source |
|---|---|---|---|
| χ | 0.9083 | 0.8631 | `kmer`, μ_m 0.3, p̂ 3e-5, clean |
| σ\* | 300 Pa (reference) | 235.8 Pa (×0.786) | `shearkmer`, sheared, μ_p 0.3 |
| δ | 1.3169 (reference) | 1.1163 (×0.848) | `shearkmer` |

Grid: x ∈ {0.85, 0.88, 0.92} × w₂ ∈ {0, 1} × U₀ ∈ {1, 2, 4} m/s = 18 cells. `JAM_FHAT` = f_c where the
suspension is inside its window (w₀ x 0.92: 0.8724; w₁ x 0.88: 0.8766; w₁ x 0.92: 0.5844), 1.0 where
it is outside (J ≡ 0 there by construction, not reported as a jam).

The three x rows are the three regimes the theory names: x = 0.85 is OUTSIDE both windows (neither
can jam; only thickening); x = 0.88 is inside the doublets' window only; x = 0.92 is inside both.

## Predictions
- **R1** x = 0.88: doublets jam at ≥ 1 U₀ (`J > 0.01`); singlets cannot (by construction). The
  doublet fraction switches the regime at fixed φ and fixed impact.
- **R2** x = 0.92: `J(w₂ 1) > J(w₂ 0)` at EVERY U₀ (f_c 0.58 vs 0.87).
- **R3** inside a window, `J` is non-decreasing in U₀ (more impact stress, more activation).
- **R4** `β_max(w₂ 1) ≤ β_max(w₂ 0)` at every (x, U₀) — to within one grid quantum (1/64).
- **R5** health: rc 0 and t/τ ≥ 5.9, else EXCLUDED and printed (jamwin: deep jams die at L9 —
  expected at x 0.92, w₂ 1, where f_c is lowest). A death is never a jam.
