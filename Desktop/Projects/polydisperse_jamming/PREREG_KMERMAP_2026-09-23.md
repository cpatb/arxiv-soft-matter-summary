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

## RESULT (graded 2026-09-22 23:10 CDT, `grade_kmermap.py` → `out/grade_kmermap.txt`; figure K3)

18/18 cells, provenance 18/18, health 12/18. All six deaths are the L9 jammed-floor failure
(five `STOP dt-collapse`, one `events.h: n < INT_MAX` after Poisson non-convergence); none is a
crash of another kind.

| gate | verdict | numbers |
|---|---|---|
| **R1** doublets jam at x 0.88, singlets cannot | **PASS** | J(w₂ 1) = 0.0579 / 0.1958 at U₀ 2 / 4 (U₀ 1 dead, ≥ 0.0583); J(w₂ 0) ≡ 0 |
| R2 J(doublet) > J(singlet) at every U₀, x 0.92 | **INCOMPLETE** (deaths) | holds on the one healthy pair: 0.636 vs 0.138 at U₀ 4 |
| R3 J non-decreasing in U₀ | **INCOMPLETE** (no row fully healthy) | — |
| **R4** β_max(w₂ 1) ≤ β_max(w₂ 0) + 1/64 | **PASS** 4/4 | |
| R5 health | 12/18 | deaths as expected at x 0.92 w₂ 1 (both), plus x 0.85/0.88 w₂ 1 at low U₀ |

**Post-hoc, labelled, verdicts above unchanged.** A dead cell's recorded max is a LOWER bound.
Every dead cell peaked early (t/τ 0.05–0.93) and was falling when it died. R2′: doublet > singlet is
DECIDED at U₀ 2 (≥ 0.386 vs 0.062) and U₀ 4, UNDECIDABLE at U₀ 1 (both dead). R3′: w₀ x 0.92 and
w₁ x 0.92 are consistent with non-decreasing. **w₁ x 0.88 is NOT**: J(U₀ 1) ≥ 0.0583 against
J(U₀ 2) = 0.0579, a 0.0004 inversion (0.7 %). The honest reading is flat between U₀ 1 and 2. Whether
0.0004 is resolvable at L9 is unmeasured.

**What the map says (reading, not a registered gate).** At the clean realistic-friction inputs
(χ, σ\*, δ all from the μ_p 0.3 DEM), the doublet fraction switches a fixed-φ drop from a regime
that cannot jam (x 0.88 is outside the singlets' window) to one that jams at U₀ 2 and 4. This is
the first DNS result for the multiplet hypothesis that does not use a predicted or excluded χ.

The SPREAD tells a different story from the jammed volume. The doublets' cut in max Rg_r is
largest where NOTHING can jam (x 0.85, U₀ 4: −31 %; healthy pair), smaller where only the doublets
jam (x 0.88: −14 % / −18 %), and ~−3 % where both jam (x 0.92). The footprint is set mostly by
thickening outside the window, not by jamming. **A reduced spread is therefore not evidence of
jamming**, and a bench test that reads arrest from the footprint cannot tell the two apart.
