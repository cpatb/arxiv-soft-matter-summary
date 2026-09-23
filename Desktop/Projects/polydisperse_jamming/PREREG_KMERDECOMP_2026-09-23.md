# PREREG — `kmerdecomp`: which channel does the doublets' work in the DNS drop?

Written 2026-09-22 ~23:55 CDT, before any cell. Launcher `scripts/kmerdecomp_2026-09-23.sh` (clone of
the kmermap launcher: same G0-gated binary `db413e87dd21`, cs_fit, L9, WC_CLOCK 1, θ₀ 5°, TEND 6,
CFL 0.10, MP_EXPO 2, PHI0 0.62). Grader `grade_kmerdecomp.py`, written before the data.

## Why
`kmermap` moved all three doublet channels together (χ 0.9083 → 0.8631, σ\* 300 → 235.8 Pa, δ 1.3169 →
1.1163), so no row isolates one (its referee, MINOR 3). Two readings depend on the split: "the
doublet fraction switches the regime at x 0.88" (is the WINDOW enough?) and "the largest footprint
cut is at x 0.85 where nothing jams" (is that cut THICKENING, i.e. the window, or the lower onset?).

## Cells (4). Reference arms are the existing kmermap cells, same binary and settings.
| label | x | φ | χ (φ_m) | σ\* | δ | JAM_FHAT | U₀ | isolates |
|---|---|---|---|---|---|---|---|---|
| win_x0.88_U2.0 | 0.88 | 0.5456 | 0.8631 (0.53512) | 300 | 1.3169 | 0.8766 | 2 | window only |
| win_x0.88_U4.0 | 0.88 | 0.5456 | 0.8631 (0.53512) | 300 | 1.3169 | 0.8766 | 4 | window only |
| win_x0.85_U4.0 | 0.85 | 0.527 | 0.8631 (0.53512) | 300 | 1.3169 | 1.0 (outside) | 4 | window only |
| act_x0.85_U4.0 | 0.85 | 0.527 | 0.9083 (0.56315) | 235.8 | 1.1163 | 1.0 (outside) | 4 | activation only |

## Predictions
- **D1** the window alone switches the regime: `J > 0.01` in at least one of win_x0.88_U2/U4
  (singlets there: J ≡ 0).
- **D2** (no direction predicted) `J(win) − J(all-channel)` at x 0.88, U₀ 2 and 4, reported with sign.
  Reason for no direction: at f_c 0.88 > 1/e the softer δ RAISES σ_jam while the lower σ\* lowers it,
  and the sheared DEM net there is ×1.05 (`out/shearkmer_sjam.txt`) — the two nearly cancel.
- **D3** thickening carries the majority of the x 0.85 footprint cut: with R = max Rg_r at U₀ 4,
  `ln(R_win/R_w0) < ln(R_act/R_w0)` (the window cut is the larger), where R_w0 = 0.5742 and the
  all-channel R = 0.3963 (kmermap). Reason: Θ rises ×17.9 from the window while at U₀ 4 the activation
  is already near saturation (⟨f̂⟩max 0.88 for singlets), so moving σ\* by 21 % should matter less.
- **D4** (reported, no bar) additivity: ln(R_all/R_w0) against ln(R_win/R_w0) + ln(R_act/R_w0).
- **D5** health: rc 0 and t/τ ≥ 5.9, else EXCLUDED and printed; a death is never a jam, and a dead
  cell's maxima are lower bounds, reported as such.
- **Provenance:** each banner must carry the table's φ, φ_m, σ\*, δ, JAM_FHAT, U₀ to 1e-6.

## RESULT (graded 2026-09-23 ~00:30 CDT, `grade_kmerdecomp.py` → `out/grade_kmerdecomp.txt`)

4/4 cells, provenance 4/4, health 3/4 (win_x0.88_U4.0: `STOP dt-collapse` at t/τ 2.6, the L9
jammed floor; its J 0.2009 is a lower bound).

| gate | verdict | numbers |
|---|---|---|
| **D1** window alone jams at x 0.88 | **PASS** (healthy cell) | J 0.0591 at U₀ 2; ≥ 0.2009 at U₀ 4 (dead) |
| D2 J(window) − J(all), x 0.88 | reported | **+0.0011** at U₀ 2 (0.0591 vs 0.0579); U₀ 4 undecidable (dead) |
| **D3** window cut > activation cut, x 0.85 U₀ 4 | **PASS** | ln(R/R_w0): window −0.3721, activation **−0.0008**, all −0.3709 |
| D4 additivity | reported | ln(all) − [ln(win) + ln(act)] = +0.0021 |
| D5 health | 3/4 | |

**Reading.** In the DNS drop the doublets' effect is carried **entirely by the window channel (χ)**
at the conditions tested. Swapping in the doublets' lower σ\* and softer δ with the singlets' χ moves
the footprint by 0.08 % (x 0.85). Swapping in the doublets' χ alone reproduces both the jamming at
x 0.88 (J to 0.001) and the whole footprint cut at x 0.85 (to 0.2 %). The channels add in log Rg to
0.002. This agrees with the sheared-DEM expectation that onset and sharpness nearly cancel near the
window edge (net ×0.94–1.05 at f_c 0.75–0.9), and extends it: here the activation channel is inert
even outside the window.

**Scope.** Tested at U₀ 2 and 4 only, where the singlets' activation is already high (⟨f̂⟩max
0.77–0.92). At low U₀, where activation is far from saturated, a 21 % lower σ\* could matter and is
untested. Stress-controlled activation (M6), L9, one seed-free DNS realisation per cell.
