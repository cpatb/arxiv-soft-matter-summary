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

## RESULT (graded 2026-09-23 ~03:20 CDT, `grade_epssig.py` → `out/grade_epssig.txt`; 4/4 healthy, eps tags verified)
| EPS_MIN (η_max/η_s) | w₂ 0.25 (J0): J, S2, Rg | w₂ 0.75 (J+): J, S2, Rg | S2 gap |
|---|---|---|---|
| 0.003 (1.1e5) | 0, 0.1719, 0.3875 | 0.0962, **0.0469**, 0.3523 | 0.125 |
| 0.01 (1e4, default) | 0, 0.1719, 0.3875 | 0.2197, 0.0469, 0.3534 | 0.125 |
| 0.03 (1.1e3) | 0, 0.1719, 0.3947 | 0.5694, **0.1406**, 0.3776 | 0.031 |

**G1 PASS, G2 PASS → the registered reading: S2 is a signature of the jammed-vs-thickened VISCOSITY
CONTRAST.** At a jammed-region viscosity of 10⁴ or 1.1e5 η_s the jammed drop does not recoil
(S2 0.047 at both, identical); at 1.1e3 η_s it recoils almost as much as the unjammed one (0.141 vs
0.172). The stiffness at which recoil switches off is bracketed in (1.1e3, 1e4) η_s, not measured.

- **Control, exact:** e003_w0.25 is BYTE-IDENTICAL to the default w₂ 0.25 cell (spread.dat md5
  30e112861cce), while its banner confirms eps_min 0.003 / η_max 60.8 code. That is the algebra of
  the prereg: this cell's ε never falls below 0.0105, so a floor moved from 0.01 to 0.003 cannot
  bind. The J0 arm is therefore cap-independent below the default. In the jammed arm, S2 is also
  unchanged at the stiffer cap and only J moves (below).
- ⚠ **New and important: J's MAGNITUDE depends on the cap.** w₂ 0.75 gives J = 0.096 / 0.220 / 0.569
  at EPS_MIN 0.003 / 0.01 / 0.03 (×5.9 over one decade). A stiffer jammed region deforms less,
  shears less and activates less. The J+/J0 CLASS held at all three caps for this pair. **Every
  quoted J value in this project's DNS (kmerdrop, kmermap, kmerdecomp) is a default-EPS_MIN number
  whose magnitude is regularisation-dependent; the classes and orderings are what has been
  tested.** The x 0.88 kmermap J (0.058 at U₀ 2) is the one most at risk of flipping class at a
  stiffer cap, and is untested.
- Reading for the bench (REASONING, not measured): a real jammed suspension is a yield-stress solid,
  i.e. stiffer than any finite cap, so the recoil contrast should be at least as large on the bench
  as at the default cap. That is the direction G1 tested.

## Amendment 1 — does the HEADLINE (kmermap R1) survive the cap? (written 2026-09-23 ~03:30 CDT, before its cells)
Four cells re-run `kmermap` w1_x0.88 (pure doublets, clean μ_p 0.3 inputs: φ 0.5456, φ_m 0.53512,
σ\* 235.8 Pa, δ 1.1163, JAM_FHAT 0.8766) at U₀ 2 and 4 with EPS_MIN 0.003 and 0.03:
`r1_e003_U2.0`, `r1_e003_U4.0`, `r1_e03_U2.0`, `r1_e03_U4.0`.
- **H1:** R1 survives the cap: at EACH EPS_MIN, J > 0.01 in at least one healthy U₀ cell (the
  default gave 0.058 / 0.196). J is expected LOWER at 0.003 (the w₂ 0.75 pattern), so H1 at 0.003 is
  the risky half.
- A dead cell's J is a lower bound and counts toward H1 only if it exceeds 0.01, as in kmermap's
  post-hoc reading (declared here, in advance, for this amendment).
- J reported per cell beside its default value; singlets at x 0.88 remain J ≡ 0 by construction.
