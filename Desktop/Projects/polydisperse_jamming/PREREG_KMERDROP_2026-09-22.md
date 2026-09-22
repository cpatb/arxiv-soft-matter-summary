# PREREG — `kmerdrop`: what the doublet fraction does to an impacting drop, in the DNS

Written 2026-09-22 (fourth block), before launch. Theory: `derivations/multiplet_knobs_2026-09-22.md`.
Launcher: `scripts/kmerdrop_2026-09-22.sh` (reuses the G0-gated `jamwin` binary `db413e87dd21`).

## 1. Design

The DNS knows particles only through `(x, χ, σ*, δ)`. Here ONLY `χ` moves, taking the five
DEM-measured / H_mix values of the singlet→doublet ladder (μ_m = 3 convention, p̂ 3e-5):
χ = 0.886 / 0.859 / 0.834 / 0.811 / 0.788 at w₂ = 0 / 0.25 / 0.5 / 0.75 / 1, with `PHI0` = 0.62,
`PHIM` = χ·0.62, `PHI` = 0.527 (x = 0.85, the jamwin arm-B position), σ\* = 300 Pa, δ = 1.3169,
everything else as `PREREG_JAMWIN` §1 (cs_fit base, U₀ = 2, L9, WC_CLOCK = 1). `JAM_FHAT = f_c(w)`
(1.0 where x < χ, so J ≡ 0 there by construction and is not reported as a jam).

⚠ Holding σ\* and δ fixed isolates the WINDOW channel. The DEM predicts doublets also LOWER σ\*
(and δ); both would push further the same way, so this ladder UNDER-states the doublet effect.

## 2. Predictions

- **D1** `β_max` falls monotonically with w₂ across all five (continuous thickening Θ rises as χ
  falls even where no jamming is possible — K5).
- **D2** `J(f_c)` is 0 at w₂ = 0, 0.25 (by construction) and **> 0 and non-decreasing** over
  w₂ = 0.5, 0.75, 1 — the drop starts to jam once the window has opened past x = 0.85.
- **D3** health: every cell reaches t/τ ≥ 5.9 with rc = 0; a death is printed, never read as a jam.
  (jamwin: arm-B cells at σ\* = 300 survived; deeper jams died.)

## 3. What this cannot show

That doublets move χ — that is the DEM's (`PREREG_KMER`); this takes it as input.

---

## RESULT (graded 2026-09-22 17:31 CDT, `out/grade_kmerdrop.txt`)

D3 health 5/5 (all rc 0, t/τ = 6, provenance 5/5).

- **D2 PASS.** Jammed volume `J(f̂ > f_c)` = 0 / 0 / **0.049 / 0.220 / 0.385** at w₂ = 0 / 0.25 /
  0.5 / 0.75 / 1: the drop at x = 0.85 does not jam as singlets and does jam, increasingly, once
  the doublet fraction has opened the window past x.
- **D1 FAIL as registered** (strict monotone β_max): 1.391 / 1.172 / 1.109 / **1.078 / 1.078** —
  the last two TIE. β = 2·r_max is grid-quantised (one quantum = 1/64 ≈ 0.016 at L9; lessons.md
  2026-08-29), so the pair differs by less than one quantum. The FAIL stands. Printed beside it,
  NOT graded: the continuous radial radius of gyration at its maximum falls strictly,
  0.4596 / 0.3875 / 0.3631 / 0.3534 / 0.3477.
- **Magnitude:** from singlets to pure doublets at the same weighed φ, the maximum spread falls
  from β 1.39 to 1.08 (−22 %) through the WINDOW channel alone — σ\* and δ were held fixed, and
  the DEM predicts doublets also lower σ\*, so this under-states the doublet effect on that channel.
