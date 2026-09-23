# Singlet/doublet figures (2026-09-22, fourth block)

Built by `python make_kmer_figures.py` (`--check` rebuilds K1 and compares md5).

## K1 — the doublet window (PREREG_KMER)
(a) χ = φ(μ_m)/φ(μ 0) against doublet volume fraction w₂ at p̂ 3e-5, seed-paired, for μ_m = 3
and 0.3; dashed = ideal specific-volume mixing of the campaign's own endpoints. (b) σ\*(w₂)/σ\*(0)
(squares, left) and δ₁₀/₉₀ (triangles, right) at μ_p 0.3 from the per-contact dumps, where every
contact is primary–primary (R\* = 1/4); shaded = the registered σ\* band [0.55, 0.95].
⚠ Isotropic packs, not sheared.

## K2 — the drop (PREREG_KMERDROP)
🔴 Scope (referee pass 3): the χ ladder here is the μ_m 3 PREDICTION (0.886 → 0.788). At the clean
μ_p 0.3 χ (0.9083 → 0.8631), x = 0.85 is outside the window at every w₂, so this panel shows what
the solver does IF χ moves that far. It is not a realistic-friction result; that test is `kmermap`.

The Wyart–Cates DNS at a FIXED x = 0.85 with χ taken from the doublet ladder (σ\*, δ fixed):
jammed volume fraction (fraction with f̂ > f_c) and the maximum radial radius of gyration.
Spreading falls most between w₂ = 0 and 0.25, where NO jamming is possible (x < χ): that part is
continuous thickening (Θ rises as χ falls). Jamming switches on once χ < x (w₂ ≳ 0.33 by H_mix).
β_max is grid-quantised at L9 and ties at w₂ 0.75/1; the radius of gyration is continuous.

⚠ In K1(a) the whole μ_m = 3 series is drawn hollow as K4-flagged, although its w₂ = 0 point
(singlets) passes K4; only the w₂ > 0 points are excluded.
