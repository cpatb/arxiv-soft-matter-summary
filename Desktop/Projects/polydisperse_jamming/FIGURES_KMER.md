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

## K4 — apex recoil against jammed volume (PREREG_JAMSIG)
Every healthy default-cap DNS cell: the maximum rise of the apex height h/D0 after its minimum,
against the jammed volume fraction J. Hollow: the 23 cells the observable was found on (in-sample).
Filled: 5 never-seen singlet cells from Amendment 2 (out-of-sample). Dashed: the threshold 0.1, fixed
before those cells ran. Non-jammed drops re-rise 0.125–0.59; jammed drops at most 0.078.
⚠ h/D0 is quantised at 1/64, so 0 means "no rise of one quantum". The out-of-sample JAMMED side is
one cell (3 of 4 died at L9). The contrast is set by the jammed-state viscosity: it holds at a 10×
stiffer cap and falls to a quarter at a 10× softer one (`PREREG_EPSSIG`). A candidate bench
observable, not a validated one.

## K5 — the regime map at EPS_MIN 0.03 (PREREG_KMERMAPE03)
Same layout and inputs as K3, with a 10× softer viscosity cap. All 18 cells survive, which completes
the orderings the default-cap deaths left open (R2′, R3′ PASS), and every cell keeps its J+/J0 class
(X1 12/12). J magnitudes are NOT comparable with K3 (they depend on the cap); classes and orderings
are. Singlet cells at x 0.85/0.88 are byte-identical to K3's (the cap cannot bind there).
