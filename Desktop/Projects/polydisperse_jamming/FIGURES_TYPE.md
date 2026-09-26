# Figure T1 — the thickening-type map (`make_type_figures.py`, `figures/T1_thickening_type.{pdf,png,csv}`)

Regenerate: `python make_type_figures.py --check` (PNG byte-identical on regeneration).

**(a) Where a Wyart–Cates suspension is CST, DST or shear-jamming**, in x = φ/φ₀ against χ = φ_m/φ₀ (closed form (T),
`derivations/thickening_type_2026-09-25.md`). Solid black: x = χ, the lower edge of shear jamming (SJ above it).
Broken curves: the DST onset x_DST(χ, nδ) for nδ = 2.147 (blue dashed; the CV 0.5 activation with n = 2), 2.634
(orange dash-dot; CV 0 activation, n = 2) and 6.58 (green dotted; CV 0 activation with the cornstarch fit's n = 5).
Between a broken curve and the black line the steady flow curve is S-shaped (DST); below it thickening is continuous
(CST). The DST band shrinks to nothing as χ → 1 and widens with nδ; it exists only for nδ > 1.

**(b) The 24 cvmap DNS drops** (EPS_MIN 0.03, `out/grade_cvmapE03.txt`; POST HOC classification
`derivations/cvmap_by_type_2026-09-26.txt`): peak spread Rg_max/D₀ against each drop's distance from ITS OWN DST
boundary, x − x_DST(χ, nδ) with its own χ and δ (φ₀ 0.62, n 2). Squares: SJ (x ≥ χ); diamonds: DST; circles: CST.
Every drop on the DST/SJ side spreads to 0.36–0.40; every CST drop to 0.49–0.67. ⚠ Only two DST drops, and no drop
lies within −0.022 … +0.017 of the boundary — the registered test `PREREG_DSTTYPE_2026-09-26` fills that gap and pairs
drops of matched thickening ratio across it.
