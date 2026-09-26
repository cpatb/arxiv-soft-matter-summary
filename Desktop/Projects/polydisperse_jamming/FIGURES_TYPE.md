# Figure T1 — the thickening-type map (`make_type_figures.py`, `figures/T1_thickening_type.{pdf,png,csv}`)

Regenerate: `python make_type_figures.py --check` (PNG byte-identical on regeneration).

**(a) Where a Wyart–Cates suspension is CST, DST or shear-jamming**, in x = φ/φ₀ against χ = φ_m/φ₀ (closed form (T),
`derivations/thickening_type_2026-09-25.md`). Solid black: x = χ, the lower edge of shear jamming (SJ above it).
Broken curves: the DST onset x_DST(χ, nδ) for nδ = 2.147 (blue dashed; the CV 0.5 activation with n = 2), 2.634
(orange dash-dot; CV 0 activation, n = 2) and 6.58 (green dotted; CV 0 activation with the cornstarch fit's n = 5).
Between a broken curve and the black line the steady flow curve is S-shaped (DST); below it thickening is continuous
(CST). The DST band shrinks to nothing as χ → 1 and widens with nδ; it exists only for nδ > 1.

**(b) The registered DNS test (`PREREG_DSTTYPE_2026-09-26`, 20/20 healthy):** peak spread Rg_max/D₀ against the
thickening ratio Θ = [(1−x)/(1−x/χ)]^n (n 2) for two input sets (δ 1.07/σ* 220 and δ 1.32/σ* 300) at U₀ 4 (large
markers) and 2 (small). Blue circles: CST (x < x_DST); orange diamonds: DST (x_DST ≤ x < χ). All points fall on one
monotone curve per U₀ with no break at the type boundary: the drop follows Θ, not the steady-shear type (D1 THETA,
D2 FAIL). The earlier post-hoc reading of the 24 cvmap drops ("DST drops spread like jammed ones",
`derivations/cvmap_by_type_2026-09-26.txt`) was Θ in disguise — DST-band drops are the ones with large Θ.
