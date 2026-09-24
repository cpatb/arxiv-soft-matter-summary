# FIGURES_CVMAP — the CV-at-fixed-weighed-φ drop (PREREG_CVMAP_2026-09-24)

`python make_cvmap_figures.py` regenerates everything (PNG byte-identical across runs, checked 2026-09-24).

## C1_cvmap_spread
Radial gyration radius of the drop, Rg_r/D0, against t/τ at x₀ = 0.932 (x of the CV-0 suspension) and
U₀ = 4 m/s. The weighed φ is the same for all curves. Black: CV 0, which is inside its window and jams
on impact. Blue: CV 0.2, still inside its window with J 0.04 (default cap) or 0.17 (EPS_MIN 0.03).
Orange: CV 0.5, which is outside its window and cannot jam. Dotted green: the window-only arm W,
which uses CV 0.5's x and χ with CV 0's σ\* and δ. Left panel: default EPS_MIN 0.01. Right panel:
EPS_MIN 0.03, where every cell of the grid survives. L9, cs_fit base, DEM inputs at μ_p 0.3.
Reading: the jammed drops arrest within t/τ ≈ 0.1. The CV-0.5 drop spreads to Rg ≈ 0.62, retracts
by t/τ ≈ 2, and ends 6.5 % (default cap: 19 %) above the jammed footprint, still creeping up at t/τ = 6.
Across the whole grid the terminal gap is sign-mixed (−8.2 % … +14.4 %); only the peak is one-signed. W overlies the full CV-0.5 curve, so the window
carries the effect. Data in `C1_cvmap_spread.csv`.
