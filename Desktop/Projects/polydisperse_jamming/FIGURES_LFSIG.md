# Figure L1 — LF_DEM flow curves (PREREG_LFSIG)

`python lfdem/make_lfsig_figures.py` → `figures/L1_lfsig_flowcurves.{pdf,png,csv}`.

Lubricated DEM (LF_DEM, N = 500, Derjaguin-scaled repulsion λ = 0.02 a, Coulomb μ = 0.3), stress-controlled, lognormal
CV 0.2 (blue) and 0.5 (orange), at φ = 0.90 (open) and 0.95 (filled) of the dry-DEM φ_m for that CV. Points are block
means over strain 1.5–4, bars one block s.e. Stress in σ₀ = F_R/(6πa²), a = number-mean radius.
**(a)** relative viscosity; it thins slightly before it thickens, which is why the registered midpoint readout was
undefined. **(b)** contacts per particle; the four arms switch contacts on over the same decade (σ₀ 3–10), and the high-stress
contact number is lower at CV 0.5. Whether CV shifts the onset depends on how Z_c is normalised (CV ratio 0.96–1.17,
`PREREG_LFSIG` RESULT); no reading comes near the dry-DEM 0.768.
