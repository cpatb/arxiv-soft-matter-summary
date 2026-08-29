# bio-fluids-bridge — blood rheology from the suspension-closure framework

A first step from the drop-impact / shear-thickening work into **biological fluids**,
showing that the same reduced-order **structural state-variable closure** used for
silica shear thickening also captures blood's shear thinning — one framework, opposite
sign of the structure→viscosity coupling.

Built 2026-07-16 (overnight, with Claude) as the entry rung of a possible biological-fluids track.

## The bridge in one line
- **Silica (drop-impact work):** state variable = fraction of *frictional contacts* `f`,
  RISES with stress → shear **thickening / jamming** (`tier0_wyart_cates_closure`).
- **Blood (here):** state variable = fraction of *aggregated red cells* (rouleaux) `λ`,
  FALLS with shear rate → shear **thinning**.

Both are first-order-kinetics structural closures. Only the sign of the coupling differs.
That transferability is the whole point — it's the concrete version of "my suspension /
jamming / closure toolkit ports to biological fluids."

## What it shows (validated, all gates PASS)
- **Fig 1** `fig1_blood_vs_silica.png` — apparent viscosity vs shear rate: blood thins from
  ~500 → ~4 mPa·s (rouleaux break up); silica thickens with a sharp DST jump. Same closure skeleton.
- **Fig 2** `fig2_thixotropic_loop.png` — blood thixotropic hysteresis loop under an up/down
  shear ramp (structure can't rebuild instantly on the down-ramp). This is the biological
  analogue of the silica *retention* hysteresis in the drop-impact reversibility demo.

## Model (see `blood_rheology.py`)
Structural / thixotropic rheology of the Moore–Coussot lineage (standard for blood):
`λ_ss(γ̇)=1/(1+γ̇/γ_c)`, `η(λ)=η_inf·(1+(η0/η_inf−1)λ)`, Casson-style yield offset,
transient `dλ/dt = k_build(1−λ) − k_break·γ̇·λ`. Physiological-ballpark params (Hct≈45%),
illustrative not fitted.

## Run
```
python run_blood_bridge.py   # PASS/FAIL validation report + 2 PNGs
```
Requires numpy, scipy, matplotlib.

## Honest caveats / next steps
- 0-D constitutive toy: homogeneous shear, single relaxation time, no Fåhræus–Lindqvist
  (confinement), no margination, no cell deformation. Casson is the textbook alternative law.
- **Calibrate** to real blood viscometry (η0, η_inf, γ_c, yield, kinetics).
- Natural extensions: margination of stiff cells (sickle/malaria) as a segregation problem;
  clotting as a *jamming* transition (reuse the silica jamming framework directly); confined
  flow (Fåhræus–Lindqvist) in a tube. See the vault's "Biological Fluids — Bridges" note.

## Portfolio note
Pairs with `tier0-closure` to tell one story: *a single state-variable rheology framework
spanning an industrial suspension (silica DST) and a biological fluid (blood) — calibrated,
validated, visual.* Strong, distinctive portfolio narrative for a complex-fluids-to-bio pivot.
