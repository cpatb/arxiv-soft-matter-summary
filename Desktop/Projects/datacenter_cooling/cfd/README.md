# Tier 1 — CFD calibration of the recirculation matrix

This closes the one soft spot in the Tier-0 model: the recirculation influence `B` was a
geometric placeholder. Here a **validated 2D buoyant-convection CFD** *measures* the
recirculation fraction φ vs containment effectiveness η, and that curve replaces the
placeholder in `../dc_thermal.py`.

## What's here
```
buoyant2d.py        streamfunction-vorticity Boussinesq solver (fast DST Poisson)
test_buoyant2d.py   VALIDATION vs de Vahl Davis (1983): Ra=1e3 -> 0.3% err, 1e4 -> 1.8% err
aisle.py            vented-partition hot/cold-aisle cell; sparse-LU Poisson with the
                    rack+containment baffle as an internal solid; phi = inlet warming
run_calibration.py  sweeps phi(eta) over containment and two heat loads -> cfd_recirc.json
calibrate.py        fits phi(eta), re-optimizes the hall, compares to the placeholder
fig_cfd_calibration.png   the calibration + its effect on optimized PUE
```

## Run
```
pip install numpy scipy matplotlib --break-system-packages
python test_buoyant2d.py     # solver validation (de Vahl Davis)
python run_calibration.py    # ~2 min: writes cfd_recirc.json (timeout-safe, resumable)
python calibrate.py          # fit + re-optimization + figure
```

## The result (calibrated 2026-07-23)
The de Vahl Davis-validated solver, applied to a hot-aisle/cold-aisle cross-section, gives

    phi(eta) = 0.007 + 0.228 * exp(-eta / 0.790)      # recirculation fraction vs containment

with only weak dependence on heat load (Ra 4e3 vs 8e3). Two findings the placeholder missed:

1. **A recirculation floor.** Even at full containment (η→1) φ ≈ 0.07, not zero. The Tier-0
   placeholder `0.28(1-η)` drove recirculation to zero at η=1, over-crediting containment.
2. **Concave, diminishing returns.** Most of the benefit is in the first ~0.4 of containment;
   past η≈0.5 the optimized-PUE curve flattens (right panel of the figure). Practically: the
   ROI on chasing the *last* increment of containment is much smaller than the reduced-order
   model assumed — a decision-relevant, non-obvious result and exactly the kind of thing a CFD
   calibration is for.

φ(0) ≈ 0.235 from CFD, close to the placeholder's 0.28 — so the near-field guess was fine; the
*shape* was wrong. `dc_thermal.py` now uses this curve by default (`recirc_model="cfd"`); pass
`recirc_model="linear"` for the old placeholder.

## Honest scope / limitations
- 2D cross-section: captures the dominant self-recirculation (the diagonal of `B`) and its
  containment dependence. Rack-to-rack (off-diagonal) coupling and along-row effects need a 3D
  or plan-view run — that is what the OpenFOAM case in `../openfoam/` is for.
- The partition's thermal boundary is treated as adiabatic with a first-order immersed method;
  φ is the normalized inlet warming, which converges with the flow field. Absolute φ carries the
  usual coarse-grid / 2D caveats (validated trend, ~10-20% on absolute values).
- Heat-load (Ra) dependence is weak in the tested range; a wider sweep belongs in the 3D study.
