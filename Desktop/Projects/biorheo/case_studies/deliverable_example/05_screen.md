# Injectability assessment — model monoclonal antibody at 164.9 mg/mL with 25 mM L-arginine (formulation Arg_25)

Device: **1 mL syringe, 0.210 mm bore x 12.7 mm needle, 1 mL/min manual push**. Screen sample `mab_prasnikar2025_Arg_25` (NOT a registry record), model `newtonian`. Generated 2026-09-08 by `biorheo.apps.report`.

## The answer

**1.87 N** at 1 mL/min, [1.85, 1.88] N — an interval from the **measurement's own repeatability** (0.42 % SD), **not** a parameter interval. This is screen data: one viscosity, no flow curve, so there is no fit, no identifiability and nothing to say about model form. **The device runs 0.74 decades outside the measured range**, so this is an extrapolation of the Newtonian behaviour the source verified over that range.

| quantity | value |
|---|---|
| plunger force | 1.868 N |
| pressure drop across the needle | 1.124 bar |
| wall shear stress | 464.6 Pa |
| wall shear rate | 1.833e+04 1/s |
| wall Reynolds number | 4.22 |

## Where your data stops and the device starts

The range Newtonian behaviour was verified over runs **960 to 3329 1/s**. This device imposes a wall shear rate of **1.833e+04 1/s**.

That is **0.74 decades outside the measured range.** Every number above is a one-parameter Newtonian law evaluated where nothing was measured. Newtonian behaviour was verified over the range above and NOT here; nothing in a single viscosity constrains what the fluid does at the device's rate. Measuring one point at the device's own rate is the cheapest thing that would change this page.

## What the answer depends on — and what it does not

`d ln F / d ln p` at the fit: a value near zero means the force does not care about that parameter at this operating point.

| parameter | value | 95 % interval | verdict | sensitivity |
|---|---|---|---|---|
| `mu` | 0.025343 Pa*s | — | — | +1.000 |

**No parameter is insensitive here** — every one of them moves the force by more than 0.05 per decade of itself.

The dominant parameter is `mu` at +1.000.

## A finer needle

`d ln F / d ln R` for this fluid at this operating point is **-4**, against **-4** for a Newtonian fluid. At a bore 30 % narrower the force goes to **7.779 N**; a Newtonian `R^-4` sizing rule would have predicted **7.779 N**, an over-prediction of **6.22e-12 %**.

The Newtonian rule is safe for this fluid here, to within 10 %.

A **5 %** error in the bore moves the force by **20 %** (|b| x the declared tolerance). Compare that with the 1.68 % parameter hull above: **the bore tolerance is the larger term**, and tightening the fit will not help until the geometry is pinned down.

## Against a force budget

No derating factor: no interval: a derating factor needs one

## What this calculation leaves out

Laminar, fully developed, no-slip flow of a generalised-Newtonian fluid in a circular needle. Not included, and they do not all point the same way:

- the barrel's own pressure drop
- entrance and exit losses
- the kinetic term
- the transient at plunger start
- wall slip

Syringe friction is a property of the device, not the fluid; it is added only if the caller supplies it. Wall slip in particular is a boundary condition, not a constitutive extrapolation, and no flag on this page covers it — detecting it needs two needle radii at the same wall stress.

## Provenance

| field | value |
|---|---|
| record | `mab_prasnikar2025_Arg_25` |
| constitutive model | `newtonian` |
| licence | CC BY 4.0 |
| commercial use permitted | yes |
| source (dataset) | https://doi.org/10.5281/zenodo.14409284 |
| needle | radius 0.105 mm, length 12.7 mm |
| syringe barrel radius | 2.3 mm |
| flow rate | 1 mL/min |
| flags | extrapolated, screen_single_viscosity |

A gauge number is a range, not a dimension; the radius above is the one supplied with the request and no standard was consulted.
