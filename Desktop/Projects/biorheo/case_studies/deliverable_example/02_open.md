# Injectability assessment — human whole blood, donor A (Horner 2020), hematocrit 42.6 %, fibrinogen 186 mg/dL

Device: **1 mL syringe, 0.210 mm bore x 12.7 mm needle, 1 mL/min manual push**. Fluid record `blood_horner2020_donor_A_carreau_yasuda`, model `carreau_yasuda`. Generated 2026-09-08 by `biorheo.apps.report`.

> **NOT LICENSABLE AS A DELIVERABLE.** The source data for this record is licensed *not stated (author-posted workbook; Horner 2020 dissertation citation requested)* and does not permit commercial use. This page is a demonstration of the method. A paid study must run on the client's own curve or on a curve whose licence permits it.

## The answer

**0.192 N** at 1 mL/min, with a one-at-a-time interval of [0.173, 29] N that is a **BOUND, not a limit** — the profile for a, lam, n ran into a box bound. **The device runs 1.42 decades outside the measured range**, so this is an extrapolation of the fitted law.

| quantity | value |
|---|---|
| plunger force | 0.1924 N |
| pressure drop across the needle | 0.1158 bar |
| wall shear stress | 47.86 Pa |
| wall shear rate | 1.85e+04 1/s |
| wall Reynolds number | 41.4 |

## Where your data stops and the device starts

The fitted envelope runs **0.1 to 700 1/s**. This device imposes a wall shear rate of **1.85e+04 1/s**.

That is **1.42 decades outside the measured range.** Every number above is the fitted law evaluated where nothing was measured. The law was chosen to fit the data you have; nothing in that data constrains its shape here. Measuring one point at the device's own rate is the cheapest thing that would change this page.

## What the answer depends on — and what it does not

`d ln F / d ln p` at the fit: a value near zero means the force does not care about that parameter at this operating point.

| parameter | value | 95 % interval | verdict | sensitivity |
|---|---|---|---|---|
| `eta_0` | 0.029909 Pa*s | [0.02799, 48.93] | measurable | +0.093 |
| `eta_inf` | 0.0023877 Pa*s | [0.002116, 0.003016] | measurable | +0.907 |
| `lam` | 7.9675 s | [6.404, 1e+06]* | degenerate | -0.035 |
| `a` | 10 dimensionless | [0.1, 10]* | degenerate | +0.000 |
| `n` | 0.58599 dimensionless | [0.05, 0.6129]* | bounded_only | +0.581 |

\* that interval ran into a box bound on at least one side. It is a **bound, not a confidence limit**.

**The force is insensitive to `lam`, `a`** (|d ln F / d ln p| < 0.05). Effort spent pinning those parameters down will not move this answer.

The dominant parameter is `eta_inf` at +0.907.

## A finer needle

`d ln F / d ln R` for this fluid at this operating point is **-3.894**, against **-4** for a Newtonian fluid. At a bore 30 % narrower the force goes to **0.7769 N**; a Newtonian `R^-4` sizing rule would have predicted **0.8014 N**, an over-prediction of **3.15 %**.

The Newtonian rule is safe for this fluid here, to within 10 %.

A **5 %** error in the bore moves the force by **19.5 %** (|b| x the declared tolerance). Compare that with the 1.5e+04 % parameter hull above: **the parameter uncertainty is the larger term**, so a better fit — or more data — buys more than a tighter bore.

## Against a force budget

Budget **25 N**. At the fitted parameters the budget is reached at **140.4 mL/min**; against the upper end of the parameter interval it is reached at **0.7755 mL/min**.

**Derating factor 0.005524** — a factor of 181 in deliverable flow rate, set by `eta_0` at its hi end. The **ratio** is the usable number; the absolute rate inherits every omission listed below.

This factor uses at least one interval that runs into a bound (`lam`, `a`, `n`), so it is itself a bound.

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
| record | `blood_horner2020_donor_A_carreau_yasuda` |
| constitutive model | `carreau_yasuda` |
| licence | not stated (author-posted workbook; Horner 2020 dissertation citation requested) |
| commercial use permitted | **NO** |
| source (paper) | https://doi.org/10.1007/s00397-023-01402-2 |
| needle | radius 0.105 mm, length 12.7 mm |
| syringe barrel radius | 2.3 mm |
| flow rate | 1 mL/min |
| flags | non_commercial, open_interval, extrapolated |

A gauge number is a range, not a dimension; the radius above is the one supplied with the request and no standard was consulted.
