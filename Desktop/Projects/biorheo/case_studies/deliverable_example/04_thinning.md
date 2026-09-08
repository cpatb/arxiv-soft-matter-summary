# Injectability assessment — Carbopol ETD2020 0.5 wt%, 40 uL triethanolamine per 100 g, in water

Device: **1 mL syringe, 0.210 mm bore x 12.7 mm needle, 1 mL/min manual push**. Fluid record `carbopol_etd2020_0p5wt_40uL`, model `herschel_bulkley`. Generated 2026-09-08 by `biorheo.apps.report`.

> **NOT LICENSABLE AS A DELIVERABLE.** The source data for this record is licensed *CC BY-NC-SA 4.0* and does not permit commercial use. This page is a demonstration of the method. A paid study must run on the client's own curve or on a curve whose licence permits it.

## The answer

**3.55 N** at 1 mL/min, [3.5, 3.59] N over the 95 % one-at-a-time parameter hull (2.68 % wide). **The device runs 2.08 decades outside the measured range**, so this is an extrapolation of the fitted law.

| quantity | value |
|---|---|
| plunger force | 3.548 N |
| pressure drop across the needle | 2.135 bar |
| wall shear stress | 882.4 Pa |
| wall shear rate | 2.422e+04 1/s |
| wall Reynolds number | 2.94 |
| unyielded plug, fraction of the bore | 0.0269 |

## Where your data stops and the device starts

The fitted envelope runs **0.000999 to 200 1/s**. This device imposes a wall shear rate of **2.422e+04 1/s**.

That is **2.08 decades outside the measured range.** Every number above is the fitted law evaluated where nothing was measured. The law was chosen to fit the data you have; nothing in that data constrains its shape here. Measuring one point at the device's own rate is the cheapest thing that would change this page.

## What the answer depends on — and what it does not

`d ln F / d ln p` at the fit: a value near zero means the force does not care about that parameter at this operating point.

| parameter | value | 95 % interval | verdict | sensitivity |
|---|---|---|---|---|
| `tau_y` | 23.756 Pa | [23.72, 23.79] | measurable | +0.033 |
| `k` | 8.903 Pa*s**0.4526 | [8.85, 8.957] | measurable | +0.967 |
| `n` | 0.4526 dimensionless | [0.4511, 0.454] | measurable | +4.229 |

**The force is insensitive to `tau_y`** (|d ln F / d ln p| < 0.05). Effort spent pinning those parameters down will not move this answer.

The dominant parameter is `n` at +4.229.

## A finer needle

`d ln F / d ln R` for this fluid at this operating point is **-2.313**, against **-4** for a Newtonian fluid. At a bore 30 % narrower the force goes to **8.12 N**; a Newtonian `R^-4` sizing rule would have predicted **14.78 N**, an over-prediction of **82 %**.

**The Newtonian rule is not safe for this fluid here.** Sizing a needle with `R^-4` will mis-state the force by more than 10 %.

A **5 %** error in the bore moves the force by **11.6 %** (|b| x the declared tolerance). Compare that with the 2.68 % parameter hull above: **the bore tolerance is the larger term**, and tightening the fit will not help until the geometry is pinned down.

## Against a force budget

Budget **25 N**. At the fitted parameters the budget is reached at **79.72 mL/min**; against the upper end of the parameter interval it is reached at **76.35 mL/min**.

**Derating factor 0.9578** — a factor of 1.04 in deliverable flow rate, set by `n` at its hi end. The **ratio** is the usable number; the absolute rate inherits every omission listed below.

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
| record | `carbopol_etd2020_0p5wt_40uL` |
| constitutive model | `herschel_bulkley` |
| licence | CC BY-NC-SA 4.0 |
| commercial use permitted | **NO** |
| source (dataset) | https://doi.org/10.25835/hjn2qcvy |
| source (paper) | https://doi.org/10.3390/ma18133164 |
| needle | radius 0.105 mm, length 12.7 mm |
| syringe barrel radius | 2.3 mm |
| flow rate | 1 mL/min |
| flags | non_commercial, extrapolated |

A gauge number is a range, not a dimension; the radius above is the one supplied with the request and no standard was consulted.
