# 58 high-concentration monoclonal antibody formulations through a needle

Written by `case_studies/mab_injectability/run.py`. Pre-registration: `data/mab_injectability_prereg.md`. Source and its own words: `data/mab_prasnikar2025/SOURCE.md`. Every number here is printed by that script.

**This is the first result in this project computed on the fluid the venture is actually about.** Everything else here is human blood, a Carbopol microgel, or cell culture medium.

## The headline

At matched mAb concentration the best and worst formulations differ in viscosity by **7.33x**. Through the Mooney slope that is worth **63 mg/mL of mAb concentration** — i.e. the right excipient buys you that much extra concentration at the same viscosity, against formulations that themselves span 53 mg/mL. Viscosity doubles every **21.8 mg/mL**.

| lever | size | in mg/mL of mAb |
|---|---|---|
| excipient, at matched concentration | 7.33x in viscosity | 63 |
| excipient, raw (CONFOUNDED — mixes both levers) | 3.47x | 39 |
| the formulations' own concentration span | — | 53 |

The raw row is reported because it is what a spreadsheet gives you and it is the wrong number: these formulations are not at one concentration, viscosity is exponential in concentration, so the raw spread mixes the excipient with the concentration. The normalised row is the answer.

## G25.3 — the force, on the blood study's own device

Needle inner diameter 0.21 mm, length 12.7 mm, barrel radius 2.3 mm, 1 mL/min — reused verbatim from `case_studies/blood_injectability`, so no new unverified constant enters and the two fluids sit on one axis for the first time.

| | min | median | max |
|---|---|---|---|
| mAb formulations (N = 58) | 1 N | 1.98 N | 3.47 N |
| blood donors, same device | 0.223 N | 0.28 N | 0.509 N |

The median mAb formulation needs **7.1x** the force of the median blood donor; 0 of 58 fall below it.

**No injectable-force limit is quoted, from a standard or anywhere else.** This project retracted an ISO 9626 comparison once for quoting a standard from memory, and a human-factors force limit would be the same mistake. The scale here is this repo's own blood study.

## G25.4 — how far the device extrapolates past what was measured

The source measured, and verified Newtonian behaviour over, **960–3329 1/s**. At this delivery rate the needle wall sees **1.833e+04–1.833e+04 1/s**, i.e. **0.74 decades** past the top of the measured range. Verdict: **PASS** against the pre-registered 1.5-decade bar.

| fluid | decades of extrapolation |
|---|---|
| **this mAb** | **0.74** |
| blood (item 5) | 1.42 |
| Carbopol in propylene glycol (item 24) | measured to 1e3 1/s |
| DMEM + FBS (item 22) | measured to 1 1/s, four to five decades |

*How far a device extrapolates its own rheology is a property of the DEVICE* — item 5's finding, now applied to the fluid the venture is for.

## G25.6 — the measurement is not the limiting term here, unlike blood

Median relative force uncertainty from the measured viscosity SD: **0.65 %**. **1525 of 1653** formulation pairs (92.3 %) have disjoint force intervals at 2 SD.

Item 14 asked the same question of the blood donors and got 36 of 190 (18.9 %) — there the uncertainty swamped the differences and two donors had to be 17 places apart to be told apart. **Here the formulations are separable**, so a formulation decision made on these numbers is supported by them in a way a patient-ranking decision on the blood data was not.

## G25.1 — a column of the source that does NOT reproduce, and is therefore not used

The paper defines its `viscosity reduction [%]` as the relative decrease against the starting formulation *at the same mAb concentration*, computed from an exponential (Mooney) calibration curve. Back-solving that baseline and fitting one exponential gives **R² = 0.4575** on the 49 rows of the main series (0.7855 over all 58). A single calibration curve would collapse this to numerical precision.

Either more than one curve was used, or the column carries processing the Methods do not describe. **The consequence is the point: this study uses the measured viscosities and never the reduction column.** A number that cannot be reproduced cannot be extended.

What does survive is the **slope**: 0.03131 mL/mg on the main series and 0.032334 on the higher-concentration one, agreeing to 3.2 %. The scatter that sinks G25.1 is in the intercept. That slope is what makes the headline's mg/mL conversion possible.

**That slope has no external check, and cannot have one from this source.** The paper fits the calibration curve `eta = a * exp(b * c)` over 110-200 mg/mL (Methods 2.2.4) by non-linear regression in R (2.2.7) and draws it with a 95 % confidence band (Fig. 1), but the fitted `a`, `b` and R-squared are printed nowhere in the accepted PDF or the deposited Table S1 (read 2026-09-28; `data/external_anchor_audit_2026-09-08.md` row 5). The authors' fitted exponent is therefore unpublished, and the doubling concentration quoted on this page is a number this study derived from the reduction column, not one the source reports.

## G25.7 — why these are NOT registry records, demonstrated rather than asserted

Fitting a one-parameter Newtonian model to the single measured viscosity of `Arg_25` **FAILED**:

```
ValueError: 1 points cannot determine 1 parameters
```

Zero residual degrees of freedom, so there is no profile interval, nothing for `write_identifiability` to write, and the loader would reject the record. This is a real boundary of the schema — and it is worth stating plainly because **it is also what a client's data usually looks like**: a table of formulation viscosities, not a flow curve. This study therefore reads the deposited table directly, which is the arriving-data path a paid study would actually take.

## Every formulation

| formulation | excipient | c(mAb) mg/mL | η mPa·s | ± SD | force N | η at matched c |
|---|---|---|---|---|---|---|
| S32_200_200 mAb | S32 | 212.2 | 20.34 | 0.421 | 1.5 | 5.44 |
| S36_200_200 mAb | S36 | 212.2 | 35.76 | 0.789 | 2.64 | 9.57 |
| Arg_200_200 mAb | Arg | 199.3 | 24.70 | 0.705 | 1.82 | 9.97 |
| Arg_100_200 mAb | Arg | 201.5 | 33.07 | 0.872 | 2.44 | 12.43 |
| S09_100mM | S09 | 180.3 | 19.52 | 0.139 | 1.44 | 14.41 |
| S36-S11 | S36-S11 | 179.8 | 19.57 | 0.093 | 1.44 | 14.68 |
| S36_100_200 mAb | S36 | 197.1 | 34.05 | 1.140 | 2.51 | 14.71 |
| Arg_100mM | Arg | 167.1 | 13.57 | 0.024 | 1 | 15.25 |
| Arg_200mM | Arg | 169.3 | 15.52 | 0.093 | 1.14 | 16.24 |
| S32_50_200 mAb | S32 | 191.2 | 31.41 | 1.489 | 2.31 | 16.39 |
| S36_50_200 mAb | S36 | 185.5 | 26.79 | 2.198 | 1.97 | 16.76 |
| Arg_50mM | Arg | 168.1 | 15.51 | 0.088 | 1.14 | 16.90 |
| S09_50mM | S09 | 170.5 | 17.35 | 0.219 | 1.28 | 17.48 |
| S32_200mM | S32 | 167.7 | 16.43 | 0.145 | 1.21 | 18.13 |
| S32_100_200 mAb | S32 | 190.8 | 34.38 | 1.224 | 2.53 | 18.15 |
| S36_200mM | S36 | 164.9 | 15.44 | 0.181 | 1.14 | 18.58 |
| S36_100mM | S36 | 164.7 | 15.92 | 0.178 | 1.17 | 19.29 |
| S36-S09 | S36-S09 | 166.9 | 17.41 | 0.157 | 1.28 | 19.71 |
| S36-S05 | S36-S05 | 166.6 | 17.31 | 0.164 | 1.28 | 19.79 |
| Arg_50_200 mAb | Arg | 190.7 | 38.04 | 0.479 | 2.8 | 20.19 |
| S32_50mM | S32 | 165.3 | 18.27 | 0.101 | 1.35 | 21.75 |
| S32-S09 | S32-S09 | 167.3 | 20.38 | 0.156 | 1.5 | 22.74 |
| S09_200mM | S09 | 176.5 | 27.65 | 0.046 | 2.04 | 23.05 |
| S36_50mM | S36 | 159.6 | 16.20 | 0.048 | 1.19 | 23.09 |
| S32-S17 | S32-S17 | 181.8 | 33.25 | 0.040 | 2.45 | 23.44 |
| S36-Pro | S36-S02 | 160.1 | 16.92 | 0.068 | 1.25 | 23.80 |
| S36-S17 | S36-S17 | 171.5 | 24.86 | 0.169 | 1.83 | 24.29 |
| S36_25 | S36 | 173.9 | 26.87 | 0.029 | 1.98 | 24.30 |
| S32_100mM | S32 | 162.8 | 18.97 | 0.175 | 1.4 | 24.43 |
| S32-S13 | S32-S13 | 176.9 | 30.19 | 0.156 | 2.23 | 24.85 |
| Pro_200mM | Pro | 170.9 | 25.27 | 0.234 | 1.86 | 25.19 |
| S32-S05 | S32-S05 | 168.1 | 23.22 | 0.091 | 1.71 | 25.29 |
| Pro_50mM | Pro | 175.0 | 29.48 | 0.608 | 2.17 | 25.80 |
| S32-S10 | S32-S10 | 167.9 | 23.69 | 0.671 | 1.75 | 25.96 |
| S32-S15 | S32-S15 | 175.2 | 30.15 | 0.107 | 2.22 | 26.23 |
| S32-Pro | S32-S02 | 177.0 | 32.21 | 0.089 | 2.37 | 26.43 |
| S32-S06 | S32-S06 | 162.0 | 20.14 | 0.077 | 1.48 | 26.65 |
| S17_25 | S17 | 173.4 | 29.02 | 0.114 | 2.14 | 26.68 |
| S36-S06 | S36-S06 | 177.2 | 32.93 | 0.470 | 2.43 | 26.84 |
| S36-S10 | S36-S10 | 172.4 | 29.16 | 0.127 | 2.15 | 27.67 |
| S09_25 | S09 | 177.1 | 34.17 | 0.013 | 2.52 | 27.97 |
| S32-S11 | S32-S11 | 169.7 | 27.05 | 0.646 | 1.99 | 28.00 |
| S36-S15 | S36-S15 | 168.9 | 26.62 | 0.332 | 1.96 | 28.21 |
| S11_25 | S11 | 181.7 | 40.91 | 0.029 | 3.01 | 28.89 |
| S36-S13 | S36-S13 | 179.6 | 38.34 | 0.244 | 2.83 | 28.94 |
| Pro_100mM | Pro | 164.5 | 23.97 | 0.226 | 1.77 | 29.25 |
| S36-S31 | S36-S31 | 171.9 | 31.02 | 0.149 | 2.29 | 29.97 |
| S06_25 | S06 | 160.7 | 22.04 | 0.059 | 1.62 | 30.39 |
| Arg_25 | Arg | 164.9 | 25.34 | 0.106 | 1.87 | 30.55 |
| S32-S36 | S32-S36 | 184.2 | 47.07 | 1.105 | 3.47 | 30.70 |
| S13_25 | S13 | 170.7 | 32.61 | 0.042 | 2.4 | 32.71 |
| S05_25 | S05 | 170.1 | 34.04 | 0.062 | 2.51 | 34.79 |
| Pro_25 | Pro | 159.3 | 24.16 | 0.234 | 1.78 | 34.80 |
| S10_25 | S10 | 166.5 | 30.85 | 0.049 | 2.27 | 35.36 |
| S31_25 | S31 | 162.9 | 27.73 | 0.084 | 2.04 | 35.64 |
| S32_25 | S32 | 163.4 | 28.28 | 0.055 | 2.08 | 35.75 |
| S15_25 | S15 | 166.5 | 31.52 | 0.007 | 2.32 | 36.08 |
| S32-S31 | S32-S31 | 173.9 | 44.08 | 0.290 | 3.25 | 39.91 |

## What this page does NOT say

- **Nothing about whether any of these is injectable by a human.** No force limit is quoted. The scale is this repo's own blood study on the same device.
- **Nothing about stability.** The source carries SEC monomer and aggregate columns; stability is not rheology and this package has nothing to say about it. An excipient that fixes viscosity and ruins the protein is not a win, and that trade-off is in the source and out of scope here.
- **Nothing that uses the source's `viscosity reduction` column.** G25.1.
- **Nothing about shear-thinning at the needle wall.** The source verified Newtonian over 960–3329 1/s and that is the model used; whether it still holds 0.74 decades higher is what that number is for, and it is not assumed either way.
