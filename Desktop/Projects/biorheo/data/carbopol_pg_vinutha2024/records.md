# Carbopol 974p in propylene glycol: two licensable records, and the first real test of the temperature shift

Written by `scripts/records_from_carbopol_pg.py` from `data/carbopol_pg_prereg.md`. Source, licence and traps: `SOURCE.md`. Every number here is printed by that script.

## The records

| T | N | rates (1/s) | `tau_y` (Pa) | `k` | `n` | rmse_log |
|---|---|---|---|---|---|---|
| 10 C | 48 | 0.002-100 | 8.855 [8.682, 9.027] | 35.48 | 0.5389 [0.5339, 0.5438] | 0.01795 |
| 50 C | 60 | 0.00126-1000 | 9.176 [9.049, 9.303] | 12.28 | 0.5325 [0.5277, 0.5373] | 0.02166 |

## G24.1 — the fitted shift factor against the published viscosity ratio

`a_pred` = 0.012 / 0.102 = **0.11765** — the ratio of the medium viscosities the DEPOSIT publishes, which this repo does not compute.

`a_fit` = **0.13984**, from `collapse_shift` on the two flow curves alone. It never sees `a_pred`.

Ratio **1.1887**, `|ln ratio|` = 0.1728 against a bar of 0.2624. Verdict: **PASS**. Collapse residual 0.03385 on 98 of 120 rates.

## G24.2 and G24.3 — the test is shown informative, and the trap is measured

Below yield both curves tend to `tau_y`, so `eta -> tau_y/gd` and `a*eta_ref(a*gd) = tau_y/gd` for **any** `a`. The low-rate half of the curve therefore contributes nothing to the objective whatever the answer is. That was pre-registered, and here it is measured rather than argued.

Split at the geometric mean of the overlap, 0.4472 1/s.

| half | residual at `a_fit` | at `2 a_fit` | at `a_fit/2` |
|---|---|---|---|
| high rate | 0.02027 | 0.3356 | 0.2982 |
| low rate | 0.04945 | 0.08376 | 0.1576 |

**G24.2**: the high-rate residual rises by a factor of **16.6** when `a` is doubled, against a bar of 3.0. **PASS** — so G24.1 is informative.

**G24.3 FAILED**: the low-rate residual moves by 69.38 % over the same change, against a pre-registered bar of 10%. **The bar is not moved and the prediction is recorded as refuted.**

### G24.3ʺ — why G24.3 could not test what it was written to test

The section 1.1 algebra is right: below yield `eta -> tau_y/gd` and the shift is the identity there. But it only acts where the material is actually **on** the plateau, and G24.3 split at the geometric mean of the overlap, 0.4472 1/s, where this material is far ABOVE yield. Its 'low-rate half' was the power-law region.

Measured, not argued. Taking the plateau as `sigma < 1.2 * tau_y(10 C)` with `tau_y` = 8.855 Pa, it holds at **4 points, up to 0.00398 1/s**. But `collapse_shift` crops the bottom of the range by its own shift factor, so the shiftable window starts at **0.0143 1/s**.

**The plateau lies entirely OUTSIDE the shiftable window** (0.00398 < 0.0143). The trap pre-registered in section 1.1 is real algebra and **inoperative on this dataset**, because the plateau is narrower than the range the shift itself removes. That is a fact about the data, not about the shift layer — and it is why G24.2's separation is so large.

The probes were also re-run on the INTERSECTION of the three point sets, because comparing the residual at `a` with the residual at `2a` on different samples compares two things at once. On the common set the high-rate separation is 16.1 and the low-rate move is 129.82 %; both verdicts are unchanged.

## G24.4–G24.6 — the parameter map the shift predicts

Item 21 derived `tau_y -> tau_y`, `k -> a^n k`, `n -> n` before this dataset was found. The authors state the first of those from their own data: *"The yield stress is independent of the solvent viscosity."*

| gate | quantity | 10 C | 50 C | verdict |
|---|---|---|---|---|
| G24.4 | `tau_y` (Pa) | 8.855 [8.682, 9.027] | 9.176 [9.049, 9.303] | FAIL |
| G24.5 | `n` | 0.5389 [0.5339, 0.5438] | 0.5325 [0.5277, 0.5373] | PASS |

Separations, as the difference over the mean interval width: `tau_y` 1.07, `n` 0.653. An overlap of two WIDE intervals is a weak pass and is labelled one.

**G24.4 FAILED, and the bar is not moved.** The two `tau_y` profile intervals miss each other by **0.022 Pa**; the point estimates differ by 0.321 Pa, i.e. **3.56 %**. Item 21's derivation and the authors' own statement both say this parameter should be invariant, so a prediction is refuted here rather than a bar being adjusted.

**The scale that matters is beside it, and it is why this is a narrow failure rather than a broken model.** The gap is a fifth of a percent of the value, and a profile interval is a STATISTICAL interval on one curve: it contains no protocol-to-protocol term. G24.8 measures that term directly from the independent replicate and finds 2.93 % — larger than the difference being called significant. So the honest reading is that `tau_y` is invariant to within everything this dataset can resolve **except its own profile interval**, which is too narrow to be the right yardstick here.

**G24.6**: `k(50)/k(10)` = 0.34597 against the predicted `a_pred^n` = 0.31778 (with `n` = 0.5357), `|ln ratio|` = 0.0850 against 0.3365. **PASS**. This is a second route to G24.1 through the fitted parameters rather than the curve, so it could have disagreed.

## G24.7 — the 20 C block is refused, on signals recomputed here

| T | N | on the programmed log grid | rate reversals | stresses at <=3 sig figs |
|---|---|---|---|---|
| 10 C | 48 | 37 (77%) | **0** | 0 |
| 20 C | 67 | 3 (4%) | **2** | 67 |
| 50 C | 60 | 50 (83%) | **0** | 0 |

A rate that goes **backwards** inside a monotonically decreasing programmed ramp is not something the instrument does: 32.4 -> 33.6 1/s, 24.4 -> 26.4 1/s. Combined with the block being off the programmed grid and rounded to three significant figures, it cannot be certified as a raw curve at the grade of the other two, and this repo does not add a record it cannot certify. **No record is written from it.**

## G24.8 — an independent replicate, in the same units as the fit's uncertainty

The cessation files carry a second, unannounced steady flow curve (`preshear-stress` at 16 rates per temperature): the same material on the same instrument under a different protocol. It is **the only independent replicate this project has for any record**.

| T | matched rates | median ratio | range |
|---|---|---|---|
| 10 C | 5 | 1.0246 | 1.0078-1.0659 |
| 50 C | 5 | 1.0293 | 0.9748-1.1157 |

Worst replicate scatter 2.93 % against a worst relative parameter-interval width of 4.37 %. **The larger term is the parameter interval.** G15.4 had to approximate this comparison with a model misfit because no replicate existed.

## G24.9 — the published 8.9 Pa is NOT graded against these `tau_y`

The deposit publishes one yield stress, **8.9 Pa at all three temperatures**. It is the authors' fit to the **three-component model** (`sigma = sigma_y + k1 gd^0.5 + eta_bg gd`), which this package does not implement. biorheo fits Herschel-Bulkley.

A yield stress is a fit-objective convention before it is a material constant — Horner's Casson `tau_y` and item 23's flow index both say so — and comparing across two model forms would call the form difference an error. So the two numbers are reported side by side as a **model-form comparison and not a disagreement**: 10 C -> 8.855 Pa, 50 C -> 9.176 Pa, against the published 8.9 Pa.

**That the authors quote ONE yield stress for all three temperatures is used, and it is used as a prediction rather than as a number**: it is the same statement item 21 derived, and G24.4 tests it on this repo's own fits in this repo's own model.
