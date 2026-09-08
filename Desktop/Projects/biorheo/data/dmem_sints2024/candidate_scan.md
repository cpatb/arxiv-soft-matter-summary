# DMEM + FBS (Sints et al. 2024) — candidate scan, NOT a registry record

Written by `scripts/dmem_candidate_scan.py`. **No record is created by this file and no verdict is shipped.** It exists so that a later session decides on evidence rather than on a first impression. Read the script's docstring before the table.

## Why it was scanned

This is the first openly-licensed (CC BY 4.0) rheometer dataset found for this project. Its source paper states that DMEM + FBS *"requires a more complex treatment (i.e., is a non-Newtonian fluid)"*. In every interval the **shear stress column is close to affine in shear rate with a non-zero intercept**, and a viscosity formed as `tau/gd` from such a series falls like `tau_0/gd + eta_p` — which looks like shear thinning on a log-log plot whatever `tau_0` is. So the narrow, answerable question is whether a two-parameter **Bingham** form describes these curves as well as a two-parameter shear-thinning **power law**.

**A Bingham fluid is also non-Newtonian, so this does not contradict the paper's headline.** What it would change is the mechanism, and everything downstream: a stress intercept and a thinning exponent extrapolate to a channel completely differently.

## Result: **the hypothesis is REFUTED.** Bingham beats the power law on only 1 of 8 intervals

The scan was run to test whether an affine stress column was masquerading as shear thinning. It is not. **The two-parameter power law fits better on 7 of 8 intervals**, at the same parameter count and on the same error measure, with `n` between 0.30 and 0.70. **The paper's shear-thinning reading is the one this data supports.** The hypothesis that motivated the scan is kept in the script's docstring rather than deleted, because a scan is only worth reading if what it was looking for is visible.

All three descriptions are compared on the same quantity — rms error in log viscosity — so the ranking is not a units artefact. `newtonian` is the one-parameter null; the other two have two parameters each. Both two-parameter forms beat the Newtonian null on 8 of 8, which is the least this data could have shown and does not by itself establish a form.

| interval | pts | rate range (1/s) | eta range (Pa s) | Newtonian rmse | Bingham rmse | power-law rmse | `tau_0` (Pa) | `eta_p` (Pa s) | `n` |
|---|---|---|---|---|---|---|---|---|---|
| `DMEM_and_10pct_FBS #1` | 14 | 0.0714–1 | 0.0286–0.185 | 0.552 | **0.060** | 0.035 | 0.01415 | 0.01531 | 0.301 |
| `DMEM_and_10pct_FBS #2` | 15 | 0.01–1 | 0.0266–0.609 | 1.067 | **0.209** | 0.040 | 0.009827 | 0.01996 | 0.302 |
| `DMEM_and_10pct_FBS_37C #1` | 5 | 0.01–1 | 0.00552–0.0657 | 1.005 | **0.250** | 0.151 | 0.001041 | 0.004537 | 0.492 |
| `DMEM_and_10pct_FBS_37C #2` | 20 | 0.01–10 | 0.00208–0.0493 | 1.081 | **0.458** | 0.492 | 0.001299 | 0.002009 | 0.700 |
| `DMEM_and_1pct_FBS #1` | 15 | 0.01–1 | 0.0147–0.109 | 0.690 | **0.340** | 0.096 | 0.002676 | 0.01393 | 0.507 |
| `DMEM_and_1pct_FBS #2` | 15 | 0.01–1 | 0.0159–0.0938 | 0.602 | **0.456** | 0.202 | 0.002819 | 0.01506 | 0.507 |
| `DMEM_and_5pct_FBS #1` | 15 | 0.01–1 | 0.0388–0.792 | 0.978 | **0.268** | 0.108 | 0.01281 | 0.03201 | 0.345 |
| `DMEM_and_5pct_FBS #2` | 15 | 0.01–1 | 0.0338–0.873 | 1.110 | **0.147** | 0.044 | 0.01194 | 0.02564 | 0.310 |

`1` point(s) were excluded across all intervals by a stated rule — shear rate below 1e-3 1/s — and the count is printed rather than the exclusion being silent. One file records a first point at 2.24e-7 1/s with an apparent viscosity of 8710 Pa s.

## What the numbers say, and what they do not

- The fitted intercept `tau_0` runs **0.001041 to 0.01415 Pa** — millipascals, which is small in absolute terms and is most of the stress at the bottom of a rate range that starts at 0.01 1/s.
- The plastic viscosity left after removing it runs **0.002009 to 0.03201 Pa s**. Water is about 1e-3 Pa s, so the coldest, most serum-rich interval still sits about 32x water while the warmest sits about 2x.
- `tau_0` is **not constant across the four samples**: it tracks FBS content and falls with temperature. That would have argued against a fixed instrument offset — but since the Bingham form is the WORSE description here, `tau_0` is a fitted parameter of a rejected model and should not be quoted as a yield stress at all.

**The magnitude is the PAPER'S OWN HEADLINE, not an oddity this scan found.** The room-temperature intervals sit 14–32x water and the 37 °C ones 2–5x, for a fluid that is ~99 % water. Reading the source (2026-09-08) settles it: the paper reports that *"while the measured values do converge to the literature value as the shear rate increases, the viscosity at lower shear rates is much grater [sic] than the literature value"*, against a baseline of 9.3e-4 Pa s. That gap is the result the paper exists to report.

**The scan's fitted indices reproduce the paper's own, without having seen them.** The authors fit the Ostwald-de Waele power law and quote n = 0.32 (10 % FBS at 20 °C), 0.31 (10 % and 5 % at 20 °C) and 0.5–0.54 (10 % at 37 °C, and 1 % at 20 °C). This scan gets 0.301/0.302 for 10 % at 20 °C, 0.345/0.310 for 5 % at 20 °C, and 0.492 and 0.507 for the two the paper groups at 0.5–0.54. That agreement is the best available evidence that the parse, the interval splitting and the zero-rate exclusion are all correct. **The one that agrees least is named rather than omitted**: the 37 °C 20-point interval gives 0.700, and it is the only interval running to 10 1/s rather than 1 1/s, so it is fitted over a decade the others do not cover and a single power law is not expected to hold across both.

**What a later session should do:** the form is settled (`power_law`, the family the authors themselves fit). Fit it through `biorheo.identify` for intervals and an identifiability verdict, pre-register a gate against the paper's published indices, and write the record. See `SOURCE.md` for the ordered list and for the pre-shear protocol any record must inherit.

## The other reason this dataset cannot carry an injectability claim

The measured rates stop at 1 1/s (10 1/s on one interval). A fine needle imposes 1e4–1e5 1/s. Any force computed from a record built here would be **four to five decades outside the measured range**, against 1.42 decades for the blood records. That makes it an excellent demonstration of the extrapolation flag and a poor basis for a number anyone should act on.
