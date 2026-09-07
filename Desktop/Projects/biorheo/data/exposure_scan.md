# Item 16: which parameters were exposed to the interval defect, and by how much

Written by `scripts/exposure_scan.py`. Gate pre-registered in
`data/exposure_prereg.md` (commit `1e33b56`) before this was run. **That file
opens by retracting a sentence item 15 shipped; read it first.**

## The quantity

The defect is a step floor of `span/5000` from the parameter's declared BOX,
applied as `step = max(wald/4, span/5000)`. It acted when the floor won, and the
amount by which it overreached is

```
inflation = (span / 5000) / (wald / 4)
```

**a ratio, not a boolean.** `inflation = 1.09` means the floor was applied and
happened to be the right size; `inflation = 1.2e+04` means the first step
overshot by four orders of magnitude. Saying only *whether* the box was consulted
is what made item 15's explanation wrong.

| gate | claim | bar | result | verdict |
|---|---|---|---|---|
| G16.1 | the grid signature implies the box floor won | 0 counterexamples | **0** over 609 closed interval ENDS (201 on the grid) | **PASS** |
| G16.2 | the converse is NOT claimed | reported | 8 inflated but not on the grid, of 209 inflated | reported |
| G16.3 | the ratio tracks the overreach | reported | Spearman 0.784 (p = 5.2e-43, n = 201) | reported |

## What G16.1 establishes, and what it does not

Every parameter whose stored endpoint sits on the bisection grid has
`inflation > 1`, so **the mechanism written on four pages of this repo is the
right one**. The converse is deliberately not claimed: a crossing found on a
LATER step lands on `step0 * 1.18^j / 1024`, a different grid that
`profile_interval_diagnosis.py` does not test for, so a parameter can be inflated
without being detectably on the grid. That cell holds **8**
ends and is a limit of the DETECTOR, not evidence against the mechanism.

## How exposed the registry is

- closed interval ENDS scanned: **609**. Each parameter contributes up
  to two, and the two differ: a viscosity at `best = 21` has `span_hi/5000 = 200`
  against `span_lo/5000 = 0.0043`, so an end must be graded against its own bound.
- with `inflation > 1` (the floor won): **209**
- with a detectable grid endpoint: **201**
- `inflation` quantiles: min 2.19e-09, q25 0.00821, median 0.0334, q75 2.22e+04, max 1.16e+07

## Worst and best exposed (model, parameter)

| record | model | parameter | end | wald/4 | span/5000 | inflation | on grid |
|---|---|---|---|---|---|---|---|
| `blood_horner2020_donor_U_carreau_yasuda` | `carreau_yasuda` | `eta_inf` | hi | 1.728e-05 | 200 | **1.157e+07** | yes |
| `blood_horner2020_donor_U_cross` | `cross` | `eta_inf` | hi | 1.867e-05 | 200 | **1.071e+07** | yes |
| `blood_horner2020_donor_C_carreau_yasuda` | `carreau_yasuda` | `eta_inf` | hi | 2.199e-05 | 200 | **9.095e+06** | yes |
| `blood_horner2020_donor_I_carreau_yasuda` | `carreau_yasuda` | `eta_inf` | hi | 2.258e-05 | 200 | **8.857e+06** | yes |
| `blood_horner2020_donor_S_cross` | `cross` | `eta_inf` | hi | 2.29e-05 | 200 | **8.735e+06** | yes |
| `blood_horner2020_donor_B_cross` | `cross` | `lam` | lo | 6.338e+07 | 200 | **3.156e-06** | no |
| `blood_horner2020_donor_L_carreau_yasuda` | `carreau_yasuda` | `eta_0` | lo | 1.805e+06 | 0.006663 | **3.692e-09** | no |
| `blood_horner2020_donor_B_carreau_yasuda` | `carreau_yasuda` | `eta_0` | lo | 1.204e+06 | 0.004273 | **3.548e-09** | no |
| `blood_horner2020_donor_L_carreau_yasuda` | `carreau_yasuda` | `lam` | lo | 9.998e+09 | 22.54 | **2.255e-09** | no |
| `blood_horner2020_donor_B_carreau_yasuda` | `carreau_yasuda` | `lam` | lo | 1.263e+10 | 27.63 | **2.187e-09** | no |

## What this is NOT

**Not an impact estimator.** `inflation` says how far the first step overshot for
one parameter. It does not say how wrong a downstream number was, because that
also depends on which parameter the number is sensitive to (item 15: a needle
force set by a flow index barely moved while its `k` interval moved 2.94 %) and on
which END of the interval the comparison binds (item 13's G13.9: the defect
inflated the upper ends and the comparison used the lower ones). Neither is
derivable from this ratio and neither should be inferred from it.
