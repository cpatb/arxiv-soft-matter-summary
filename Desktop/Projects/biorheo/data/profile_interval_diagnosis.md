# Item 13, step 1: the profile intervals sit on a bisection grid

Written by `scripts/profile_interval_diagnosis.py`. Reads the committed sidecars
and each model's own declared bounds; refits nothing.

## The signature

`dst_fit.profiling.profile_likelihood` sizes its outward step as
`max(wald/4, span/5000)` with `span = |box_bound - best|`, then bisects the
bracket 10 times and returns the OUTER end. If the statistic crosses on the
first step, the returned endpoint is `best + m * span / 5.12e6` for an integer
`m <= 1024`. Counted here as `|m - round(m)| < 1e-06`.

## What is on the grid

- sidecars read: **115**
- closed interval ends (an end at a box bound is not a crossing and is excluded): **609**
- of those, on the grid: **219** (36.0 %)
- closed UPPER ends: **296**, on the grid: **203**

## Worst case, by how far the endpoint sits above the value it describes

`blood_horner2020_donor_A_carreau_yasuda` · `carreau_yasuda` · `eta_0` (hi)

- best fit: `0.029909`
- endpoint: `49.0533`
- box bound: `1e+06`
- endpoint - value: `49.0234360338` = `251 x 1e+06/5.12e6` (grid index 251.0000000000)
- that is **1639x the best-fit value** above it

## Which (model, parameter) slots are affected

- `('casson', 'tau_y')`: 20
- `('casson', 'eta_c')`: 20
- `('carreau_yasuda', 'eta_inf')`: 20
- `('cross', 'eta_inf')`: 20
- `('newtonian', 'mu')`: 20
- `('power_law', 'k')`: 20
- `('carreau_yasuda', 'eta_0')`: 19
- `('cross', 'eta_0')`: 18
- `('carreau_yasuda', 'lam')`: 16
- `('herschel_bulkley', 'tau_y')`: 15
- `('herschel_bulkley', 'k')`: 15
- `('cross', 'lam')`: 13
- `('carreau_yasuda', 'a')`: 3

## What this does NOT show

It does not show the interval is too WIDE by a specific amount, and it is not a
measurement of where the F statistic actually crosses. The grid arithmetic says
only that the returned endpoint is a bisection endpoint. Whether each one is far
from the crossing is the direct check, and that needs the profile re-evaluated at
the returned endpoint -- gate G13.2.
