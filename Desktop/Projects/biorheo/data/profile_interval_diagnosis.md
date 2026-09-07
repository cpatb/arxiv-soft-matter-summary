# Is the registry still on the bisection grid? (regression check)

Written by `scripts/profile_interval_diagnosis.py`. Reads the committed sidecars
and each model's own declared bounds; refits nothing.

**Item 13 regenerated every sidecar with the corrected profiler, so this should
report ZERO.** The historical measurement -- how bad it was before that -- is
`data/profile_interval_prefix.md`, re-graded from the row set committed at
`d38bc77`. Do not quote this file for the size of the defect; it is the check that
the defect is gone.

## The signature

`dst_fit.profiling.profile_likelihood` sizes its outward step as
`max(wald/4, span/5000)` with `span = |box_bound - best|`, then bisects the
bracket 10 times and returns the OUTER end. If the statistic crosses on the
first step, the returned endpoint is `best + m * span / 5.12e6` for an integer
`m <= 1024`. Counted here as `|m - round(m)| < 1e-06` **and**
`1 <= m <= 1024`.

**The range half of that test was added on 2026-09-06** (NEXT_STEPS item 16).
Without it the test accepted any `delta/span` that happened to be a round
fraction -- including every end whose `delta` EQUALS its `span`, which is
`m = 5.12e+06` and means the endpoint sits AT the box bound rather than
at a crossing. **0** ends are integer-indexed but out of range and
are now excluded; 5 of them sit at exactly 5120000. The bound is derived from
`_REFINEMENTS = 10`, not chosen to make anything pass.

## What is on the grid

- sidecars read: **115**
- closed interval ends (an end at a box bound is not a crossing and is excluded): **601**
- of those, on the grid: **0** (0.0 %)
- closed UPPER ends: **292**, on the grid: **0**

## Which (model, parameter) slots are affected


## What this does NOT show

It does not show the interval is too WIDE by a specific amount, and it is not a
measurement of where the F statistic actually crosses. The grid arithmetic says
only that the returned endpoint is a bisection endpoint. Whether each one is far
from the crossing is the direct check, and that needs the profile re-evaluated at
the returned endpoint -- gate G13.2.
