# The pre-fix interval measurement, re-graded

Written by `scripts/profile_interval_regrade.py` from
`data/profile_interval_prefix_rows.json`, which is the row set committed at
`d38bc77`. **The working tree can no longer produce this**: item 13
regenerated every sidecar with the corrected profiler in that same commit, so
`profile_interval_diagnosis.py` run today reads the CORRECTED registry and finds
zero grid points. That is the right answer for today's registry and it is a
different question from the one this file answers.

## What the re-grade changed

The original test asked only whether `m = delta * 5.12e6 / span` is an integer.
With `_REFINEMENTS = 10` a bisection can only return `m` in
`[1, 1024]`, so an integer outside that range is not a grid
point. The range test was forced by NEXT_STEPS item 16, whose G16.1 failed on
exactly these rows.

- rejected as out of range: **18**
- of those, sitting at `m = 5.12e6`, i.e. `delta == span` and the endpoint AT the
  box bound rather than at any crossing: **5**

## The corrected counts

| quantity | as first reported | re-graded |
|---|---|---|
| grid points, all ends | 219 | **201** |
| closed UPPER ends on the grid | 203 of 296 | **201 of 296** |

The worst case is **unchanged**: `blood_horner2020_donor_A_carreau_yasuda` `eta_0`
(hi), whose endpoint sits **1639x the
best-fit value** above it, at a legitimate grid index `m = 251`.

## What this does and does not change

**It does not change item 13's conclusion or any of its gates.** G13.2 evaluated
the profile statistic AT the shipped endpoints and found a median of 710 times
`f_crit`; that measurement does not go through this detector at all. The
correction moves a COUNT by about one per cent and removes a class of false
positive that would otherwise have been quoted as evidence.

It does change the number that appears in prose: **201 of 296**, not 203.

## Which (model, parameter) slots are affected

- `('casson', 'tau_y')`: 20
- `('casson', 'eta_c')`: 20
- `('carreau_yasuda', 'eta_inf')`: 20
- `('cross', 'eta_inf')`: 20
- `('newtonian', 'mu')`: 20
- `('power_law', 'k')`: 20
- `('cross', 'eta_0')`: 16
- `('herschel_bulkley', 'tau_y')`: 15
- `('herschel_bulkley', 'k')`: 15
- `('carreau_yasuda', 'eta_0')`: 13
- `('carreau_yasuda', 'lam')`: 12
- `('cross', 'lam')`: 10
