# Item 13 pre-registration — correcting the profile-likelihood interval

**Committed before any fix code was written.** Written to be graded, not to be passed.
Nothing below is a bar that may be moved; a bar that turns out to be the wrong bar is
recorded with its number and a *stricter* test is added beside it, as items 3, 6, 11 and
12 did.

## 0. What is already known, and what is new here

`dst_fit.profiling.profile_likelihood` walks outward from the best fit with

```
step = max(wald/4, span/5000),   span = |box_bound - best|
```

then, once the F statistic passes `f_crit`, runs `_REFINEMENTS = 10` bisections on the
bracket and returns its **outer** end. Two separate defects ride on that:

- **The step floor is sized from the parameter's declared box bound**, which is a
  guard-rail chosen so the fitter cannot wander, not a statement about the data. For a
  viscosity bounded at `1e6` the floor is `step = 200` on a parameter whose value is
  `0.0038`.
- **Ten bisections leave a resolution of `step/1024`.** When the crossing happens on the
  first step, the endpoint is `best + m * span/5.12e6` for an integer `m <= 1024` — a
  point on a bisection grid, not a likelihood crossing.

### 0.1 The reproduction, run before this file was written

`scripts/profile_interval_diagnosis.py` reads the committed sidecars and each model's own
declared bounds, and refits nothing. Result:

- 609 closed interval ends (an end sitting at a box bound is not a crossing and is excluded)
- **219** of them on the grid to one part in `1e-6` of an integer index
- **203 of 296 closed UPPER ends**
- worst: `blood_horner2020_donor_A_carreau_yasuda` `eta_0`, whose upper end sits
  **1639x the best-fit value** above it, at grid index `m = 251.0000000000`

**This corrects item 13's own numbers.** `NEXT_STEPS.md` says 117 of 296 and a worst of
82x, and flags that count as conservative because it tested only three candidate bounds.
It was conservative: reading each model's declared bound gives 203 of 296 and 1639x.
The item's diagnosis is confirmed and its magnitude was understated.

### 0.2 What the reproduction does NOT establish

The grid arithmetic shows only that the endpoint is a bisection endpoint. It is not a
measurement of where the F statistic crosses, and an endpoint could sit on the grid and
still be close to the crossing. **G13.2 is the gate that decides whether this matters at
all, and it is written so it can come out the other way** (§3).

## 1. Where the fix goes, decided in advance

**In `biorheo`, not in `dst_fit`.** The repo rule is that `dst_fit` is consumed unchanged,
and `tests/test_gate_0_1_dstfit_unchanged.py` exists to detect the *absence* of change —
its figure hash is pinned. A patch to someone else's library needs its own justification
and its own release; item 13 is not that justification. The wrapper lives at
`biorheo/identify/profiling.py` and `fit_rate_curve` calls it.

If the wrapper turns out to need a private name from `dst_fit`, that is recorded as a
coupling and pinned by a test that fails loudly when the private name moves.

## 2. What a corrected endpoint must satisfy

A profile confidence limit is, by definition, the parameter value at which the profile F
statistic equals `f_crit`. So the correctness criterion is *self-checkable and needs no
reference implementation*: evaluate the profile at the returned endpoint and ask whether
F is `f_crit`.

## 3. The gates

| id | claim | bar | what a failure means |
|---|---|---|---|
| **G13.1** | the stored endpoints sit on the bisection grid | at least one endpoint reproduces `m` to `1e-6`; count reported | already run, §0.1 |
| **G13.2** | the stored endpoint is *not* the crossing | `F(hi_stored) >= 2 * f_crit` on **>= 90 %** of the graded sample | **below 90 %, item 13 downgrades to cosmetic** and the sidecars are not regenerated |
| **G13.3** | the corrected endpoint *is* the crossing | `|F(hi_new) - f_crit| / f_crit <= 0.01` on **100 %** of ends the new routine closes | the fix does not work; do not ship it |
| **G13.4** | the corrected routine is right on a case with a closed-form answer | matches the analytic crossing to **<= 1e-6 relative** on all 20 blood `newtonian` records | the fix is wrong in a way G13.3 cannot see, because G13.3 uses the same machinery it is grading |
| **G13.5** | the correction can only *narrow* an interval | `hi_new <= hi_old` and `lo_new >= lo_old` on **>= 95 %** of ends both routines close | a violation means the two walks followed different branches of the profile; report it per case, do not average it away |
| **G13.6** | items 9 and 10 re-graded against **unchanged** bars | report `G9.1/9.2/9.3`, `G10.1/10.2/10.3` verdicts | **whatever comes out is recorded**; see §4 |
| **G13.7** | the loader's containment rule survives | every record value still inside its corrected sidecar interval | **a finding, not a bug to paper over**; see §5 |
| **G13.8** | `dst_fit` is untouched | `git status --porcelain` empty in the sibling checkout and the pinned figure md5 unchanged | the repo rule was broken |

### G13.2's sample, fixed here so it cannot be chosen after the fact

Every closed, on-grid **upper** end belonging to the 20 `blood_horner2020_donor_*`
`casson` records and the 15 Carbopol `herschel_bulkley` records. Both families are
two-parameter, so the sample is at most 70 ends; it is defined by a rule over the record
ids, not by inspection of which ones look bad.

### G13.4's closed form, derived here rather than reverse-engineered from the answer

The `newtonian` model is `eta(gd) = eta_0`, and `fit_rate_curve` minimises squared
residuals in `log eta`. Write `y_i = log eta_i` and `u = log eta_0`. Then

```
SSR(u) = sum_i (y_i - u)^2
```

which is minimised at `u_hat = ybar`, so

```
SSR(u) - SSR_min = n (u - ybar)^2
```

exactly — there are no nuisance parameters to re-minimise, so the profile is the full
likelihood. With `s^2 = SSR_min / (n - 1)` the statistic is
`F(u) = n (u - ybar)^2 / s^2`, and `F = f_crit` at

```
u = ybar +/- s * sqrt(f_crit / n),        eta_0 = exp(u)
```

with `f_crit = F_{1, n-1}(0.95)`. This touches neither the outward walk nor the root
finder, so it is an anchor and not a round trip.

## 4. What happens if items 9 or 10 flip — decided now, not after seeing the result

The item 9/10 seed-scatter metric divides a seed-to-seed span by the parameter's own 95 %
profile interval. An inflated denominator makes that metric too **small**, so those gates
passed more easily than they should have. The reported worst is `2.972e-05` against a
`1e-3` bar, a margin of **33.6x**, and the inflation on the worst-affected parameters is
of that order or larger.

- The bars `G91_BAR`, `G92_BAR`, `G93_BAR`, `G101_BAR`, `G102_BAR` are **not touched**.
- If a gate now FAILS, it is recorded FAILED on `multistart_scatter.md` /
  `carbopol_multistart_scatter.md` with the new number beside the old one, the way item
  9's own retraction is recorded, and the Status card and `NEXT_STEPS.md` say so.
- The old (inflated-denominator) value stays printed next to the new one. A number is not
  deleted because it turned out to be wrong; it is labelled.

## 5. What happens if containment fails — decided now

The loader requires a record's value to lie inside its sidecar interval. Narrowing every
interval can only make that harder, and the **Carbopol records are the pre-registered
risk**: their shipped values come from a different fitter (`carbopol/fit_hb.py`) and sit
0.1–4 % from the sidecar's log-eta fit, which is why the loader checks containment rather
than equality.

If a corrected interval excludes its record's shipped value, that is a **finding** — the
two fitters disagree by more than the corrected uncertainty — and it is reported as one.
The interval is **not** widened to accommodate the value, and the containment rule is
**not** relaxed.

## 6. Predictions, so that being wrong is visible

1. G13.2 passes at a high rate — most stored upper ends are far past `f_crit`.
2. Corrected intervals are **narrower**, in most cases dramatically (G13.5).
3. Items 9 and 10 still pass. *Reason:* the numerator of that metric is a seed-to-seed
   span of order `1e-5` of an interval, and the records that move are the ones already
   carrying `degenerate` / `not_identifiable`. But the margin is only 33.6x and the
   inflation is of that order, **so this prediction is the one most likely to be wrong**,
   and §4 governs if it is.
4. Containment holds on the blood records and is at risk on Carbopol.

Prediction 3 being refuted is a result. It is written down here so that it cannot be
quietly reread as having been expected.
