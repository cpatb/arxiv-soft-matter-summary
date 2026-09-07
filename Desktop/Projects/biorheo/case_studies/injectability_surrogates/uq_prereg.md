# Item 15 pre-registration — the Carbopol study's intervals predate the item-13 fix

**Committed before anything was recomputed.** Bars are not moved; a bar that turns out to
be the wrong bar is recorded with its number and a stricter test is added beside it. A
pre-registration is never rewritten; corrections go in an addendum at the foot.

## 0. Item 15's stated premise is FALSE, and checking it took one grep

`NEXT_STEPS.md` item 15 — written earlier the same day — says *"the Carbopol case study
ships point forces too"* and asks for a UQ column to be added.

**It does not ship point forces.** `run.py` has called `force_with_uq` since it was
written, and `results.md` has carried a *"one-at-a-time profile envelope 95 %"* column, a
*"covariance 95 %"* column, a per-parameter sensitivity line and an `Insensitive to:` line
all along. The figure already draws both bands.

This is the fourth time in this repo that an item's stated premise has been wrong (item 6's
`BirdCarreau` `a`, item 11's `structural_steady` counterpart, item 13's own understated
diagnosis, and now this), and it is recorded rather than worked around.

**What is actually true is worse and more useful.** Those intervals were computed by
`dst_fit.profiling` before NEXT_STEPS item 13, so **the study's shipped error bars are
bisection-grid artefacts**, exactly as the registry's were. This item is a CORRECTION to a
shipped page, not an extension of one.

### 0.1 What is not in doubt

- `run.py` REFITS each gel rather than reading registry parameters, so the correction
  reaches it through `fit_rate_curve` with no change to this case study's own code.
- All 15 Carbopol records were graded `measurable` on every parameter by
  `identify_carbopol.py`, so an OPEN profile end is unlikely here. The "an open interval
  is never printed as a number" branch in the sibling study has still never fired.

## 1. Two gates that CANNOT fail, named here so they are not later quoted as evidence

Item 14's lesson was that the obvious gate was unfalsifiable. Two here are:

- **The ordering reversal survives.** At the 10 s dose the two gels' shipped hulls are
  `[19.7, 20.8]` N and `[54.8, 113.8]` N — disjoint by a factor of nearly three, against a
  3.9x effect. Narrowing the intervals can only make them more disjoint. This is reported
  as a number and **is not counted as a gate**.
- **The nominal forces are unchanged.** Item 13 touched only the interval walk, not the
  fit. Reported, not counted.

## 2. The gates

| id | claim | bar | what a failure means |
|---|---|---|---|
| **G15.1** | the correction reaches this page at all | the nominal forces agree to **<= 1e-9** relative AND at least one interval endpoint moves by **> 1 %** | if no interval moves, this study was never affected and the item closes as a no-op |
| **G15.2** | the corrected intervals NARROW | new hull inside old hull on **>= 5 of 6** rows | a widening means the two walks followed different branches; report per case |
| **G15.3** | no end is silently open | any parameter with `lo_at_bound` or `hi_at_bound` is reported as OPEN, not as a number; count reported | a number would be quoted for a limit the data does not set |
| **G15.4** | **the parameter interval is still the dominant error term** | for each gel, the hull's relative half-width vs `100 x rmse_log` | see §3 — **this is the gate that can flip** |

### G15.4, stated precisely

`rmse_log` is the fit's own root-mean-square residual in log viscosity, i.e. a relative
viscosity error, and `d ln F / d ln k` is 0.979 and 1.000 for the two gels — the force is
very nearly proportional to the consistency. So a viscosity misfit of `rmse_log`
propagates to approximately `rmse_log` in force, and the two quantities are directly
comparable as percentages.

Shipped (pre-correction) values:

| gel | hull half-width at 10 s | `rmse_log` | which dominates |
|---|---|---|---|
| 0.5 wt% 50 uL | ±2.6 % | 0.4 % | parameter interval |
| 0.1 wt% 40 uL | ±37 % | 16.6 % | parameter interval |

**The question is whether that survives.** The correction can only narrow the intervals,
so a large enough narrowing flips one or both rows to *model-form dominates* — which would
mean the page's error bars stop being the honest total error and it must say so.

**Bar:** report the ratio `hull_half_width_pct / (100 * rmse_log)` per gel, before and
after. A ratio below 1 means the model-form term is larger than the quoted interval.

## 3. What happens if G15.4 flips — decided now

The interval is **not widened** and `rmse_log` is not reinterpreted. The page gains a
sentence naming the model-form term as the larger one and saying that the quoted interval
is a parameter interval, not a total uncertainty. `results.md` already says this
qualitatively in its Reading point 1 (*"they bound the parameter uncertainty, not the
model-form uncertainty"*); what would change is that it becomes quantitative and states
which is bigger.

## 4. Predictions

1. G15.1 passes: the forces are unchanged and the intervals move. **If they do not move at
   all, item 15 closes as a no-op and that is a finding about which callers item 13
   reached.**
2. G15.2 passes at 6/6.
3. G15.3 finds **zero** open ends, so that branch remains untested. Recorded as untested
   rather than counted as verified.
4. **G15.4 does NOT flip on the 0.1 wt% gel** (±37 % against 16.6 % is a factor 2.2 of
   headroom and the corrected hull would have to shrink by more than that), **and DOES
   flip on the 0.5 wt% gel** (±2.6 % against 0.4 % is a factor 6.5, but that gel's
   intervals are the narrow ones and therefore the ones the bisection grid could most
   easily have inflated).

Prediction 4 is deliberately split so it cannot be read as right whichever way it comes
out. The figure's pinned md5 will change and must be updated with the run that changes it.
