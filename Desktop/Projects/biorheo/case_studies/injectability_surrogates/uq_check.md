# Item 15: did the item-13 interval correction reach this study?

Written by `uq_check.py`. Gate pre-registered in `uq_prereg.md` (commit
`850a647`) before anything was recomputed.

**The parameter intervals moved and the FORCE did not.** The worst interval
endpoint moves **2.94 %**, so G15.1 passes its
1 % bar -- but the force hull moves at most
**0.0081 %**, which is why the shipped
table barely changed. The reason is the whole finding and it is in the next
section.

| gate | claim | bar | result | verdict |
|---|---|---|---|---|
| G15.1 | the correction reaches this page | an endpoint moves > 1 % | worst 2.9449 % | **PASS** |
| G15.2 | the corrected intervals NARROW | >= 2 of 2 | 2/2 | **PASS** |
| G15.3 | no end is silently open | reported | 0 open ends | **PASS** |

## The correction moved the parameters the force does not depend on

Item 13's defect is a step floor of `span/5000` taken from the parameter's
declared BOX bound, so **its size scales with how wide that box is**. In
Herschel-Bulkley the three boxes are not alike: `tau_y` is bounded at `1e5` and
`k` at `1e6`, while the flow index `n` is bounded at O(1). So the floor can
inflate the first two and essentially cannot touch the third.

That is exactly what happened, and the force is set by the one it could not
touch: `d ln F / d ln n` is 4.87 and 7.44 for these two gels, and `n` is the
widest contributor to the push-through hull on BOTH of them.

| gel | parameter | value | dst_fit interval | corrected interval | worst end move | step rule |
|---|---|---|---|---|---|---|
| `carbopol_etd2020_0p5wt_50uL` | `tau_y` | 41.1114 | [41.0194, 41.209] | [41.0194, 41.203] | 0.0145 % | `value` |
| `carbopol_etd2020_0p5wt_50uL` | `k` | 13.7321 | [13.5954, 13.9275] | [13.5954, 13.8697] | 0.4143 % | `wald` |
| `carbopol_etd2020_0p5wt_50uL` | `n` | 0.428261 | [0.425939, 0.430591] | [0.425939, 0.43059] | 0.0001 % | `value` |
| `carbopol_etd2020_0p1wt_40uL` | `tau_y` | 0.838648 | [0.74338, 0.955835] | [0.743403, 0.937988] | 1.8671 % | `wald` |
| `carbopol_etd2020_0p1wt_40uL` | `k` | 4.45963 | [4.04351, 5.04556] | [4.04354, 4.89698] | 2.9449 % | `wald` |
| `carbopol_etd2020_0p1wt_40uL` | `n` | 0.644484 | [0.613303, 0.676593] | [0.61331, 0.676592] | 0.0011 % | `wald` |

Read the `worst end move` column against the box bound: `k` (box `1e6`) moves
2.94 % and 0.41 %, `tau_y` (box `1e5`) moves 1.87 % and 0.01 %, and `n` (box
O(1)) moves **0.0011 % and 0.0001 %** -- three orders of magnitude less. The
`step rule` column says which term set the first step: `wald` on 4 of 6 parameters and the 0.1 %-of-value floor on 2. **Neither is the box**
-- `box_fallback` fires on none of them, which is the direct statement that the
defect could not act here. These curves carry 54 points against the blood
records' 16, so the Wald errors are large.

**This is the same shape as item 13's own G13.9 failure**, one item later: a
correction reaches a downstream number only through the parameters that number
is sensitive to. There it was the wrong END of an interval; here it is the wrong
PARAMETER. Neither is a defect in the fix, and neither was predicted.

## G15.4: is the parameter interval still the dominant error term?

`rmse_log` is the fit's own relative viscosity error and `d ln F / d ln k` is
0.979 and 1.000 here, so a viscosity misfit of `rmse_log` propagates to about
`rmse_log` in force and the two are comparable as percentages.

| gel | force at the 10 s dose | hull half-width | `rmse_log` | ratio | which dominates |
|---|---|---|---|---|---|
| `carbopol_etd2020_0p5wt_50uL` | 20.240 N | ±2.64 % | 0.35 % | 7.53 | parameter interval |
| `carbopol_etd2020_0p1wt_40uL` | 78.545 N | ±37.53 % | 16.58 % | 2.26 | parameter interval |

**Nothing flipped** (no gel changed side). The
pre-registration's prediction was deliberately SPLIT -- flip on the 0.5 wt% gel,
not on the 0.1 wt% -- and **the half that predicted a flip is refuted**. The
reasoning behind it was that the 0.5 wt% gel's intervals are the narrow ones and
therefore the ones a bisection grid could most easily have inflated. That is
backwards: the grid's size is set by the parameter's BOX, not by how narrow the
interval is, and this gel's hull is set by `n`, whose box is O(1).

The comparison is still worth having on the page, because it is the first time
these two error terms have been put in the same units here: on the 0.1 wt% gel
the quoted interval is only about twice the model's own misfit, so the error bar
and the form error are the same order and the interval is not a total error.

## What this says about item 13's blast radius

The defect was real and large where it bit -- 203 of 296 closed upper ends in the
registry, worst 1639x the best-fit value -- and it moved this page's forces by
**less than a hundredth of a per cent**. Both statements are about the same bug.
Two things decide whether a given number is affected, and neither is about the
fix: **how wide the parameter's declared box is** (the floor is `box/5000`), and
**how many points the curve has** (a large Wald error means `wald/4` wins the
`max` and the box never enters). A quantity dominated by an exponent bounded at
O(1), fitted to 54 points, is about as immune as it is possible to be.

**So `[UNVERIFIED]` is the wrong label for an untouched number here, and so is
`corrected`.** The right one is that the two routines agree on this study to
0.0081 %, which is a measurement and
is now on the page.

The two gates named UNFALSIFIABLE in the pre-registration are reported and not
counted: the ordering reversal survives (the hulls were already disjoint by a
factor of three against a 3.9x effect) and the nominal forces are unchanged.

## No open ends, so that branch is still untested

All parameters of both gels close at both ends. The sibling study's
"an open interval is never printed as a number" branch has therefore still
never executed. It is recorded as untested code, not as verified behaviour.
