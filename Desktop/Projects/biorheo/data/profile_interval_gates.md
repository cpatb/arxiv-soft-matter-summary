# Item 13, step 2: the pre-registered gates

Written by `scripts/profile_interval_gates.py`. Bars are in
`data/profile_interval_prereg.md` and were committed before the fix existed.

| gate | claim | bar | result | verdict |
|---|---|---|---|---|
| G13.2 | the SHIPPED endpoint is not the crossing | `F >= 2 f_crit` on >= 90 % | 70/70 = 100.0 %, median `F/f_crit` = 711, worst 8.13e+03 | **PASS** |
| G13.3 | the corrected endpoint IS the crossing | `|F - f_crit|/f_crit <= 0.01` on 100 % of ends the routine CLOSES | worst 2.784e-04 over 325 profiles | **PASS** |
| G13.4 | matches a closed form that uses none of this code | <= 1e-06 relative | worst 2.378e-09 over 20 records | **PASS** |
| G13.5 | the correction can only NARROW | >= 95 % | 591/601 = 98.3 %, median narrowing 1x, worst 1227x | **PASS** |
| G13.8 | `dst_fit` untouched | clean tree, pinned figure md5 | `clean`, `e9d3c0430ce5b8fd09ba463901062218` | **PASS** |

## What G13.2 says

On the 70 shipped upper ends in the pre-registered sample, the profile F
statistic evaluated AT THE SHIPPED VALUE has a median of
**711 times `f_crit`** and a worst of
**8.13e+03**. A correct confidence limit has `F = f_crit`
there, by definition. This is the gate that could have downgraded item 13 to
cosmetic and did not.

## The ends the routine REFUSED to close, which G13.3 does not count

G13.3's population is *ends the routine closes*, as pre-registered. That wording
would be a hole if the routine could close fewer ends to make the number look
better, so what it declined is reported here: **5** of
325 profiles had an end rejected, at a worst
`|F - f_crit|/f_crit` of **0.951**
-- i.e. the rejected candidates were nowhere near being confidence limits, which
is the direction that says the check is catching real discontinuities rather than
discarding good endpoints.

A further **11** apparent crossings were refuted by a second nuisance
start during the walk (`false_crossings`): the warm-started solve had landed on a
worse local optimum, which is indistinguishable from a rise in the profile.

- `blood_horner2020_donor_B_carreau_yasuda` `carreau_yasuda` `eta_0`: lo -, hi [123.34063675622386, 0.736356122688396] (value, F there)
- `blood_horner2020_donor_B_carreau_yasuda` `carreau_yasuda` `eta_inf`: lo -, hi [0.0035870525595866187, 4.844335389068836] (value, F there)
- `blood_horner2020_donor_B_carreau_yasuda` `carreau_yasuda` `lam`: lo [66.57149055937975, 8.167243212861475], hi - (value, F there)
- `blood_horner2020_donor_B_carreau_yasuda` `carreau_yasuda` `n`: lo [0.35473851180231947, 4.847928100309653], hi - (value, F there)
- `blood_horner2020_donor_L_carreau_yasuda` `carreau_yasuda` `eta_0`: lo -, hi [198.5851785731262, 0.23905878896363014] (value, F there)

On these the profile is not a function of the parameter in any usable sense: the
nuisance re-minimisation lands on different local optima at neighbouring values,
so `F` moves 0.239 -> 10.06 between points 0.1 % apart. A sign change exists and
a root finder returns a number from it; that number is a discontinuity, not a
confidence limit. Both records already carry `degenerate` / `not_identifiable`
verdicts, which is the layer that was meant to catch exactly this.

## What G13.4 says, and why it is the load-bearing one

G13.3 grades this module with this module's own machinery: it re-evaluates the
statistic through the same nuisance minimiser the root finder used, so a
consistent error would pass it. The `newtonian` closed form does not go through
the outward walk or the root finder at all -- for a one-parameter fit in log
space the profile is exactly `SSR(u) - SSR_min = n (u - ybar)^2`, so the crossing
is `exp(ybar +/- s sqrt(f_crit/n))` in closed form. The routine reproduces it to
**2.378e-09** on all 20 donors.

## How much the intervals moved

Over the 601 ends both routines close, the corrected end is a median
**1x** nearer the best fit and at worst
**1227x**.

### G13.5 violations, reported per case rather than averaged away

- `blood_horner2020_donor_C_carreau_yasuda` `eta_0` (hi): old `269.925`, new `270.317`, best `0.0933417`
- `blood_horner2020_donor_D_carreau_yasuda` `a` (lo): old `9.99033`, new `1.38931`, best `10`
- `blood_horner2020_donor_F_carreau_yasuda` `lam` (lo): old `6.99592`, new `1.87965`, best `19.6724`
- `blood_horner2020_donor_K_carreau_yasuda` `eta_0` (lo): old `132.986`, new `0.045272`, best `133.116`
- `blood_horner2020_donor_K_carreau_yasuda` `n` (hi): old `0.417286`, new `0.539103`, best `0.362807`
- `blood_horner2020_donor_M_carreau_yasuda` `a` (lo): old `4.7849`, new `4.01044`, best `10`
- `blood_horner2020_donor_P_carreau_yasuda` `eta_0` (hi): old `436.27`, new `1685.98`, best `0.269652`
- `blood_horner2020_donor_U_carreau_yasuda` `eta_0` (hi): old `436.117`, new `1424.44`, best `0.116925`
- `blood_horner2020_donor_U_carreau_yasuda` `lam` (lo): old `24.3593`, new `9.02379`, best `24.3826`
- `blood_horner2020_donor_U_carreau_yasuda` `n` (lo): old `0.289491`, new `0.266067`, best `0.448022`

A widening means the two walks followed different branches of the profile:
`dst_fit`'s first step is large enough to jump over a local feature that the
finer walk resolves. It is a difference in the WALK, not in the objective --
both routines call the same nuisance minimiser.

## What these gates do NOT establish

Nothing about whether the MODEL is right. A profile interval describes the cost
surface of one model on one curve; a correct interval on a misfitting model is
still a wrong uncertainty, and `residual_structure` is the field that says so
(Casson is flagged `structured` on 19 of 20 donors).

They also say nothing about items 9 and 10, which are re-graded by re-running
`scripts/multistart_seed_scatter.py` and `..._carbopol.py` against the corrected
intervals with their bars unchanged -- gate G13.6.
