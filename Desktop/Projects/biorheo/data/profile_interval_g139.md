# G13.9 — is item 12's force disagreement resolved at 95 %?

Written by `scripts/profile_interval_g139.py`. Pre-registered in
`data/profile_interval_prereg.md` addendum A.2 (commit `eac96db`) before it was run.

| gate | claim | bar | result | verdict |
|---|---|---|---|---|
| G13.9a | the disagreement is NOT covered by Casson's own parameter uncertainty | Cross force outside the hull on >= 11/20 | **9/20** | **FAIL** |
| G13.9b | the OLD intervals could not have refuted anything | Cross force inside the old hull on >= 19/20 | **12/20** | **FAIL** |

## Both gates FAILED and the prediction was wrong. Here is why.

The pre-registration predicted 20/20 on both. Neither bar has been touched and
neither number is restated as a success.

- median |Cross - Casson| force difference: **12.46 %**
- median Casson push-through hull, CORRECTED intervals: **29.06 %** of the force
- median Casson push-through hull, dst_fit's intervals: **5.16e+03 %** of the force
- median covariance-SAMPLED (joint) interval: **25.23 %** of the force

**The reason both gates failed is the same one, and it is a sign error in my
reasoning, not in any code.** The bug inflated the UPPER interval ends only --
the lower ends were never on the bisection grid in the direction that matters,
because the step floor bites where the box bound is far away and the lower bound
of a viscosity is near zero. But the force difference is SIGNED: item 12 measured
Cross BELOW Casson on 16 of 20 donors (paired sign test, p = 0.0118). So the
comparison turns on the LOWER end of the hull, which the bug never corrupted.

Measured rather than argued: of the 20 donors, the corrected hull's binding end
is the LOWER one on **8** and the upper one on **1**; on dst_fit's own intervals it is lower on **8** and upper on **0**.

Two consequences, both of which cut against the reason this gate was written:

1. **G13.9b could not have come out the other way for the right reason.** The
   one-at-a-time hull on the old intervals CONTAINS the hull on the corrected
   ones -- narrowing an interval can only shrink a push-through range -- so
   `old_contains >= new_contains` always. The bar of 19/20 was reachable in
   principle, but only if the excluded donors were excluded from ABOVE, and they
   are not. This is the 'a control that shares an identity with its target'
   failure in a new costume, and it was not caught by writing the gate down.
2. **Item 13 did not change item 12's conclusion.** That is the honest headline
   here, and it is the opposite of what the pre-registration expected.

## What IS established, and it is the useful part

On **11 of 20** donors the Cross force lies inside the Casson record's
own one-at-a-time 95 % push-through, and on the sharper covariance-sampled
interval it lies inside on **11 of 20**.
The Casson uncertainty at this needle is a median **29.1 %** of the
force against a median family difference of **12.5 %** -- the two
effects are the same size, with parameter uncertainty the larger of them.

So the answer to *is the plateau disagreement resolved at 95 %* is **no, on a
majority of donors**, and item 12's caveat is strengthened rather than lifted:
a needle force from this case study should not be quoted to better than about
29 %, and the family choice is inside that.

The sampled interval is reported beside the pre-registered one, not instead of
it. It is narrower because it is a joint region rather than a box hull, but it
assumes the Wald normal -- which is the approximation the profile likelihood
exists because it distrusts. Neither is simply the better number.

## What this does NOT say

It does not say Cross and Casson agree. Item 12's mechanism is untouched: Casson
misfits these curves on 20/20 donors in `rmse_log` and `residual_structure` flags
19/20 against 2/20, and none of that is a statement about parameter uncertainty.
It does not say Cross is right either -- both families are read a median 1.42
decades past their own envelope at this needle.

## Per donor

| donor | Casson F (N) | Cross F (N) | diff | corrected hull | verdict | old hull | verdict |
|---|---|---|---|---|---|---|---|
| A | 0.2873 | 0.2006 | -30.18 % | [0.2426, 0.3386] | **excludes** | [0.2426, 14.7] | **excludes** |
| B | 0.2621 | 0.249 | -5.01 % | [0.2429, 0.2826] | contains | [0.2429, 14.69] | contains |
| C | 0.225 | 0.2311 | +2.72 % | [0.2095, 0.2414] | contains | [0.2095, 14.66] | contains |
| D | 0.2313 | 0.2641 | +14.19 % | [0.2106, 0.2536] | **excludes** | [0.2106, 14.67] | contains |
| F | 0.2231 | 0.189 | -15.28 % | [0.1964, 0.2529] | **excludes** | [0.1964, 14.64] | **excludes** |
| G | 0.2938 | 0.2294 | -21.93 % | [0.2326, 0.3679] | **excludes** | [0.2326, 14.7] | **excludes** |
| H | 0.301 | 0.2732 | -9.23 % | [0.2418, 0.3717] | contains | [0.2418, 14.73] | contains |
| I | 0.2855 | 0.2632 | -7.81 % | [0.2343, 0.3457] | contains | [0.2343, 14.71] | contains |
| J | 0.3348 | 0.3391 | +1.29 % | [0.2895, 0.3858] | contains | [0.2895, 14.76] | contains |
| K | 0.2752 | 0.2041 | -25.84 % | [0.2346, 0.3215] | **excludes** | [0.2346, 14.7] | **excludes** |
| L | 0.2962 | 0.2574 | -13.11 % | [0.26, 0.3365] | **excludes** | [0.26, 14.74] | **excludes** |
| M | 0.2721 | 0.273 | +0.34 % | [0.2182, 0.3364] | contains | [0.2182, 14.71] | contains |
| N | 0.5086 | 0.38 | -25.27 % | [0.4179, 0.6147] | **excludes** | [0.4179, 14.96] | **excludes** |
| O | 0.2575 | 0.2251 | -12.56 % | [0.2226, 0.2968] | contains | [0.2226, 14.68] | contains |
| P | 0.2935 | 0.2483 | -15.41 % | [0.2587, 0.3321] | **excludes** | [0.2587, 14.72] | **excludes** |
| Q | 0.2715 | 0.238 | -12.36 % | [0.2279, 0.3218] | contains | [0.2279, 14.7] | contains |
| R | 0.2741 | 0.2602 | -5.06 % | [0.2337, 0.3201] | contains | [0.2337, 14.71] | contains |
| S | 0.2907 | 0.2747 | -5.50 % | [0.2507, 0.3359] | contains | [0.2507, 14.73] | contains |
| T | 0.2651 | 0.2385 | -10.03 % | [0.2317, 0.3024] | contains | [0.2317, 14.69] | contains |
| U | 0.3031 | 0.2481 | -18.15 % | [0.2652, 0.3453] | **excludes** | [0.2652, 14.73] | **excludes** |
