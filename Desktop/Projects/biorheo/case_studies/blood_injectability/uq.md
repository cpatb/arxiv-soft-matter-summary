# Injection force with uncertainty, and whether donors can still be ranked

Written by `case_studies/blood_injectability/uq.py`. Gate pre-registered in
`uq_prereg.md` (commit `3d267db`) before any of this was computed.

| gate | claim | bar | result | verdict |
|---|---|---|---|---|
| G14.2 | the borrowed intervals belong to these forces | <= 1e-05 relative on 20/20 | worst 1.225e-06 | **PASS** |
| G14.3 | both modules describe the SAME device | exact equality | all four constants equal | **PASS** |
| G14.4 | an open interval is never printed as a number | 0 donors | 0 | **PASS** |
| G14.5 | donors can still be ranked at 95 % | >= 50 % of 190 pairs disjoint | **36/190 = 18.9 %** | **FAIL** |
| G14.6 | insensitivity reproduced by an independent route | `tau_y` insensitive on >= 19/20 | 20/20, worst |d ln F / d ln tau_y| = 1.48e-02 | **PASS** |

## The size of the uncertainty, and what kind of object it is

- one-at-a-time profile push-through hull: median **29.1 %** of the force
- covariance-sampled joint interval: median **25.2 %**
- donor-to-donor spread across all 20: **128 %**

**They are not interchangeable and neither is simply the better number.** The
hull moves one parameter at a time to each end of its own 95 % profile interval
and takes the force range, so it is WIDER than the joint region it is drawn from
-- the joint region is an ellipse inscribed in that box. The covariance sample IS
a joint region, but it assumes the Wald normal, which is the approximation the
profile likelihood exists because it distrusts. Every interval below says which
it is.

Every interval here is a CASSON interval. Reading the same donor through the
better-fitting Cross family moves the force by a further median 12.5 % (up to
30.2 %); see `crosscheck.md`. That is a separate term and it is not added in.

## G14.5 FAILED, and the ordinal claim is restated rather than deleted

**36 of 190** donor pairs have disjoint push-through hulls
(18.9 %, bar 50 %), and
39 of 190 on the sampled interval. The bar was not moved.

**What this does and does not overturn.** `results.md`'s arithmetic is
untouched: the donor spread is 128 % and the family-choice
effect is 12.5 %, so the spread is still the bigger term. What the arithmetic
does not license is the PAIRWISE reading of it -- 'which patient the blood came
from matters more' is a statement about telling two patients apart, and that is
a comparison of one donor's interval with another's, not of a range with an
effect size.

Ranked by force, two donors are separated at 95 % once they are
**17 places apart** in the order; nearer than that they
overlap. That is the honest form of the claim: **the extremes are
distinguishable, the middle of the donor distribution is not.**

## Which parameters the force is insensitive to

Phase 0's unit 0.6 asks for this list explicitly. It is a central-difference
`d ln F / d ln p` at the operating point, which is a **third and independent**
route to a claim `results.md` already makes two other ways (`x = tau_y/tau_w`
reaches only 1.7e-4, and zeroing `tau_y` moves the force under 2.9 %).

- `tau_y`: insensitive on **20 of 20** donors at the
  `INSENSITIVE = 0.05` threshold; worst
  `|d ln F / d ln tau_y|` = **1.48e-02**
- `eta_c`: worst `|d ln F / d ln eta_c|` = **0.995**, i.e. the
  force is very nearly proportional to the plateau viscosity and to nothing else

All three routes agree, which is worth stating because they could not have been
made to: the sensitivity is computed from the force solver by finite differences
and never touches the stress ratio the other two arguments are built on.

## Per donor

| donor | force (N) | push-through hull (N) | width | sampled joint (N) | width | open ends |
|---|---|---|---|---|---|---|
| F | 0.2231 | [0.1964, 0.2529] | 25.3 % | [0.1978, 0.2467] | 21.9 % | none |
| C | 0.2250 | [0.2095, 0.2414] | 14.2 % | [0.2108, 0.2389] | 12.5 % | none |
| D | 0.2313 | [0.2106, 0.2536] | 18.6 % | [0.2121, 0.2501] | 16.4 % | none |
| O | 0.2575 | [0.2226, 0.2968] | 28.8 % | [0.2250, 0.2886] | 24.7 % | none |
| B | 0.2621 | [0.2429, 0.2826] | 15.1 % | [0.2444, 0.2793] | 13.3 % | none |
| T | 0.2651 | [0.2317, 0.3024] | 26.7 % | [0.2340, 0.2956] | 23.2 % | none |
| Q | 0.2715 | [0.2279, 0.3218] | 34.6 % | [0.2304, 0.3116] | 29.9 % | none |
| M | 0.2721 | [0.2182, 0.3364] | 43.5 % | [0.2205, 0.3228] | 37.6 % | none |
| R | 0.2741 | [0.2337, 0.3201] | 31.5 % | [0.2362, 0.3112] | 27.4 % | none |
| K | 0.2752 | [0.2346, 0.3215] | 31.6 % | [0.2370, 0.3113] | 27.0 % | none |
| I | 0.2855 | [0.2343, 0.3457] | 39.0 % | [0.2369, 0.3323] | 33.4 % | none |
| A | 0.2873 | [0.2426, 0.3386] | 33.4 % | [0.2443, 0.3269] | 28.8 % | none |
| S | 0.2907 | [0.2507, 0.3359] | 29.3 % | [0.2533, 0.3274] | 25.5 % | none |
| P | 0.2935 | [0.2587, 0.3321] | 25.0 % | [0.2611, 0.3248] | 21.7 % | none |
| G | 0.2938 | [0.2326, 0.3679] | 46.0 % | [0.2348, 0.3488] | 38.8 % | none |
| L | 0.2962 | [0.2600, 0.3365] | 25.8 % | [0.2627, 0.3292] | 22.5 % | none |
| H | 0.3010 | [0.2418, 0.3717] | 43.1 % | [0.2444, 0.3558] | 37.0 % | none |
| U | 0.3031 | [0.2652, 0.3453] | 26.4 % | [0.2678, 0.3375] | 23.0 % | none |
| J | 0.3348 | [0.2895, 0.3858] | 28.8 % | [0.2923, 0.3758] | 25.0 % | none |
| N | 0.5086 | [0.4179, 0.6147] | 38.7 % | [0.4231, 0.5924] | 33.3 % | none |

## What is still not in these numbers

- The FAMILY term (median 12.5 %, worst 30.2 %), which is not added in above.
- The needle bore. `F ~ R^-4` measured, and a 10.5 % bore error reproduces the
  whole donor spread; the 27G inner diameter remains [UNVERIFIED].
- Any error in the Casson FORM itself. These are intervals on the parameters of a
  model that misfits these curves on 20/20 donors (`residual_structure` flags
  19/20). A correct interval on a wrong model is still a wrong uncertainty.
