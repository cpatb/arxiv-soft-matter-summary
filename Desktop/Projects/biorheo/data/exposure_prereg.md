# Item 16 pre-registration — which numbers were exposed to the interval defect, and how to tell

**Committed before the scan was run.** Bars are not moved; a bar that turns out to be the
wrong bar is recorded with its number and a stricter test is added beside it. A
pre-registration is never rewritten; corrections go in an addendum at the foot.

## 0. Item 15 shipped a FALSE statement about this, and it was found building item 16

`case_studies/injectability_surrogates/uq_check.md`, `README.md`, the Status card and
`NEXT_STEPS.md` item 15 all say some version of:

> `box_fallback` fires on none of them, which is the direct statement that the defect
> could not act here.

**That is wrong, and it conflates two different comparisons.**

- `step_rule` is a property of the **corrected** routine. It reports which of `wald/4` and
  `|best| * 1e-3` set the first step, and `box_fallback` fires only when **both are zero**.
- The **defect** is a property of the old routine, whose floor was `span/5000` with `span`
  from the box. It acted whenever `span/5000 > wald/4`.

Those are different tests, and `box_fallback` never firing says nothing about the second.
Measured on the six Carbopol parameters, the old box floor **won on five of six**:

| gel | parameter | `wald/4` | `span/5000` | box won? | inflation ratio |
|---|---|---|---|---|---|
| 0.5 wt% | `tau_y` | 1.139e-02 | 19.99 | yes | **1755** |
| 0.5 wt% | `k` | 1.699e-02 | 200 | yes | **11770** |
| 0.5 wt% | `n` | 2.879e-04 | 3.143e-04 | yes | **1.09** |
| 0.1 wt% | `tau_y` | 1.294e-02 | 20 | yes | **1546** |
| 0.1 wt% | `k` | 5.697e-02 | 200 | yes | **3510** |
| 0.1 wt% | `n` | 4.195e-03 | 2.711e-04 | **no** | 0.065 |

**The mechanism item 15 reported is right in its conclusion and wrong in its reason.** The
flow index really is the parameter the defect could not distort — but not because the box
was never consulted. It is because `n`'s box is `[.., 2.0]`, so `span/5000` is
`3e-4`, which is *the same size as* `wald/4` for that parameter. The floor was applied and
it happened to be the right size. For `k` it was 11 770 times too big.

**So the quantity that matters is not a boolean; it is the ratio**

```
inflation = (span / 5000) / (wald / 4)
```

and item 16's own suggested predictor — `step_rule` — is the wrong instrument, which is
the *fifth* premise in this repo to need checking before use.

The item-15 conclusions that DO survive unchanged: the parameter intervals moved up to
2.94 % while the forces moved 0.0081 %, `n` moved 0.0011 % against `k`'s 2.94 %, and the
force is set by `n`. Only the explanatory sentence changes, and it changes from a boolean
to a ratio that explains the three orders of magnitude rather than asserting them.

## 1. The gates

| id | claim | bar | what a failure means |
|---|---|---|---|
| **G16.1** | the grid signature implies the box floor won | **0** parameters with `on_grid = True` and `inflation <= 1` | the mechanism written on four pages of this repo is wrong |
| **G16.2** | the converse is NOT claimed | report the count of `inflation > 1` and `on_grid = False`; **expected to be large** | if it were near zero the ratio would be a two-way predictor, which §2 says it is not |
| **G16.3** | exposure is graded, not binary | report the distribution of `inflation` over the registry, and the correlation between `inflation` and the measured `hi` movement | a ratio that does not track the movement is not the mechanism |
| **G16.4** | the correction is surfaced where a reader will find it | `wald_stderr` is written into every sidecar so `inflation` is computable without refitting | a reader has to re-run a fit to answer "was my number exposed?" |

### G16.1 is the falsifiable one and it is ONE-DIRECTIONAL

The grid signature `best + m * span/5.12e6` can only arise if the crossing was bracketed on
the **first** step and that step was `span/5000`. So:

- `on_grid = True` **must** imply `inflation > 1`. That is G16.1, and a single counterexample
  refutes the mechanism.
- `inflation > 1` does **not** imply `on_grid = True`: if the crossing needed several steps,
  the bracket is `step0 * 1.18^j` wide and the endpoint lands on a *different* grid that
  `profile_interval_diagnosis.py` does not test for. Those cells are expected and are not
  failures. G16.2 exists to report them rather than let them be read as noise.

## 2. What this is NOT

**Not a general impact estimator.** `NEXT_STEPS` item 16 says so and it is repeated here:
the deliverable is one number per parameter that a reader can compute from the sidecar,
plus the measurement that it is the right number. No API for "how wrong was my result",
because that depends on which parameter the result is sensitive to and on which END of the
interval it binds — two things measured in items 15 and 13 respectively and neither of them
derivable from `inflation`.

## 3. Predictions

1. **G16.1 passes with zero counterexamples.** If it does not, four pages need correcting
   and the mechanism has to be rederived.
2. **G16.2 finds a large count** — most parameters have `inflation > 1` (every viscosity
   bounded at `1e6` does) while only 219 ends are on the detectable grid.
3. G16.3 finds `inflation` spans many decades and correlates positively but loosely with
   the measured movement, because movement also depends on how far the true crossing is.
4. G16.4 is mechanical.

Prediction 1 is the one that matters. It is written so that its failure is expensive and
visible.
