# What survives when you use the authors' own steady criterion

Written 2026-08-29. Gate N8, the strategic item from `PLAN_2026-08-29.md`.
Reproduce: `python scripts/gate_n8.py`. Exit 1 is the result.

**Short version: the headline survives, but the dramatic version of it does not,
and the spec has to be rewritten around the weaker claim.**

---

## 1. The question, and the rule

The product's headline is that a single flow curve does not determine `phi_m`.
Part of the evidence has been that `phi_m` drifts monotonically across three
published phi-series, with reduced chi-squared of 221, 351 and 287 against a
single common value.

If that drift lives in curves the source labs themselves say are not steady flow
curves, the evidence is weaker than written. Two admission rules, both taken from
the sources rather than chosen here — which is the only thing that stops this
being a search for a convenient exclusion:

- **Richards et al. 2019.** Their deposited `Figure1/Plotting/Figure1.py` colours
  0.38 / 0.40 / 0.42 / 0.44 blue and **0.46 / 0.48 / 0.50 red**, and shades red
  over the region above the `phi_M` flow curve — the shear-jammed regime. The
  paper's subject is that the material *oscillates* there, so those three are
  time averages of an oscillating rate. The gate parses that colour list out of
  their script and fails if it ever stops saying so, rather than hard-coding it.
- **Guy, Hermes & Poon 2015.** The paper states the flow "was unsteady for
  phi >= 0.56 in both (a) and (c), and points represent temporal averages", with
  massive wall slip above `phi_m`.

The decision rule was committed in `PLAN_2026-08-29.md` before the gate existed:
if consistency largely recovers on the steady subsets, the claim is weaker than
currently written and the spec must say so. What was expected to survive
regardless was the degeneracy mechanism itself.

## 2. What happened

| dataset | all curves | admitted subset | reduced chi2 |
|---|---|---|---|
| Richards corn starch | 0 of 7 within 2σ, chi2 **221.1** | **3 of 4**, chi2 **11.0** | 20.1x better |
| Guy PMMA 404 nm | 0 of 7, chi2 **350.9** | **3 of 4**, chi2 **3.6** | 98.4x better |
| Guy PMMA 3770 nm | 2 of 7, chi2 **287.0** | **0 of 4**, chi2 **36.9** | 7.8x better |

The restriction is decisive on all three. **The kill criterion has fired.**

The `phi_m` ranges collapse with it. Richards goes from 0.4475–0.5070 to
**0.4475–0.4593**; Guy 404 nm from 0.5537–0.6018 to **0.5537–0.5631**, which
brackets that paper's own published `phi_m = 0.558(5)`.

**Guy's 404 nm series at reduced chi-squared 3.6 is close to consistent**, and it
is the clearest counterexample to the claim as previously stated. On the curves
its authors call steady, that material's `phi_m` is very nearly one number.

## 3. But it does not vanish, and that is the surviving claim

Reduced chi-squared of 11.0 and 36.9 are still far from 1, and two of the three
series still fail a 2σ cross-curve test on their admitted curves alone. `phi_m` is
determined to about ±0.001 on individual curves and still spans 0.009 to 0.013
across curves of the same material — ten times its own error bar.

So the honest claim is no longer *"`phi_m` drifts wildly across a series"*. It is:

> Even on the curves the source lab calls steady, independent fits of one material
> return `phi_m` values that no single value reconciles — by 10 to 30 times the
> quoted uncertainty. The dramatic version of that spread, reduced chi-squared in
> the hundreds, is substantially an artefact of including curves those labs
> already disowned.

That is weaker and much more defensible, and it is still a product: the tool's job
is exactly to catch a parameter that looks measured and is not.

## 4. N8c failed, and the reason matters more than the failure

N8c pre-registered that the degeneracy mechanism must survive the restriction:
worst `r(phi, phi_m)` above 0.95, smallest `phi_m - phi` below 0.02.

| dataset | worst r | smallest `phi_m - phi` | N8c |
|---|---|---|---|
| Richards corn starch | 0.948 | 0.0169 | FAIL (marginally, on r) |
| Guy PMMA 404 nm | 0.978 | 0.0120 | PASS |
| Guy PMMA 3770 nm | 0.944 | 0.0250 | FAIL (on both) |

Against the headline's `r` running to **+0.99** and `phi_m - phi` collapsing to
**0.0013**, the admitted curves give 0.944 to 0.978 and 0.012 to 0.025. High, but
an order of magnitude off the quoted extremes.

The mechanism is **present but weaker** everywhere, not absent. Richards misses by
0.002 in correlation.

**And the reason is close to tautological, which is the real lesson.** The
degeneracy grows as `phi` approaches `phi_m`. The unsteady curves are the
highest-`phi` curves. So excluding them excludes the most degenerate curves *by
construction*, and N8c could hardly have done anything else.

That is worth saying out loud rather than scoring: **the strongest evidence for
the degeneracy comes from exactly the curves that are the least trustworthy
measurements.** The two facts are not independent, and any future version of this
argument has to hold them together. A threshold I chose from the full-dataset
numbers, applied to a subset defined by removing the extreme end, is not a fair
test, and N8c should be rewritten before it is cited.

## 5. A defect this exposed, and what fixing it did

Guy's 3770 nm 19% curve fits to `phi_m = 0.3000 +/- 2.3768`. That is the
parameter's **lower bound**, with an error bar seven times the plausible range of
the quantity. At `phi = 0.10` the suspension barely thickens, so the curve carries
no information about `phi_m` at all, and the fitter returned a bound instead of
refusing.

The authors' unsteadiness rule does not exclude it, because it is not unsteady —
it is uninformative. So a second admission rule was added, and it is a property of
the fit rather than of the answer: **a fit that lands on a parameter bound is a
constrained optimum, not a measurement.** `FitResult.parameters_at_bound()` now
reports it and `assess()` says so first, before anything else in the verdict.

Getting the tolerance right mattered. A span-relative test in linear space calls
every ordinary `sigma_star` "at the lower bound", because its bounds run 1e-4 to
1e5 and a linear tolerance of 1e-4 of that span is 10. The check is done in log
space for any parameter whose bounds span more than two decades, linearly
otherwise. A test caught the first version.

**It moved the 3770 nm row the opposite way to the one predicted.** Dropping the
bad curve took that series from 1 of 5 to **0 of 4**, and reduced chi-squared from
27.6 to 36.9. The bound-hitting fit had a huge error bar, so it was dragging the
inverse-variance weighted mean toward itself and *masking* the disagreement among
the four good curves. Removing an uninformative curve made the inconsistency
clearer, not smaller — which strengthens §3's surviving claim on that dataset.

It also removed a number that had been quoted as evidence: that series' worst
`r(phi, phi_m)` was **1.000**, and it came from the bound-hitting fit. On the four
admitted curves it is 0.944.

## 6. What to do

1. **Rewrite the spec around §3's claim**, not the old one. This was the reason
   the spec rewrite was ordered after this gate.
2. ~~Add a bound-hitting admission rule.~~ **Done** — see §5. It moved the
   3770 nm row the other way from the one predicted.
3. **Rewrite N8c.** Its thresholds came from the full datasets and its subset is
   defined by removing the extreme end, so it cannot be a fair test as written.
   The honest version compares the degeneracy at *matched* `phi_m - phi`, or drops
   the pass/fail and reports the trend.
4. Re-examine whether the `phi_m` drift figure (`figures/fig_phim_drift.png`)
   should show the steady subset separately. It currently shows the full series
   and therefore the dramatic version.
