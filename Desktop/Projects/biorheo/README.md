# biorheo

Cited fluid records, identifiable constitutive models, and their handoff to solvers,
for fluids that are not water. Phase 0 of the biofluids venture
(`../biofluids_venture/`, design spec and plan under `docs/superpowers/`).

The package exists to make three questions answerable from a flow curve before any
simulation is run:

1. Which constitutive parameters does this data actually determine? (`identify`)
2. Where did the numbers come from, and over what range were they measured? (`registry`)
3. What does a device calculation rest on when it leaves that range? (`apps`)

## Layout

```
src/biorheo/
  units.py            one pint registry; boundaries convert once, everything inside is SI
  constitutive/       RateModel (a dst_fit Model in shear-rate space) and the families:
                      gnf.py (newtonian, power_law, herschel_bulkley, casson,
                      carreau_yasuda, cross), structural.py (Moore-Coussot kinetics),
                      regularise.py (papanastasiou / bercovier_engelman / clamp, with the
                      knob named in the model and its induced error reported)
  identify/           fit_rate_curve -> IdentifiedModel: profile-likelihood intervals and a
                      per-parameter verdict (measurable / bounded_only / not_identifiable /
                      degenerate, the last from |corr| > 0.99 or a Fisher null direction);
                      synthetic.py recovery studies; record.py writes the sidecar
  registry/           YAML fluid records (parameters with units, measured envelope, DOI
                      sources, licence); identifiability lives only in a sidecar the
                      identify layer writes
  solver/             openfoam.py: export a record to an OpenFOAM v2606
                      `transportProperties` viscosity model -- six of the seven families
                      (HerschelBulkley, Casson, CrossPowerLaw, BirdCarreau, Newtonian,
                      powerLaw) -- parse it back, and refuse a record whose sidecar calls
                      a parameter degenerate or not identifiable without an explicit
                      override. `structural_steady` is refused with a named error: its
                      steady curve IS CrossPowerLaw with n = 1, but OpenFOAM's
                      viscosityModel cannot carry the structure equation, so the export
                      would be right in steady state and silently wrong in a transient
  apps/               injectability.py: Herschel-Bulkley needle flow (closed form for any
                      n, quadrature for any monotone model), plunger force with profile
                      and covariance intervals, log-sensitivities, and an extrapolation
                      flag against the fitted envelope
derivations/          SymPy scripts that regenerate their .md; each exits non-zero on a
                      failed identity or hand spot
scripts/              build the Carbopol records from the ETD2020 deposit and identify them
case_studies/         one script per study, writing results.md and a figure beside it
tests/                fast suite; `-m "not slow"` skips re-running dst_fit's own suite
```

## Install and test

```
pip install -e ../dst_fit -e .[dev]
python -m pytest -q tests/ -m "not slow"
python derivations/gnf_limits.py && python derivations/structural_kinetics.py && python derivations/buckingham_reiner.py && python derivations/couette_uniform_shear.py && python derivations/casson_tube_flow.py
```

`dst_fit` is consumed unchanged; `tests/test_gate_0_1_dstfit_unchanged.py` runs its
suite to prove it.

## What the data currently supports

Fifteen Carbopol records (doi:10.25835/hjn2qcvy, ETD2020, smooth plates, 20 C) with
measured envelopes of [1e-3, 200] 1/s and every Herschel-Bulkley parameter graded
`measurable`; 13 of the 15 are flagged `residual_structure.structured`, so the HB form
is a systematic misfit at rmse 0.2-17 % even where every parameter is determined. The
one case study so far, needle injectability of two of those gels, puts the needle wall
2.9 decades above the top of that envelope and finds the force ordering of the two gels
reverses between rheometer and needle; see
`case_studies/injectability_surrogates/results.md`.

Twenty human whole-blood records (Horner 2020 thesis workbook via Farrington,
Jariwala, Armstrong 2023, doi:10.1007/s00397-023-01402-2; 37 C, EDTA, 16 rates in
[0.1, 700] 1/s, one donor each) as Casson fits with `commercial_use: false` (author
posted, no licence stated). The authors' Casson table is reproduced to 1e-8 before
anything is built on it, and it is a convention: an unweighted line in sqrt space,
whose yield stress is 0.4-71 % (median 33 %) above biorheo's log-viscosity fit of the
same 16 points, because Casson misfits every donor but one (3 runs of residual sign
against 9 expected, under either objective). Cross and Carreau-Yasuda fit the same curves 5x better; the table
is `data/horner2020/model_comparison.md`. That table now also carries two NULL columns,
`newtonian` and `power_law`, fitted to the same curves and shipped as records so the L3
bridge has one of each to gate against. They are the floor the rest of the table is read
against, and they say two things worth knowing: Casson beats a one-parameter constant on
20/20 donors but beats a two-parameter power law on only **14/20**, at the same parameter
count; and the `residual_structure` flag does NOT separate a fit from a null here -- it
fires on 20/20 for both nulls and for Herschel-Bulkley, and the only two forms it clears
are Cross and Carreau-Yasuda. The workbook's 300 -> 1 1/s step sheet
(`data/horner2020/step_response.md`) shows the steady sheet's 1 1/s stress is
bit-identical to the 4-6 s peak of that transient on 14/20 donors (a later sample on
one more, within 1e-4 on four); the thesis confirms
this is the protocol (at or below 2 1/s the steady value is "the maximum of the
transient data", quoted in `SOURCE.md`), so the six lowest points of every blood record
are peak-hold values, 7-46 % above the 40-50 s stress at 1 1/s. Refitting with that one
point replaced moves tau_y only -1 to -7 %. Every donor's stress decays after the peak;
the thesis attributes the decay to syneresis (wall depletion), so it is not a bulk
signal for a thixotropic model to fit. `constitutive.structural` fitted to the rise
(`scripts/structural_step_horner.py`, `data/horner2020/structural_step.md`): the
timescale gate passes 20/20 (tau 1.1-2.2 s, 95 % approach at 3.4-6.6 s) and the shape
gate fails 6/20 -- a single exponential misses the rise by 1-12 % of its amplitude
(median 3.4 %) with the same residual sign pattern on 15/20, and anchoring the endpoints
to the steady sheet does not rescue it (5/20). The rise is less concave than one
exponential; one first-order structure variable cannot give that.

Whether that rise is a **lag or a second timescale** is not decidable from this workbook,
and the reason is not the one we expected (`scripts/oscillatory_capability_horner.py`,
`data/horner2020/oscillatory_capability.md`). Five criteria were fixed before any was
graded; four failed. The one that PASSED is the question we had written down -- the LAOS
sheets carry 20 to 1571 samples per cycle over 3 to 64 cycles, against a 4-6 sample rise,
so sampling was never the constraint. What decides it is that the workbook stores a
**storage modulus and nothing conjugate to it** -- no loss modulus, no phase angle, no
tan(delta) -- and a lag and a second timescale differ chiefly in phase; that **linearity
is never demonstrated**, the amplitude sweep starting at 9.95 % strain, at or above the
frequency sweep's own operating strain, with G' still rising there in both ramps on 9 of
15 donors; and that the **population is not the same** -- 12 of 20 donors appear on any
LAOS sheet and 8 (A, B, C, D, F, Q, R, T) on none, against a step-rise result on all 20.
A discriminating measurement needs a phase channel on the same donors, an amplitude sweep
that reaches into the linear region, or a step sampled an order of magnitude faster on
the rise.

The L3 solver bridge exports a record to an OpenFOAM v2606 `transportProperties` block.
The export itself is exact: OpenFOAM's own `nu` field and its `wallShearStress` function
object reproduce `record.eta_at` to 9.6e-13 relative, on 49 `simpleFoam` plane-Couette
cases spanning six records' whole measured envelopes -- one of each exportable family --
in a comparison that never touches
this repo's transcription of `calcNu` (`data/openfoam/couette_anchor.md`). Two of those
dictionaries are ones OpenFOAM must REFUSE, and does: a `powerLaw` without its required
`nuMin`, and a `Newtonian` with its `nu` moved into a `NewtonianCoeffs` block, which
OpenFOAM reports as `nu` missing from `transportProperties` itself -- i.e. it never looked
in the block. The useful
finding is what the export cannot carry. **A record's refusal to extrapolate does not
survive the export** -- `eta_at` raises outside the measured envelope, a dictionary
cannot, and every cell of a run started from rest sits decades below the lower bound.
And **the solve, not the export, is where a yield-stress hand-off breaks**: on a 20-cell
1-D problem `simpleFoam`'s outer iteration converges at a rate set by the local Bingham
number `B = tau_y/(k gd^n)`, returning a viscosity 172 % wrong at B = 29 with `rc = 0`,
no warning and a written result. That B and not the viscosity magnitude is the controller
was measured, not assumed, with a 2x2 that fixes the shear rate and moves `eta` 100x at
constant B (the error is unchanged to ten digits) against a 5.6e4x change across B. The
one diagnostic that needs no reference answer is the wall traction imbalance, which ranks
with the true error at Spearman 0.98 over 14 decades.

The needle app now has a Casson closed form beside the Herschel-Bulkley one
(`derivations/casson_tube_flow.md`), derived rather than quoted and anchored three ways:
against the package's independent quadrature (8.6e-12 on 24 cases), against a
velocity-profile route that carries the plug explicitly (7.2e-11 at x up to 0.93, and
unlike the quadrature it shares no step with the algebra), and against a hand-worked
point. Its bracket has a TRIPLE root at yield, so the flow rate leaves the no-flow branch
cubically in `1 - tau_y/tau_w` -- halving the distance to yield divides the flow rate by
eight, and the module evaluates the factored form because the expanded one goes negative
and is 6e2 times wrong within 1e-6 of yield.

Running the 20 blood donors through a nominal 27G needle at 1 mL/min
(`case_studies/blood_injectability/results.md`) gives two results worth more than the
force numbers. **All 20 are `extrapolated`, all by the same 1.42 decades**, because for
an effectively Newtonian fluid the wall rate is `4Q/(pi R^3)` regardless of viscosity:
how far a device extrapolates its own rheology is a property of the DEVICE, computable
at design time. (Scope: the decade is a device property; the third significant figure is
a property of the Casson form, which cannot shear-thin further above yield.) **The yield
stress is invisible here** -- `x = tau_y/tau_w` reaches only 1.7e-4, zeroing `tau_y`
moves the force under 2.9 %, and the force ranks with `eta_c` at Spearman 0.995 against
0.170 for `tau_y`; the parameter the published Casson tables disagree about by 71 % is
the one this device cannot see. And because `F ~ R^-4` (measured, -3.9723), **a 10.5 %
error in the needle's bore reproduces the entire donor-to-donor force spread** of the
17 donors the source study does not flag as atypical.

A claim was retracted from that study the day it was written: an ISO 9626 tolerance
comparison whose band mixed regular-wall and thin-wall figures and which reverses under
the narrower reading. The standard could not be checked offline, so the section was
rebuilt on the `R^-4` scaling, which needs no standard; the retraction, both readings and
both verdicts are on the page.

## The needle study carries uncertainty now, and it cannot rank patients

Phase 0's own exit gate asks this package for "full UQ from the L2 posterior". Pushing
each Casson parameter through its corrected 95 % profile interval gives a plunger-force
hull a median **29.1 %** of the force wide (covariance-sampled joint interval: 25.2 %),
against a donor-to-donor spread of 128 % and a family-choice effect of 12.5 %.

The useful question is not whether 29 % is smaller than 128 %. It is whether two given
patients can be told apart, and **they mostly cannot**: only **36 of 190** donor pairs
have disjoint intervals, and ranked by force two donors must be **17 places apart out of
20** to be reliably separated. The extremes are distinguishable; the middle of the
distribution is not.

That distinction cost a sentence. The page had said *"which patient the blood came from
still matters more than which family was fitted"* — true as a comparison of a range with
an effect size, and not a licence for the per-patient reading it invites. The range
comparison is kept, the pairwise one is withdrawn on the page with the old wording
quoted, and a test asserts the two sections cannot drift apart again.

`case_studies/blood_injectability/uq_prereg.md` (gate, committed before the code) and
`uq.md`.

## Every confidence interval in the registry was rebuilt on 2026-09-06

`dst_fit.profiling` sizes both its outward step floor and its ten bisections from a
parameter's declared BOX BOUND -- a guard-rail so the fitter cannot wander, not a
statement about the data. For a viscosity bounded at `1e6` the first step is 200 on a
parameter of value 0.0038, and the returned endpoint is a point on the grid
`best + m*span/5.12e6`: **201 of 296 closed upper ends in this registry sat on it**,
worst 1639x the best-fit value above it.

They were not confidence limits. Evaluating the profile F statistic at the shipped
endpoints gives a **median of 710 times `f_crit`** where a limit has `F = f_crit`
exactly. `biorheo.identify.profiling` sizes the step from the parameter's own Wald error
and locates the end with Brent's method on `F - f_crit`; `dst_fit` is consumed unchanged
and its own nuisance minimiser is reused, so only the WALK differs. The routine is
anchored against a closed form it does not touch -- a one-parameter fit in log space has
`SSR(u) - SSR_min = n (u - ybar)^2` exactly -- to 2.4e-9 on 20 donors.

Three things worth carrying:

- **Verdicts moved only toward more flagging**, and only in `carreau_yasuda`, on 4 of 20
  donors. No fit, rmse or verdict changed for the other five families, and no number in
  either case study moved.
- **On two records the profile is not a function of the parameter.** The nuisance
  re-minimisation lands on different local optima at neighbouring values, so `F` moves
  0.239 to 10.06 between points 0.1 % apart. A root finder returns a number from that
  sign change; it is a discontinuity, not a limit, and such an end is now reported OPEN
  with the rejected value and its `F` kept so the rejection is auditable.
- **A correct interval did not change the conclusion it was expected to.** The gate
  asking whether the case study's Cross-vs-Casson force difference is resolved at 95 %
  FAILED both arms: the bug inflated only the UPPER ends while that difference is signed
  downward, so exactly one donor of twenty flips. What the corrected intervals do say is
  useful and was not the question asked: the Casson record's own parameter uncertainty at
  a 27G needle is a median **29 %** of the force, against a **12.5 %** family-choice
  effect -- **the uncertainty is the larger of the two**, so a needle force from that
  study should not be quoted to better than about 30 %.

`data/profile_interval_prereg.md` (gate, committed before the code), `..._diagnosis.md`,
`..._gates.md` and `..._g139.md`.

## The same bug was 1639x in one place and 0.008 % in another

Worth stating together, because each number alone is misleading. The profile-interval
defect (above) moved the registry's endpoints by up to **1639x the best-fit value**, and
it moved the Carbopol needle study's forces by **0.0081 %**. Both are measurements of the
same bug.

Three things decide which you get, and none of them is a property of the fix:

- **How wide the parameter's declared box is.** The bad step floor is `box/5000`, so a
  viscosity bounded at `1e6` is exposed and a flow index bounded at O(1) is not. Measured
  on the two Carbopol gels: `k` (box `1e6`) moved 2.94 %, `tau_y` (box `1e5`) 1.87 %, and
  `n` (box O(1)) **0.0011 %**.
- **How many points the curve has.** The step is `max(wald/4, box/5000)`, so a large Wald
  error means the box is never consulted. These gels carry 54 points; the blood records
  carry 16.
- **Which parameter the downstream number is sensitive to.** The needle force here is set
  by `n` — `d ln F / d ln n` reaches 7.4 — which is exactly the parameter the defect could
  not touch.

The same lesson landed a second way one item earlier: correcting the intervals did not
change the blood study's family comparison either, because the bug inflated the UPPER ends
while that comparison binds on the lower ones. **A correction reaches a downstream number
only through the parameters, and the ends, that number actually uses.**

`case_studies/injectability_surrogates/uq_check.md`.

## Conventions

- Models return viscosity `eta(gd)` in Pa s; `fit_rate_curve` takes `eta = sigma/gd`.
- Parameters are SI floats internally. A record stores `value` and `unit`, never `si`.
- A flow-curve fit reports its envelope; anything computed outside it is flagged
  `extrapolated` and the flag is carried, not dropped.
- Every derivation has a hand-computed anchor that does not pass through the code
  under test.
