# The strain clock, measured

Written 2026-08-28. What gate G11 does, what it found, and what the Basilisk
drop-impact solver should do about it. Companion to `SOLVER_HANDOFF.md`, which
proposed this gate that morning, and to `DATA_INVENTORY.md` §2 on the deposit.

Reproduce: `python scripts/gate_g11.py` from the package root. Needs the
gitignored Lin deposit. Exit 1 is the result, not a regression.

---

## 1. What this measures that nothing else here does

Every other gate fits steady flow curves. At steady state the relaxation law's
left-hand side is zero, so the two clocks that govern every transient drop out of
the fit. They are not poorly determined by a flow curve; they are absent from it.

Richards et al. (2019) Eq. 4, which the solver integrates:

```
df/dt = -(gammadot/gamma_0 + 1/t_r) (f - f_ss)
```

At constant rate this integrates in strain, `gamma = gammadot t`, giving an
exponential approach with one constant per rate:

```
1/gamma_eff = 1/gamma_0 + 1/(gammadot t_r)
```

So `gamma_eff` is measurable one rate at a time, and the two clocks separate by
how it moves across rates. Note the sign the law predicts: `gamma_eff` **rises**
with rate toward `gamma_0`, because the time term matters less the faster you
shear.

Lin et al.'s panels 2b and 2d are the measurement. Reversing the flow destroys
the frictional network; the instantaneous viscosity climbs back over accumulated
strain, at four rates, on two materials. Inverting those records through the same
Wyart-Cates fit that G10 validated turns them into `gamma_0`.

## 2. The result that makes the number quotable

`f` comes from inverting Wyart-Cates, and the inversion is analytic:

```
eta_r = (1 - phi/phi_J)^-n,   phi_J = (1-f) phi_0 + f phi_m
  =>   f = (phi_0/phi - 1/u) / (phi_0/phi - phi_m/phi),   u = 1 - eta_r^(-1/n)
```

**Only the ratios `phi_0/phi` and `phi_m/phi` appear.** Those are exactly what
this project's central degeneracy preserves: `phi` and `phi_m` enter the model
only through `phi/phi_J`, so a flow curve pins the ratios and neither number, and
the fitted pair slides along that valley. The inversion is a function of the
ratios alone, so it is constant along the valley.

The invariance is in fact **stronger than that**, and the first version of this
section understated it while gating on something that could not fail. Both
corrections are due to the Drop Impact session, which reviewed this gate on the
night it was written; the algebra below was then re-derived and re-checked here
rather than taken on trust.

Write `u = phi_J/phi = 1/(1 - eta_r^(-1/n))`. Then

```
f = A + B u,   A = (phi_0/phi) / (phi_0/phi - phi_m/phi),   B = -1 / (phi_0/phi - phi_m/phi)
```

**`f` is affine in `u` for every packing triple.** Two different triples therefore
give two `f` records related by an affine map of the ordinate — and a relaxation
fit with *both* amplitudes free is exactly invariant under such a map, because
the amplitudes absorb `A` and `B` and the decay constant does not see them.

So `gamma_eff` does not depend on `phi`, `phi_m` or `phi_0` **at all**, jointly or
separately, and not merely along the degenerate ray. It depends only on `n`.
Checked over four triples spanning (0.20, 0.35, 0.62) to (0.55, 0.60, 0.63):
`f` is affine in `u` to a maximum residual of 5.7e-14 to 8.0e-13, and the fitted
`gamma_eff` on a synthetic record is **bit-identical** across all of them, not
merely close.

**What the original G11a actually tested, and why it was worthless.** It swept
`phi_0` over 0.58 to 0.68, refitting the steady curve each time, and reported a
worst spread in `gamma_eff` of 7.1e-8 on the silica and 4.6e-10 on the PMMA. Those
numbers are real but they are an identity, not a measurement: `Model.with_fixed`
moves the bounds with `phi_0`, the steady fit slides `phi` and `phi_m` along the
covariance ray (`phi/phi_0` held at 0.611341 and `phi_m/phi` at 1.299132 to six
digits across the whole sweep), `phi/phi_J` is unchanged pointwise, and so the
inverted `f` is the same array. **The gate could not have failed.** That is the
same failure mode as a gate that names the wrong quantity, caught here by review
rather than by the gate.

The replacement gates the two things that can move: the triple varied *off* the
ray, which must give exact invariance, and the exponent `n` over 1.8 to 2.5,
which is the one constant `gamma_eff` really does depend on.

**`gamma_0` inherits `sigma*`'s identifiability, not `phi_m`'s degeneracy** — and
by a wider margin than first claimed. This matters commercially as much as
scientifically. The package's headline is that `phi_m` is not identifiable from
one flow curve. G11 says which parameters escape that, and why, instead of
leaving the impression that nothing survives.

## 3. Wyart-Cates cannot reach the reversal dip at all

The model bounds the viscosity from below at `f = 0`: with no frictional contacts
the jamming point is the frictionless one, so `eta_r >= (1 - phi/phi_0)^-n`.

On **every record tried, the post-reversal minimum sits below that floor.**

| | floor eta_r | measured minima, as a fraction of it |
|---|---|---|
| 2 um silica | 6.620 | 0.743, 0.688, 0.646, 0.685 |
| 1.4 um PMMA | 19.662 | entirely below, 0.564, 0.638, 0.557 |

Lin's PMMA at `gammadot = 0.05` never climbs above the floor anywhere in its
record. `admissible_window` refuses it, and rightly; `fit_recovery` still returns
a clock for it, because the clock is a property of `u = phi_J/phi`, which is
defined and monotone for any `eta_r > 1`, and is invariant to the affine map that
turns `u` into `f`. So the record has a relaxation constant but no frictional
fraction: its inferred plateau in `f` is **-0.225**. It is reported and then
excluded from anything that tests the model — see §4.

This is not a fitting failure. Reversal does not only break contacts, it reverses
the microstructural anisotropy that carries hydrodynamic stress, and Wyart-Cates
has no state variable for microstructure. **The first ~0.7 to 1.3 strain units
after a reversal are outside the model's reach entirely**, whatever `f` does.
For the solver this is direct: a Wyart-Cates closure driven by an `f`-clock will
under-predict the softening after any flow reversal, and in a drop impact the
flow reverses at the rim on every impact.

## 4. Which clock governs, and a form that does not hold

Judged on a landmark that needs no fit: `gamma_50`, the strain from the reversal
minimum to half the plateau in `f`. `gamma_eff` is fitted over a window whose
*start* moves with rate, because a slower recovery crosses the floor later, so a
rate trend in it could in principle be a trend in the window. `gamma_50` cannot
carry that.

Over a 10x rate span, a pure strain clock gives `gamma_50` spread 1x and `t_50`
spread 10x; a pure time clock gives the reverse.

| | gamma_50 spread | t_50 spread | verdict |
|---|---|---|---|
| 2 um silica | 2.14x | 21.4x | strain |
| 1.4 um PMMA | **1.16x** | 11.6x | strain |

**G11b passes on both.** On the PMMA the strain collapse is nearly exact.

**G11c fails on both, but only one of the two failures means anything**, and
separating them was the last correction of the night.

On the **silica** the fitted clocks fall monotonically across all four rates —
1.0896, 0.8312, 0.5846, 0.5439, a span of 2.00x — and the two-clock line returns
`t_r = -1.82 s`. That is a real result: an additive time term cannot produce a
clock that *falls* with rate, at any `gamma_0` and any `t_r`.

On the **PMMA** it is not. Under the corrected estimator the clocks are
0.9690, 1.0896, 0.8775 across the three in-range rates: **not monotone**, and the
regression returns r² = -3.54, meaning a horizontal line would fit better. The
sign of `t_r` there is describing scatter. **On the PMMA the additive form is
untested, not contradicted**, and the gate now says so in those words rather than
counting it as a second refutation. The fit-free `gamma_50` is monotone on that
material (1.638, 1.536, 1.414) but over only three rates and a 1.16x span.

The first version of this file claimed the refutation held on both materials, on
the strength of the earlier estimator. It does not.

One record is excluded from G11c and from the quoted clock: Lin's PMMA at
`gammadot = 0.05`, whose inferred plateau in `f` is **-0.225**, outside [0, 1]
entirely. The rule is a property of the model rather than of the answer — a record
the model has no state to describe cannot be used to test one of the model's laws
— and it excludes exactly one record. Its clock is still reported, because the
clock is affine-invariant and is really a property of `u`, not of `f`.

This answers the question `SOLVER_HANDOFF.md` left open that morning. The
landmark version of this trend was visible in viscosity space, and the obvious
escape was that `eta` is strongly nonlinear in `f` near jamming, so the trend
might be an artefact of reading `eta` rather than `f`. It is not. Working in
`f`-space makes the trend **stronger**, not weaker.

### The interpretation has to be weaker than the gate outcome

A single exponential fitted to the whole transient is a **blend of two
processes**, so the `gamma_eff` fed to the two-clock fit is not the model's
`gamma_eff`. The Drop Impact session fitted four nested models jointly across all
four silica rates, eliminating the amplitudes exactly:

| model | par | SSR | per-rate clock |
|---|---|---|---|
| M1 one exponential per rate | 12 | 2.05e6 | 1.122, 0.837, 0.624, 0.592 — falls **1.90x** |
| M3 two exponentials, both clocks shared across rates | 14 | 1.06e5 | — |
| M4 shared slow clock, per-rate fast clock | 16 | 7.3e4 | 0.366, 0.366, 0.304, 0.215 — falls **1.70x** |

M3 beats M1 by 19x in SSR for two more parameters, and AICc prefers it decisively
even at an effective sample size of 100 taken from residual sign runs rather than
the nominal 570 bins. M4 then beats M3 by 32 AICc units.

So the honest sentence is: **the single-exponential clock falls 1.90x with rate;
correcting for a second component reduces that to 1.70x but does not remove it.**
An additive time clock still cannot describe what is left — the sign is wrong
either way — but "Richards Eq. 4 is refuted" overstates a test run on a
mis-specified quantity, and the two-component description changes the parameter
values substantially.

Reproduce with `python scripts/probe_two_clock.py`. Every number in the table
above came back identical when it was re-run here after being copied in.

**Do not read M3's two clocks as two physical processes.** They come out at 0.2205
and 0.3760 — close together, and with *opposite-sign* amplitudes, so that fit is
partly describing curvature rather than separating a contact clock from a
hydrodynamic one. A second fit on a different window gave 0.492 and 1.909 for the
same pair. A factor of five between two windows is the honest uncertainty on the
decomposition. What survives is the model *comparison* and the rate trend in M4's
fast clock.

The script also tests, and **refutes**, the natural follow-up guess that the fast
component's share tracks the steady-state frictional fraction. The share runs
0.465, 0.453, 0.405, 0.357 while `f_ss` runs 0.116, 0.411, 0.703, 0.885 — the two
move in *opposite* directions across the rate span.

Two things this does *not* establish. It does not say Richards' law is wrong on
Richards' own corn starch, which is a different material and is not tested here.
And it does not identify the right functional form; a multiplicative or
`f_ss`-dependent clock is one alternative among several, and the `f_ss`-dependent
one is now the less likely of the two.

## 5. What G11d caught, which is why it was worth pre-registering

**G11d failed on the silica under the original estimator**, moving `gamma_eff`
by 1.88x at the lowest rate and 1.64x at the next, against 1.11x and 1.04x at the
two highest. **Under the corrected estimator it passes on both**, worst 1.104x on
the silica and 1.371x on the PMMA. The gate did its job: it failed, the failure
was diagnosed, the estimator was fixed, and it now passes for a reason rather than
by having its tolerance widened.

The obvious reading — that the low-rate records had not converged — is wrong, and
checking it mattered. Over `gamma` 5 to 10 the inverted `f` drifts by +0.015 at
the low rates and +0.003 at the high ones, and each record's endpoint matches the
independently deposited steady-state file to between 0.2% and 3.5%. The records
are converged. The failure is in the *estimator*.

The specific cause, found by the Drop Impact session's review: the fit was
starting at `admissible_window.i_start`, the first point back above the
frictionless floor, which on the silica records is `gamma` = 0.70 to 1.26. That is
**at or past `gamma_eff` itself** (0.42 to 1.67), so the fit was looking almost
entirely at the plateau and the constant was being set by the plateau's shape
rather than by the recovery. Truncating at `gamma` = 3 then left barely one
e-fold of lever arm, which is why the low rates moved 1.9x and the high rates,
whose windows open earlier, did not.

The fix, now in `fit_recovery`, is to start from the binned reversal minimum
instead, log-bin the record so the plateau stops carrying 90% of the weight, and
solve the amplitudes exactly by variable projection while scanning the clock on a
log grid. That last part is not a style preference: the previous `least_squares`
call returned `gamma_eff` = 0.587 or 25.3 for the *same record* depending only on
the convergence tolerance, because the amplitudes and the log-time are scaled
decades apart. `admissible_window` is still computed and reported, because the
distance below the floor in §3 is a result in its own right.

So the fitted `gamma_eff` at the two lowest silica rates is window-sensitive under
the original estimator and should not be quoted from it. The fit-free `gamma_50`
is what the rate trend rests on, and it agrees.

## 6. The number the solver needs

Bracketed over all three estimators — the full-window fit, the truncated fit, and
the fit-free landmark:

| | gamma_0 | estimates | vs the solver's 0.2 |
|---|---|---|---|
| 2 um silica | 0.54 to 1.59 | 12 | 2.7x to 7.9x |
| 1.4 um PMMA | 0.71 to 1.64 | 9 | 3.5x to 8.2x |

**Every one of the 21 estimates on both materials exceeds 0.2**, and the two
materials' brackets overlap over 0.71 to 1.59. The claim to carry across is the
ordering, not a single value: the strain clock on these suspensions is of order
one strain unit, not two tenths.

The solver's `gamma_0 = 0.2` is borrowed from Richards' corn starch, and the sim
card's L9 row shows peak `jam_frac` = 0.912 / 0.997 / 0.671 at `gamma_0` =
0.05 / 0.2 / 0.8. A measured value near 1 sits **off the end of that sweep**, so
the sweep should be extended before anything is concluded from it. Corn starch
may genuinely differ from silica here; what is not defensible is running a
corn-starch number on silica while calling it the differentiating term.

## 7. Also done this session, and separate from G11

`assess()` fired MISSPECIFIED on fourteen of the fifteen real curves in this
project, which made it a statement about data quality rather than a verdict. The
runs test measures whether structure is *detectable*, and on a precise 25-point
sweep it nearly always is.

`systematic_misfit()` now reports how *big* the structure is: the rms of a cubic
through the ordered residuals, in log viscosity. The model is rejected only when
structure is both detectable and larger than `SYSTEMATIC_TOL = 0.05`, five per
cent, which is the reproducibility of a rheometer flow curve. That threshold was
set from the instrument before the audit was first run and has not been adjusted.

Across all 23 real curves on disk (`python scripts/audit_verdicts.py`): structure
detected on 15, large enough to reject the model on **5**, the other 10 reported
as the model being slightly wrong and the data being good enough to notice.

An external check nobody arranged: on Guy's d = 3770 nm series the two curves it
rejects are **56% and 58%**, which are exactly the two the paper itself reports as
unsteady with massive wall slip at `phi >= 0.56`. The threshold never saw that
sentence. This is a candidate answer to the open question about a pre-registered
steadiness criterion, since it is set by the instrument rather than inherited from
the source paper.

The blind spot, which is real and is not fixed. Rejection needs both a significant
runs test and a large amplitude, so a short sweep cannot be rejected at all. Guy's
404 nm 59% and 61% curves carry amplitudes of 0.077 and 0.059, above the
threshold, on 11 points, so the runs test only reaches p = 0.058 and the
conjunction never fires. The paper calls those unsteady too. A criterion for short
sweeps is an open question, and loosening `RUNS_ALPHA` after seeing that table
would be choosing the answer.

`Assessment.verdict()` also now runs `profile_likelihood` on request and says
"bounded below by X, unbounded above" where that is the truth, instead of quoting
a symmetric `+/-` the cost surface does not support.

## 8. A conjecture that died, recorded so it stays dead

On the silica, `gamma_50 * f_plateau` looked constant: four products of 0.418,
0.638, 0.674, 0.622, a spread of 1.61x against `gamma_50`'s own 2.14x. Tempting,
and the sort of thing worth an out-of-sample test before it becomes a claim.

On the PMMA the products are 0.924, 1.021, 1.203, a spread of 1.30x against
`gamma_50`'s own 1.16x. **The product is less constant than what it was meant to
improve on.** It was a coincidence of four points. The gate prints it either way.

## 9. What is left

- **`Fig3` is the control on the method and is not used.** The simulation records
  carry contact and hydrodynamic stress evaluated *directly*, so `f` can be read
  without inverting anything. Fitting the same clock both ways on the same record
  would test the inversion rather than assume it. `read_lin_transients` already
  reads `Fig_3a`; what is missing is the steady fit to pair it with, which
  `Fig_3c_total.csv` supplies at 26 rates.
- **What functional form does fit the rate dependence.** G11c rules one out and
  proposes nothing.
- **Corn starch.** Richards' own deposit has no reversal transients, so the 0.2
  cannot be checked on its own material from anything on disk.
- **A defect in `joint.py`, verified here tonight and not fixed.** `fit_series`
  frees `phi_0`, `phi_m` and a per-curve `phi` at once, and the model is
  homogeneous of degree one in all three, so the fit has an **exactly flat
  direction**. Confirmed two ways: scaling the triple by c changes every predicted
  viscosity by 2e-15, and the residual Jacobian's smallest singular value is
  1.5e-5 against a largest of 119, with a null vector proportional to
  `(phi_0, phi_m, 0, phi)` — the covariance ray, `sigma*` untouched. `np.linalg.inv`
  does not raise on that, which is why **N5 prints `phi_0 = 0.6912 +/- 924.4`**:
  those are not wide error bars, they are arithmetic noise. `n_params` also counts
  the flat parameter, so every joint AIC is 2 too high. The conclusions survive
  (the joint fit loses by 55 and 568 AIC), but the printed uncertainties should not
  be quoted. Fix: pin `phi_0` or `phi_m`, never both, and use `pinvh` with a rank
  check that raises.

**Closed the same night, by the Drop Impact session, and no longer open:**

- **`SOLVER_HANDOFF.md` §2's 179x.** Lin's solvent viscosities are *not* unstated
  after all. `Fig4a` is absolute viscosity and `Fig4b` is relative viscosity for
  the same sample, so `eta_s = eta_abs/eta_rel` gives **36.3 mPa.s for the 2 um
  silica** and **1.770 Pa.s for the 1.4 um PMMA**; `FigS3` measures the PMMA
  solvent directly at 1.765 Pa.s, agreeing to **0.29%**. The discrepancy was
  mostly my own error: reduced stress cancels `eta_s` only *within* a dataset, and
  §2 compared *across* datasets whose solvents differ by 7.66x while claiming the
  cancellation. In absolute units the gap is **24.6x, not 179x**, and
  179/7.66 = 24.6 exactly. §2 needs rewriting, and its first suspect, a wrong
  `eta_s`, is demoted: a 2.1x error there moves the fitted reduced `sigma*` by
  about 5%.
- **G12** is owned by the Drop Impact session, not by this one.
