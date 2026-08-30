# onset_clocks — what arrests an impacting suspension drop?

Started 2026-08-30. An elimination audit over the drop-impact programme's own
75-impact dataset, testing every candidate criterion the field has on the table
for why a dense-suspension drop stops spreading.

```
python run_audit.py                       # the elimination audit + its figure
python run_inversion.py                   # what phi_m the footprints would need
python -m pytest test_onset_clocks.py -q  # 27 passed
```

---

## 1. The contradiction this starts from

`impact_stress/materials.py` carries this as a caveat on the lab's own
calibration of its own silica, and it has never been resolved:

> At phi = 0.44 < phiM = 0.5176 the suspension can NEVER jam by this model, at
> any stress; it can only thicken continuously. The drops in the experiment
> demonstrably do arrest. That mismatch is real.

Every packing in the dataset is on the wrong side of it. The frictional fraction
`f` a drop would need for `phi_J(f)` to reach its own packing:

| phi | f required | can it jam? | n drops |
|---|---|---|---|
| 0.40 | 1.64 | **never** | 7 |
| 0.44 | 1.43 | **never** | 53 |
| 0.48 | 1.21 | **never** | 10 |
| 0.50 | 1.10 | **never** | 3 |

`f` cannot exceed 1. So on the lab's own fitted closure, none of these drops can
shear-jam at any stress, and all of them arrest.

Two escapes have been offered. **The first is that the rheometer's
volume-fraction axis is offset from the drops' "by an unquantified amount".**
That offset was quantified the same night, on two published datasets, in
`tier0_wyart_cates_closure/CALIBRATION_2026-08-30.md`: the inferred `phi` runs
*above* the true one, so correcting for it moves `phi_m` further from `phi`, not
closer. It makes the contradiction worse. **The second is that the impact stress
is so large the material thickens into arrest anyway.** That is what this audit
tests, along with everything else on the table.

## 2. The method: eliminate by degeneracy, not by correlation

With 75 drops and a strong velocity axis, almost any monotone function of `U0`
correlates with the response. Correlation cannot choose between the candidates.
Two things can.

**The exponent a candidate would need.** If the arrested footprint `beta_p` is a
function of a group `G` alone, then over the dataset
`d ln beta_p / d ln G ~ ln(range of beta_p) / ln(range of G)`. A group that
moves 1% while the response moves 6.4x is being asked to do work its own dynamic
range cannot do.

**Degeneracy.** If `G` sets `beta_p`, drops with the same `G` must have the same
`beta_p`. Among pairs whose `G` agree to 5%, the median `|d ln beta_p|` is
divided by the median over *all* pairs. Near 0 means `G` determines the
response; near 1 means knowing `G` has told you nothing.

Pre-registered rule: **FALSIFIED** if the 5% ties retain more than half;
**STRAINED** if the required exponent exceeds 10.

One honest note on the statistic: the strict falsifier is the *maximum* over
tied pairs, since one counterexample refutes a deterministic law. Over thousands
of pairs that is a single-pair statistic and saturates for nearly every
candidate, so the verdict uses the median instead. That choice was made after
seeing the maxima; both are printed.

## 3. The response, and the data

`Github/Drop Impact/analysis_2026-07-25/c3_master_curve_per_drop.csv`, 75 rows,
73 used (two have a bounded rather than a point `beta_p`). `beta_p` is the
arrested contact footprint over the drop diameter — a contact width, so it
starts near zero, not at one.

- phi 0.40, 0.44, 0.48, 0.50 · LiCl 0, 10, 50, 100 mM
- `U0` 0.75 to 2.59 m/s (3.5x) · `d0` 1.79 to 2.56 mm
- **`beta_p` 0.309 to 1.962, a factor of 6.35.** This is what a candidate has to
  explain.
- arrest, by the repo's own `beta_p < 0.9 beta_eq`, recomputed rather than read
  from a stored flag: 63 of 73.

**Is `beta_p` resolved well enough to carry a 6.35x claim?** Yes, measured two
ways. `beta_p_lo`/`beta_p_hi` are **not** error bars — they are a field-of-view
censoring bracket, set equal to `beta_p` for 53 of the 73 rows, so they say
nothing about uncertainty. The two real routes:

| route | n | scatter |
|---|---|---|
| replicate drops within a (phi, salt, nominal U0) cell, 12 cells with n >= 3 | 63 | median CV **6.7%** |
| the independent under-view channel against the side view | 21 | log sd **20%**, correlation 0.96 |

The response's log-range is 1.849, so it is **28x** the replicate noise and 9x
the cross-channel disagreement. The audit's degeneracy test is not resolving
noise.

## 4. What the closure's own viscosity does

Evaluate the model the drop-impact repo actually calibrated, at the stress each
impact actually delivers. Over the whole dataset its viscosity moves 42.7x —
but that is almost all `phi`, which the campaign varied deliberately. At the
single packing `phi` = 0.44, over 53 drops and a 2.26x range of impact speed:

- **the closure's viscosity moves 1.0076x**
- **the measured footprint moves 3.71x**

In log terms the response outruns the closure's stress response by a factor of
**174**.

## 5. The audit

Sixteen candidates. Ten falsified, five survive, and the survivor that orders
the data is not a stress clock.

| candidate | range | exponent needed | ties@5% retain | rho | verdict |
|---|---|---|---|---|---|
| `f` at the deviatoric gauge | 1.04x | 43 | 1.00 | −0.25 | **FALSIFIED** |
| `f` at `rho k U0^2`, k = 2 | 1.01x | 295 | 1.00 | +0.07 | **FALSIFIED** |
| `f` at `rho k U0^2`, k = 12 | 1.00x | 3533 | 1.00 | +0.07 | **FALSIFIED** |
| `f`, no added mass | 1.04x | 43 | 1.00 | +0.07 | **FALSIFIED** |
| `f`, `sigma*` x7.53 for salt | 1.11x | 18 | 0.97 | +0.07 | **FALSIFIED** |
| `Sigma/sigma*` | 11.5x | 0.76 | 0.77 | +0.07 | **FALSIFIED** |
| accumulated strain, `U0/D0` | 1.00x | — | 1.00 | −0.23 | **FALSIFIED** |
| `gammadot = U0/D0` | 3.42x | 1.50 | 0.75 | +0.17 | **FALSIFIED** |
| We | 12.2x | 0.74 | 0.74 | +0.02 | **FALSIFIED** |
| St, Pe (particle) | 3.42x | 1.50 | 0.75 | +0.17 | **FALSIFIED** |
| **jamming deficit** `(phi0-phi)/(phi0-phiM)` | 1.5x | 4.56 | 0.43 | **+0.708** | survives |
| accumulated strain on `L_c` | 6.04x | 1.03 | 0.31 | +0.29 | survives |
| Re | 7.17x | 0.94 | 0.30 | +0.19 | survives |
| `<tau>/sigma*` (deviatoric) | 4.57x | 1.22 | 0.48 | −0.25 | survives |
| `gammadot = U0/L_c` (Mobaseri) | 14.9x | 0.68 | 0.33 | +0.20 | survives |
| `f` deviatoric **with sigma* x7.53 for salt** | 2.01x | 2.65 | 0.60 | −0.25 | **FALSIFIED** |
| `kappa^-1/a` (Debye) | 100x | 0.40 | 1.05 | +0.09 | **FALSIFIED** |
| **pore Deborah** (contact length) | 3.54x | 1.46 | 0.37 | **−0.517** | survives |

The last three were added after the first pass, because nothing on the original
list could carry both the velocity and the salt axis.

- **The one live corner of the stress clock is now closed.** The first pass
  never tested the deviatoric gauge *combined* with a salt-shifted `sigma_star`,
  and at the 7.53x shift the programme's own collapse hunt inverted, `f` does
  acquire a 2.01x range — no longer flat enough to dismiss on range alone. It
  fails the degeneracy test instead (ties retain 0.60).
- **`kappa^-1/a` fails informatively**: its ties retain 1.05, i.e. *more* than a
  random pair. It takes one value per salt level, so its "ties" are "same salt",
  and within a salt level the footprint varies fully. Salt alone does not set it.
- **The pore Deborah number is the second-strongest survivor** (rho = −0.517,
  p = 2.8e-6). One correction to it, made after an over-claim: **it does not
  carry the salt axis.** The permeability is Carman-Kozeny, purely geometric;
  salt would enter through an electroviscous term or the drained modulus and
  neither is in there. What it carries is `U0` and `phi`.

## 6. The answer

**Every stress-activation clock is falsified, on both gauges.** The frictional
fraction `f` moves by 1 to 4% across the whole dataset while the footprint moves
by 6.35x. That is not a marginal call: `f` would need an exponent of 43 on the
deviatoric gauge and 295 on the inertial one. Even granting salt the largest
`sigma_star` shift the programme's own collapse hunt could invert out of the
data — a factor of 7.53 — `f` moves 11% and stays falsified.

**The accumulated-strain criterion on the field's standard shear rate is
falsified for a structural reason worth stating.** With `gammadot = U0/D0` the
strain accumulated over an impact time `D0/U0` is *identically 1*, for every
drop, at every velocity. A strain criterion built on the convention this lab's
own `gammadot_c = 632 1/s` is quoted on has no velocity dependence at all and
cannot produce a velocity onset. On Mobaseri's validated lamella length it does
vary, and it survives.

**The one candidate that both survives and orders the data strongly is a pure
packing group with no stress, rate, strain or time in it** — the jamming
deficit, at rho = +0.708, p = 2.5e-12, an order of magnitude better than
anything else on the list.

And that is also the audit's honest limit. The jamming deficit is a monotone
relabelling of `phi`, so it carries no velocity and no salt dependence — and
both of those axes move the data. So:

> No candidate on the field's table explains the velocity axis. The only one
> that explains anything explains the packing axis and nothing else. The
> stress-activated closure has no mechanism for the two knobs this experiment
> actually turns within a single material.

## 7. What this does and does not license

**Does.** Stop putting `f(sigma)` at the centre of an impact-arrest model for
this material; it is saturated and flat over the whole accessible range whatever
gauge is used. Stop quoting accumulated strain on `U0/D0`; it is a constant. If
a rate convention is used at all, use Mobaseri's `L_c`, which is the only one
that both survives here and was validated against a measured dissipation
integral rather than asserted.

**Does not.** This audit does not identify a mechanism. Four groups survive and
none of them orders the data well; the strongest survivor is a relabelling of a
knob. Ruling out the incumbents is not the same as finding the answer, and this
document does not claim to have found it.

## 8. What is untested here

- **`beta_p` may be the wrong response.** It is a contact width measured in side
  view, and the campaign's own analysis treats it as a *predictor* of the
  relaxation time rather than as a thing to be predicted. Nothing here checks
  whether it is resolved well enough to carry a 6.35x claim.
- **Salt enters no candidate at all.** Its only route into this model is through
  `sigma_star`, and the stress clock is falsified by dynamic range regardless of
  `sigma_star`. So salt's effect, which is large in the data, is unmodelled
  rather than tested.
- **The material constants are one calibration**, `Github/Drop Impact/closure/
  params_fitted.yaml`, which failed its own acceptance gate (rms(log eta) 0.374
  against a required 0.35) and whose `phi` axis is offset by the amount §1
  discusses. Every number here inherits that.
- **`d0` varies by 1.4x across the dataset** but was not controlled, so it is a
  nuisance variable in every group here rather than an axis.
- The Mobaseri lamella length is **undefined for the lowest-velocity,
  highest-`phi` drops**, where Re falls to order 1 and the correlation returns a
  lamella thicker than the drop. Those drops are reported as undefined rather
  than extrapolated, which is itself a finding: the field's best-validated
  shear-rate convention does not cover this experiment's own regime.

---

## 9. The inversion, and why its headline was withdrawn

`run_inversion.py`. The idea: instead of testing candidates, take the closure
the repo already fitted, treat `phi_m` as the single unknown, and ask what value
each drop needs to reproduce **its own** footprint. `phi_m` is a material
constant, so if the required value holds still the closure is fine and the
fitted 0.5176 was simply wrong.

It appeared to work. The required `phi_m` came back at a median of **0.441**
against the fitted 0.5176, constant to ~2x the volume-fraction tolerance across
salt and across velocity, and insensitive to the unknown prefactor over four
decades. That looked like a measured replacement value and a pre-registered
bench prediction.

### It is a tautology, and here is the check that shows it

If the method measures a material constant, the answer must not depend on which
packing the drop was made at. It does:

| `phi` | n | median required `phi_m` | minus `phi` |
|---|---|---|---|
| 0.40 | 6 | 0.4249 | +0.0249 |
| 0.44 | 44 | 0.4410 | **+0.0010** |
| 0.48 | 10 | 0.4782 | **−0.0018** |
| 0.50 | 3 | 0.4955 | −0.0045 |

**The required `phi_m` tracks `phi` to within 0.005 over three of four
packings.** The inversion is pinned: an arrested drop implies a large viscosity,
a large viscosity forces `phi_J` up against `phi` because the divergence is
steep, and `phi_m` is then read off as approximately `phi`. The prefactor
insensitivity that looked like robustness is the *same* pinning seen from
another angle.

So the "median required `phi_m` = 0.441" was the median of a dataset in which
44 of 63 drops are at `phi` = 0.44. It is `phi`, not `phi_m`.

**Withdrawn:** the replacement value 0.441, the claim that it is constant across
salt and velocity (those constancies are also pinned — within a cell they say
only that `beta_p ~ U0^(1/p)` at fixed material, which is the viscous spreading
law holding, not a statement about `phi_m`), and the bench prediction built on
them. The `phi_m` in §1 remains unmeasured by anything here.

### What does survive

1. **The viscosity gap, which needs no inversion.** The closure's ceiling at
   `phi` = 0.44 with `f` = 1 is **0.211 Pa.s**. The footprints imply 0.36 to
   1321 Pa.s. **63 of 63 arrested drops are above the ceiling**, by up to 1000x.
   That is the §1 contradiction restated in the units a solver uses, and it does
   not depend on the pinning: it is a comparison of two viscosities, not an
   inference of a packing.
2. **The arrested state is marginally jammed at its own local packing.** That is
   what the pinning is telling us, read forwards instead of backwards:
   `phi_J` at arrest sits within 0.005 of `phi`, at every packing from 0.40 to
   0.50. The drops do not overshoot into a deeply jammed state and they do not
   stop short — they arrest at marginal stability. That is a physical statement
   and it is the one thing this analysis genuinely adds.
3. **Nothing above 0.44 fully jams.** A true `phi_m` below 0.48 would send
   `eta` to infinity and `beta_p` to zero for the 0.48 and 0.50 drops. Observed
   `beta_p` there is 0.31 to 0.48, finite and well resolved. So there is no
   single `phi_m` in this range at which the material simply locks.

## 10. What salt actually does here

Worth stating because it removes a worry §6 raised. At `phi` = 0.44, among
drops that arrested, at matched velocity around 2.3 m/s, the footprint across
0 to 100 mM varies by only **1.28x** and non-monotonically, against a 6.35x
total range. At **1.15 m/s** it varies by 2.1x — but that is the onset moving:
at 0 mM the drops have not arrested yet, at 100 mM they have.

> **Salt shifts the onset velocity. It does not change the arrested state once
> arrested.**

So "no candidate carries the salt axis" was the wrong worry — in *this*
observable there is very little salt axis to carry. Salt's large effect in this
programme is on the un-jamming time `tau_u` (a factor 5.8x at 100 mM), which is
a different clock and a different paper. That is the programme's own two-clocks
picture, confirmed from the arrest side.

## 11. Where this leaves it

**Settled, no bench needed:**
1. The response is sound — 6.7% replicate noise against a 6.35x range.
2. Every stress-activation clock is falsified, including the salt-shifted
   deviatoric corner that the first pass had not tested.
3. The accumulated-strain criterion on `U0/D0` is identically 1 and cannot
   produce a velocity onset at all.
4. Every arrested drop demands 1.7x to 1000x more viscosity than the calibrated
   closure can produce at any stress.
5. Drops arrest at marginal jamming: `phi_J` at arrest is within 0.005 of the
   local `phi`, at every packing tested.
6. Salt moves the onset, not the arrested state.

**Not settled, and one thing withdrawn:**
- **`phi_m` is still unmeasured.** The inversion cannot supply it (§9), the
  flow-curve fit that supplied 0.5176 failed its own gate, and the calibration
  work in `tier0_wyart_cates_closure` shows a flow curve returns it only as a
  band 6x to 24x the tolerance `phi` is set to. It needs a measurement that is
  neither — compression or confining pressure. **No numerical prediction for it
  is offered here; the earlier one was withdrawn above.**
- The residual velocity dependence. Whether anything real is left once the
  pinning is accounted for is unresolved, and testing it needs a spreading law
  with an explicit `f(t)` rather than the steady one used here.
- Refitting the flow curves with `phi_m` pinned, to see whether `sigma_star`
  survives, is a session's desk work and has not been done.

## 12. A note on how this file was built

Three claims in this project were overstated and then withdrawn within the same
session: that no candidate carried the salt axis (§10 — there is barely a salt
axis in this observable), that the pore Deborah number carried both axes (§5 —
its permeability is purely geometric), and the `phi_m` = 0.441 headline above.
Each was caught by a check that could have been run first. They are recorded
rather than quietly edited out, because the pattern — a quantity that looks
measured and is actually pinned by the construction that produced it — is the
same failure the audit in §2 was designed to catch in *other people's* criteria.

---

## 13. The calibration behind every number here is withdrawn

Item 1 of the remaining desk work was "refit the flow curves with `phi_m`
pinned, and see whether `sigma_star` = 44.2 Pa survives". It cannot be done, and
finding out why is the more important result.

`lab_silica()` ships the `si_fit` preset from
`Github/Drop Impact/closure/params_fitted.yaml`, fitted 2026-07-20 to
`E:\rheology`. Since then:

- **the eleven flow curves were retracted, 2026-07-29.** `OPEN_PROBLEMS_MASTER.md`:
  *"No in-house rheometer; the eleven flow curves were retracted 2026-07-29 ...
  **No sigma\* ground truth on any axis.**"*
- **`si_fit` itself was withdrawn, 2026-08-26.** `WC_DNS_MASTER_PLAN.md`, entry
  A4 → S1: *"No cleared silica rheology exists (E:\rheology RETRACTED; 26 FEB 26
  in-house barred by the 2026-08-24 policy; nothing published on disk)."*
- the named defect: `eta_s` = 4.736 mPa.s is **4.736x water** for a water-based
  suspension, and the thickening ceiling is *exactly* proportional to `eta_s`,
  so the fit bought 4.736x more achievable thickening than water allows,
  invisibly to any low-stress measurement. Held at water, that entry reports
  `phi_0` collapsing 0.70 → 0.5306.

So `sigma_star` = 44.2 Pa is not a number with an uncertainty. It is an output
of a withdrawn fit to retracted data, and **§1's contradiction is stated between
the impact data and a parameter set nobody stands behind.**

### What survives it: a sweep instead of a fit

`run_robustness.py`. Sweep the four constants over defensible ranges — `eta_s`
1.0 to 4.736 mPa.s (water to the fit), `phi_0` 0.5306 to 0.70, `phi_m` 0.44 to
0.58, `sigma_star` 4.4 to 442 Pa — and ask which conclusions hold on every
corner. 224 valid corners of a 4⁴ grid.

| claim | holds on |
|---|---|
| **C2. The stress-activation clock cannot span the response** | **224 / 224 (100%)** |
| C1. No packing in this dataset can shear-jam | 96 / 224 (43%) |
| C3. Every arrested drop is above the closure's viscosity ceiling | 96 / 224 (43%) |

- **C2 is fully robust.** Over the entire box `f` spans 1.00x to 2.81x, never
  approaching the response's 6.35x. **The audit's central result does not depend
  on the withdrawn calibration at all** — which it could not have been shown to
  do by refitting.
- **C1 and C3 are conditional**, and on exactly one thing: they fail only where
  `phi_m` <= 0.487, i.e. where the model would jam the `phi` = 0.48 and 0.50
  drops outright and there is no finite ceiling to be above. So the §1
  contradiction is real *provided* `phi_m` > 0.487, and that is precisely the
  quantity nobody has measured.
- **C4 was over-tightened and is corrected.** "`phi_J` at arrest within 0.005 of
  `phi`" was the median at the three well-populated packings; it is 0.025 at
  `phi` = 0.40, and up to **0.057** over the parameter box. Marginal jamming
  holds qualitatively; the 0.005 figure does not.

## 14. Does the strain clock rescue the activation half?

Item 2. The steady closure sets `f = f_ss(sigma)` instantly. Richards et al.'s
Eq. 4 does not: `df/dt = -(gammadot/gamma_0)(f - f_ss)`, so from an unsheared
drop, after accumulating strain `gamma_acc`,

    f = f_ss(sigma) * (1 - exp(-gamma_acc / gamma_0)).

This is the one route that puts velocity dependence into the activation half
that the steady audit could not see, because `gamma_acc` on the Mobaseri lamella
convention grows with `U0` (3.47 at 0.75 m/s to 4.65 at 2.6 m/s). `gamma_0` is
**measured**, 0.54 to 1.64, by `dst_fit` gate G11 from Lin et al.'s reversal
transients — not fitted here.

| `gamma_0` | range of `f` | ties@5% retain | verdict |
|---|---|---|---|
| 0.54 (bottom of the measured range) | 1.20x | 0.93 | **FALSIFIED** |
| 1.64 (top) | 2.14x | 0.39 | survives |

**The strain clock rescues the activation half only at the top of the measured
`gamma_0` range, and not at the bottom.** That makes `gamma_0` a discriminating
quantity rather than a detail: the mechanism is available if `gamma_0` >~ 1.6
and unavailable if it is ~0.5, and the measurement spans both.

But see §15 before crediting the surviving row.

## 15. The pooled correlations are a Simpson's paradox

The audit's `rho` column pools all 73 drops across four packings. Restricting to
`phi` = 0.44 and arrested — 44 drops, which isolates the **velocity** axis —
**five of the six survivors reverse sign**:

| candidate | pooled | within `phi` = 0.44 | flips |
|---|---|---|---|
| `f_transient` (`gamma_0` = 1.64) | +0.289 | −0.263 | **yes** |
| accumulated strain on `L_c` | +0.286 | −0.264 | **yes** |
| Re | +0.190 | −0.261 | **yes** |
| `gammadot = U0/L_c` | +0.202 | −0.161 | **yes** |
| We | +0.015 | −0.270 | **yes** |
| **pore Deborah** | −0.517 | −0.261 | no |

The pooled sign is set by the packing axis; within a packing it reverses. The
within-cell sign is the physical one (more activation → more arrest → smaller
footprint), so **the audit's `rho` column should not be read as a mechanism
ranking.** The degeneracy test is unaffected — it compares drops that are
already tied on the group — but `rho` needs this warning attached.

Two consequences:

1. **The pore Deborah number is the only survivor whose sign is stable**, which
   promotes it from "second-strongest" to "the only one that is consistently
   pointing the same way".
2. **On the velocity axis alone, nothing is significant.** Best of the field at
   n = 44: We at rho = −0.270 (p = 0.077), strain-on-`L_c` −0.264 (0.084),
   `f_transient` −0.263 (0.085), Re and pore-De −0.261 (0.087). Five candidates
   within 0.01 of each other, none reaching p = 0.05. **The velocity dependence
   of the arrested footprint is real and this dataset cannot say what carries
   it.** The jamming deficit is constant within the cell and explains none of it.

## 16. Standing summary after the second desk pass

**Robust:**
- The stress-activation clock cannot span the response, on **every** corner of
  the withdrawn calibration's plausible parameter box.
- Accumulated strain on `U0/D0` is identically 1 and cannot produce an onset.
- Salt shifts the onset, not the arrested state.
- `beta_p` is resolved to 6.7% against a 6.35x range.

**Conditional on `phi_m` > 0.487, which is unmeasured:**
- That no packing here can shear-jam, and that every arrested drop exceeds the
  closure's viscosity ceiling.

**Withdrawn or downgraded in this pass:**
- `phi_m` = 0.441 from the inversion (§9, a tautology).
- "`phi_J` within 0.005 of `phi`" → within 0.057 over the parameter box.
- The audit's `rho` ranking, which is a pooling artefact for five of six
  survivors.

**Open, and not answerable from this dataset:**
- What carries the velocity axis. Five candidates tie at p ~ 0.08.
- `phi_m`, `sigma_star`, `eta_s`, `phi_0` — all of them. **The blocking item is
  not analysis; it is that this material has no cleared rheology at all.**
