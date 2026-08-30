# Calibrating the tier-0 closure against real flow curves

Written 2026-08-30. Closes the README's first open item, and closes it mostly
with a negative.

**Revised the same night after three independent reviews.** Two of the four
headline claims in the first draft did not survive them. The first draft is kept
as `CALIBRATION_2026-08-30.draft1.md`; §7 says what changed. Read §7 before
quoting anything from it.

Reproduce:

```
python calibration/calibrate.py          # fits, tables, out/calibration.json
python calibration/report.py             # the verdict + out/calibration_transfer.png
python calibration/synthetic_control.py  # the control that clears the fitter
python -m pytest calibration/test_calibration.py -q   # 30 passed
```

---

## 1. What was asked, and what came back

The README's open item read:

> Calibrate `(phi, phi_m, sigma_star, sigma_yield, f_res, tau)` to real flow curves.

| tier-0 parameter | can a steady flow curve set it? | notes |
|---|---|---|
| `sigma_star` | **yes** | the one clean win. Its size ratio across a 9.3x change is reproduced to 2%; it checks against an independent shear-reversal measurement at 0.73x and 1.83x (`dst_fit` gate G10, `dst_fit/scripts/gate_g10.py`) |
| `phi` | **probably not the sample's** | the inferred value runs above the packing the sample was made at, slope 1.02 to 1.25 across ten regressions — but only two of eight testable regressions exclude the identity. §2 |
| `phi_m` | **no** | however four curves of one material are chosen, its `phi_m` spans at least 0.012 and typically 0.048, against a per-curve error of 0.001 |
| `tau` | **no, structurally** | at steady state the relaxation law's left side is zero. The clock is absent from the curve, not ill-determined by it (`dst_fit` G11) |
| `f_res` | **no, structurally** | retention is a property of unloading; a steady curve never unloads |
| `sigma_yield` | not tested | only corn starch has one, and it is confounded with the thinning branch the model has no term for |

## 2. Does a flow curve report the volume fraction of the sample?

This is what matters for a solver, because a drop-impact solver advects the
volume fraction the experimentalist weighed out and needs `phi_J` on that scale.

Fit Wyart-Cates to a published flow curve with that paper's own measured
`eta_s`, `phi_0`, `n` and `beta`, and regress the inferred `phi` against the one
the authors made the sample at. Two routes: **free fit** (fit `phi, phi_m,
sigma_star`) and **plateau**, which is fit-free —
`phi = phi_0 (1 - eta_r^(-1/n))` read off the lowest-stress viscosity, no
optimiser, no starting guess, no branch cut.

| dataset | route | all usable curves | steady subset |
|---|---|---|---|
| Guy PMMA 3770 nm | free fit | 1.240 [1.093, 1.518] **excl.** | 1.115 [0.894, 1.565] |
| Guy PMMA 3770 nm | plateau | 1.020 [0.976, 1.265] | 0.990 [0.947, 1.317] |
| Richards corn starch | free fit | 1.097 [0.989, 1.204] | 1.246 [1.037, 1.384] **excl.** |
| Richards corn starch | plateau | 1.124 [1.040, 1.181] **excl.** | 1.044 [0.887, 1.203] |
| Guy PMMA 404 nm | free fit | 1.133 [0.737, 1.260] | 1.174 [0.681, 1.818] |
| Guy PMMA 404 nm | plateau | **not defined**, §3 | **not defined**, §3 |

Percentile bootstraps resampling *curves*, the independent unit. "excl." means
neither the slope's interval contains 1 nor the intercept's contains 0.

**Supported, not established.** All ten slope estimates are at or above 0.99 and
eight exceed 1.04, so the sign is consistent across two labs, two materials and
two routes. But only **two of eight testable regressions** exclude the identity,
and they are different datasets by different routes. With three materials a sign
test alone reaches only *p* = 0.125.

Two weaknesses a reader must carry:

- **Richards' corn starch is in weight fractions**, as its authors use, because
  the granules are porous. Four of the ten regressions are therefore
  weight-fraction statements, and weight-to-volume is not itself affine. The
  cleanest fit-free exclusion in the table is one of them.
- **The steady subsets hold four curves each** and their bootstrap intervals
  roughly triple; a third of the n = 4 resamples have only two distinct
  abscissae and are effectively two-point lines.

**The gauge is handled but it does not defend the headline.** Wyart-Cates is
exactly invariant under `(phi, phi_m, phi_0) -> c(...)`, so a `phi_0` that is
`c` times too large scales the inferred `phi`, the slope and the intercept
together. What survives is their ratio:

    phi_x = -intercept / slope = 0.032 (Guy 3770, plateau)
                               = 0.042 [0.013, 0.063] (corn starch, plateau)

`phi_x` is the packing at which the inference would agree with the truth. **The
slope itself is gauge-dependent**: it reaches 1 at a `phi_0` about 13% below
Guy's and 11% below Richards'. Guy's `phi_0` = 0.6123(10) is measured to 0.16%
and survives that comfortably; Richards' 0.546 is a `phiRCP` taken from a
plotting script, in weight units, with no quoted uncertainty, and does not
obviously. **Quote `phi_x`; treat the slope as gauge-bound.**

### It is not the fitter

Curves generated from the model itself, on the real stress grids, at the real
volume fractions, fitted back with the same constants:

| noise (sd of log eta) | slope | intercept |
|---|---|---|
| 0.00 | 1.0000 | +0.0000 |
| 0.03 | 1.0093 | −0.0045 |
| 0.10 | 1.1494 | −0.0804 |

The 10% row is the honest caveat: at that scatter, with the branch cut on, noise
alone manufactures 62% of the real slope excess. With the cut off the control is
clean at every level, which locates the mechanism in the cut. Rerun at each real
curve's *own* measured residual scatter (0.016 to 0.324 over six curves), 20
independent draws, at two truths, the control returns slope 0.999 and the real
1.240 lies outside its spread.

## 3. Where the fit-free route is not defined, and why that matters

The plateau inversion is the `f -> 0` limit of the model, which is the
**hard-sphere** viscosity. Wyart-Cates has no Brownian physics, so it is
meaningless where the material is Brownian at the lowest stress sampled.

The definition-free evidence, and this is the one to quote: at `phi` = 0.5363
Guy's 404 nm sweep has a lowest-stress relative viscosity of **28 904**, against
the **111.5** the frictionless branch predicts at that packing — a factor of
260. Inverting anyway returns `phi` = 0.602 for a sample made at 0.536. Pinning
`phi` there is worse: the fit then returns `sigma_star` of 16 kPa and 33 kPa
against a published 200 Pa.

The single-particle Peclet number at that point is 2×10⁻⁴, against 2 to 157 at
3770 nm — a four-decade separation, which is why the threshold is left at the
definitional `Pe = 1` rather than tuned. **That separation is
definition-dependent.** Multiplied by the crowding factor `eta_r`, as the
suspension Peclet number requires, the 404 nm series runs 0.047 to 59 and four
of its five curves would be admitted. So Pe is reported as a *diagnostic* and
the 260x viscosity discrepancy is the gate. Guy et al. separately report wall
slip at all stresses below `phi_m`, which would produce the same failure and
does not scale with Pe at all; this document cannot separate the two.

## 4. Is the half a flow curve DOES measure in play during an impact?

**The answer depends entirely on the stress gauge, and the first draft used the
wrong one.**

`Sigma = rho k U0^2` is a *total* stress. Under undrained loading — which this
silica is by five to six decades — the pore fluid carries the isotropic part and
only the deviatoric part generates particle pressure. The lab's own
`impact_stress/STRESS_AND_DILATANCY_BRIEF.md` establishes this for exactly this
geometry and states the consequence: `rho U0^2` runs 40 to 800x above
`sigma_star` across the whole accessible velocity range and so has no reachable
threshold, while the area-averaged Wagner wall shear is 5 to 10x smaller and
*does* cross `sigma_star` near the bottom of the window.

At `U0` = 1.15 m/s, `rho` = 1440 kg/m³: inertial 3.8–22.9 kPa, deviatoric
0.19–0.38 kPa.

| material | `sigma*` | inertial gauge | deviatoric gauge |
|---|---|---|---|
| Guy 3770 nm | 9.9 Pa | f = 0.997–1.000, **saturated** | f = 0.949–0.974, **not saturated** (29% of the band) |
| Richards corn starch | 5.9 Pa | f = 0.998–1.000, **saturated** | f = 0.963–0.980, **not saturated** (30%) |
| Guy 404 nm | 613 Pa | f = 0.851–0.974, not saturated | f = 0.040–0.200, barely activates |

Saturation is judged against the calibration's own resolution: what the residual
`1 - f` buys is a shift `(1-f)(phi_0 - phi_m)` in the jamming point, and that
matters only if it is comparable to the `phi_m` band. On each material's own
band, the boundary `sigma_star` above which the activation half stops being
saturated is 36–81 Pa on the inertial gauge and **1.8–4.0 Pa on the deviatoric
gauge** — below every `sigma_star` fitted here.

> **On the gauge this programme's own analysis says is correct, both halves of
> Wyart-Cates are in play during a drop impact.** The first draft concluded the
> opposite, from the total-stress gauge, and was wrong. The choice of gauge is
> the single most consequential unstated convention in this transfer.

## 5. How wide is `phi_m`, and does it matter to anyone?

| material | subset | `phi_m` band | width | random-4 null (mean / min / max) |
|---|---|---|---|---|
| Guy 3770 nm | steady, n=4 | [0.5297, 0.5429] | 0.0132 | 0.0484 / **0.0132** / 0.0579 |
| corn starch | steady, n=4 | [0.4475, 0.4593] | 0.0118 | 0.0473 / **0.0118** / 0.0595 |
| Guy 404 nm | steady, n=4 | [0.5537, 0.5631] | 0.0094 | 0.0082 / 0.0035 / **0.0094** |

**"Restricting to the steady curves narrows the band fivefold" is not evidence
about steadiness.** `phi_m` rises with `phi`, and the steady cut removes the
highest-`phi` curves, so the steady block is the tightest four-subset available
by arithmetic — for two of the three it is *exactly* the minimum over all 15 and
all 35 subsets of that size, and for the third it is the maximum.

What is defensible is the null's floor: **however four curves of one material
are chosen, its `phi_m` band is at least about 0.012 and typically 0.048, i.e.
6x to 24x the ±0.002 tolerance a suspension's volume fraction is set to.**

Does it change a verdict? Only for someone working inside a band, and
**`phi` = 0.44, this programme's own packing, is below every one of them.** The
band bites between 0.53 and 0.59 for the PMMA and 0.45 and 0.51 for the corn
starch. Evaluating the model at a band's own midpoint, as the first draft did,
makes the jamming verdict flip by construction and tests nothing.

## 6. What to actually do

1. **Take `sigma_star`,** and state the definition: the fitted value is where
   `f = 1/e`, 3.4x to 5.0x above the value a paper reads off where thickening
   visibly begins. Ratios between systems are reliable; the absolute is
   definition-bound.
2. **State the stress gauge before deciding whether `sigma_star` matters.** On
   the deviatoric gauge it does.
3. **Do not take `phi_m` as a number.** Take the band, and check whether the
   experiment's `phi` sits inside it.
4. **Do not take `phi` from the fit** without checking the low-stress viscosity
   against the frictionless branch first; a factor of 260 there is not a fitting
   problem.
5. **Do not fit `tau`, `gamma_0` or `f_res` to a flow curve.** `gamma_0` comes
   from reversal transients, measured at 0.54 to 1.64 against the 0.2 the solver
   borrows (`dst_fit` G11).
6. **What would fix `phi_m`:** a measurement that is not a steady flow curve —
   the oscillation branch Richards et al. reach, a confining-pressure jamming
   point, or compression. A bench item, not an analysis item.

## 7. What the reviews changed

Three reviewers ran on the first draft: a units audit, a science-code review and
an adversarial physics referee. Five findings changed a conclusion.

1. **The fit-free plateau route was being selected by the optimiser.** Its
   sample came from a filter on the *free fit's* bound flag, which dropped Guy's
   19% curve. With protocol C given its own admission rule, that series' plateau
   slope falls from 1.152 to **1.020** and no longer excludes the identity. The
   claim "the fit-free route excludes the identity in both places it can be run"
   is withdrawn; it excludes in one.
2. **The saturation conclusion reverses on the correct stress gauge.** §4.
3. **The "verdict flips, three for three" table was a tautology** — evaluated at
   each band's own midpoint. Replaced with whether any real packing lands in a
   band. None does, for this programme.
4. **"The steady restriction narrows the band fivefold" is arithmetic**, not
   steadiness. §5, with the null.
5. **`relative_viscosity` was hardcoded at n = 2** while Guy's `phi_m` was
   fitted at n = 2.2, understating both PMMA `eta_r` values by 2.4x. The
   exponent is now a parameter.

Fixed without changing a conclusion: three curves the fitter refused were
counted nowhere (Guy 404's 40% sheet among them, the lowest-`phi` curve of a
five-point regression); the matched synthetic control was regressed on a
different curve set than the real map it was compared against.

## 8. Still open

- **The missing low-`phi` term.** Nothing here tests whether the affine map is
  an Einstein/Batchelor correction the single-pole divergence lacks at small
  `phi`. The referee could not rule it in or out and neither can this document.
  It is the most likely benign explanation and it has not been excluded.
- The `k` = 2 to 12 added-mass coefficient comes from an object penetrating a
  suspension bath, not a drop spreading on a substrate. Its transfer across
  geometries is assumed.
- No golden-output regression on `out/calibration.json`; the pipeline is
  deterministic and should have one.
- The bootstrap itself is untested. Every interval test in the suite feeds exact
  noiseless affine data, which tests `np.polyfit` rather than the coverage.
