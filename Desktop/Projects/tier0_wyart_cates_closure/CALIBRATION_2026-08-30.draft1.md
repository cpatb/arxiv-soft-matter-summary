# Calibrating the tier-0 closure against real flow curves

Written 2026-08-30. Closes the README's first open item, and closes it with a
negative for most of it.

Reproduce:

```
python calibration/calibrate.py          # fits, tables, out/calibration.json
python calibration/report.py             # the verdict + out/calibration_transfer.png
python calibration/synthetic_control.py  # the control that says the fitter is not at fault
python -m pytest calibration/test_calibration.py -q   # 30 passed
```

---

## 1. What was asked, and what came back

The README's open item read:

> Calibrate `(phi, phi_m, sigma_star, sigma_yield, f_res, tau)` to real flow curves.

Three of those six are not in a steady flow curve at all, and the other three do
not all survive the trip into a solver. The result, in one table:

| tier-0 parameter | can a steady flow curve set it? | notes |
|---|---|---|
| `sigma_star` | **yes** | the one clean win. Ratio across a 9.3x size change reproduced to 2%; checks against an independent shear-reversal measurement at 0.73x and 1.83x |
| `phi` | **no, not the sample's** | the fit returns `1.10` to `1.24` times the packing the sample was made at, minus `0.04` to `0.14` |
| `phi_m` | **no** | one material's own curves return a band 0.023 to 0.060 wide against a per-curve error of 0.001 |
| `tau` | **no, structurally** | at steady state the relaxation law's left side is zero. The clock is absent from the curve, not ill-determined by it |
| `f_res` | **no, structurally** | retention is a property of unloading; a steady curve never unloads |
| `sigma_yield` | not tested here | only corn starch has one, and it is confounded with the thinning branch the model has no term for |

The `tau` and `f_res` rows are not new; `dst_fit`'s gate G11 established the
first and the second follows from the same argument. They are recorded here
because a calibration document that omits them invites someone to fit them
anyway.

## 2. The finding: a flow curve does not report the volume fraction of the sample

This is the part that matters for a solver, because a drop-impact solver advects
the volume fraction the experimentalist weighed out, and it needs `phi_J` on the
same scale.

Fit Wyart-Cates to a published flow curve, with that paper's own measured
`eta_s`, `phi_0`, `n` and `beta`, and regress the inferred `phi` against the
`phi` the authors made the sample at. The map is affine and its slope exceeds 1
in **every one of the ten regressions run**, across two labs, two materials, two
inference routes and two admission rules:

| dataset | route | all usable curves | steady subset only |
|---|---|---|---|
| Guy PMMA 3770 nm | free fit | 1.240 [1.093, 1.518] **excl.** | 1.115 [0.894, 1.565] |
| Guy PMMA 3770 nm | plateau, fit-free | 1.152 [1.060, 1.343] **excl.** | 1.124 [0.968, 1.793] **excl.**¹ |
| Richards corn starch | free fit | 1.097 [0.989, 1.204] | 1.246 [1.037, 1.384] **excl.** |
| Richards corn starch | plateau, fit-free | 1.124 [1.040, 1.181] **excl.** | 1.044 [0.887, 1.203] |
| Guy PMMA 404 nm | free fit | 1.133 [0.737, 1.260] | 1.174 [0.681, 1.818] |
| Guy PMMA 404 nm | plateau, fit-free | **not defined**, §3 | **not defined**, §3 |

Intervals are percentile bootstraps resampling *curves*, the independent unit.
"excl." means neither the slope's interval contains 1 nor the intercept's
contains 0. ¹ the steady-subset plateau row excludes on the intercept alone.

**The honest strength of this.** On all usable curves the fit-free route
excludes the identity in both places it can be run, which is the strongest
single statement here: it involves no optimiser, no starting guess and no branch
cut, so it cannot be an artefact of any of them. Restricting to the four curves
each source lab calls steady costs most of the statistical power — the intervals
roughly triple — and only one of four testable rows still excludes. The point
estimates do not move much (1.04 to 1.25 throughout), but with three materials a
sign test alone would only reach *p* = 0.125. **So: established on two datasets
over their full curve sets by a fit-free route, consistent everywhere, and not
established from the steady subsets alone.**

**The gauge is handled.** Wyart-Cates is exactly invariant under
`(phi, phi_m, phi_0) -> c (phi, phi_m, phi_0)`, so assuming a `phi_0` that is
`c` times too large scales every inferred `phi` by `c`, and with it both the
slope and the intercept. No choice of `phi_0` can turn this map into the
identity, because their ratio does not move:

    phi_x = -intercept / slope = 0.088 [0.057, 0.153] (Guy 3770, plateau)
                               = 0.042 [0.013, 0.063] (corn starch, plateau)

`phi_x` is the packing at which the inference would agree with the truth. Above
it the fit under-reads the volume fraction; the error shrinks toward zero as
`phi` approaches jamming, which is exactly the regime where the curve carries
information about the distance to jamming and nothing else.

### It is not the fitter

Generate curves from the model itself on the real stress grids, at the real
volume fractions, with the same constants, and fit them back:

| noise (sd of log eta) | slope | intercept |
|---|---|---|
| 0.00 | 1.0000 | +0.0000 |
| 0.01 | 1.0004 | −0.0003 |
| 0.03 | 1.0093 | −0.0045 |
| 0.10 | 1.1494 | −0.0804 |

The 10% row is the honest caveat: **at 10% scatter with the branch cut on, noise
alone manufactures 62% of the real slope excess.** The mechanism is real, so the
control was rerun at each real curve's *own* measured residual scatter (0.016 to
0.324 across the six curves), over 20 independent noise draws:

| truth | curves generatable | slope, median [95% of draws] | reproduces |
|---|---|---|---|
| `phi_m` = 0.558, published | 5 of 7 | 0.9992 [0.9972, 1.0007] | −0% |
| `phi_m` = 0.600 | 7 of 7 | 0.9985 [0.7223, 1.0168] | −1% |

Run at two truths because a synthetic curve only exists for `phi < phi_m`: at the
published value the two highest-`phi` curves cannot be generated, and those are
the two noisiest. At both, the real (1.2401, −0.1443) lies outside the control's
spread. With the branch cut off the control is clean at every noise level, which
locates the 10% mechanism in the cut rather than in the optimiser.

## 3. The admission rule this needed, and why it is not a convenience

Reading `phi` off the low-stress viscosity is the only fit-free handle a flow
curve offers. It is the `f -> 0` limit of the model, `eta_r = (1 - phi/phi_0)^-n`,
inverted. That limit is the **hard-sphere, high-Peclet** viscosity. Wyart-Cates
contains no Brownian physics, so the inversion is meaningless where the material
is Brownian at the lowest stress sampled.

Single-particle Peclet number at the lowest sampled stress, `Pe = 6 pi eta_s a^3
gamma_dot / (k_B T)`:

| dataset | Pe range | plateau admitted |
|---|---|---|
| Guy PMMA 3770 nm | 2.1 to 157 | 6 of 6 |
| Richards corn starch | 21 to 647 | 7 of 7 |
| **Guy PMMA 404 nm** | **2.3e-4 to 3.2e-3** | **0 of 5** |

Same chemistry, same lab, same figure, 9.3x apart in size, and four decades
apart in Pe. Because the separation is four decades, no threshold anywhere
between 1e-3 and 1e2 changes the verdict; the threshold is left at the
definitional `Pe = 1` rather than tuned.

The size of the error this rule prevents: at `phi` = 0.5363 the 404 nm sweep's
lowest-stress relative viscosity is 28904, against the 111.5 the frictionless
branch predicts at that packing. Inverting it anyway returns `phi` = 0.602 for a
sample made at 0.536. **Pinning `phi` on that series is worse still: the fit
then returns `sigma_star` of 16 kPa and 33 kPa against a published 200 Pa.**

That is the general warning. Protocol B (pin the known `phi`, fit `phi_m`) is
the right protocol *only where the model's assumptions hold at the low-stress
end*. On the 3770 nm series it narrows the `phi_m` band from 0.058 to 0.023. On
the Brownian 404 nm series it widens it from 0.009 to 0.075 and wrecks
`sigma_star` by two orders of magnitude.

## 4. What this does to a drop impact, which is the point

The impact stress scale is `Sigma = rho k U0^2`, with `k` = 2 to 12 the measured
added mass of a jamming front (Waitukaitis & Jaeger 2012; Peters, Majumdar &
Jaeger 2016; Han et al.). At the rig's own onset velocity `U0` = 1.15 m/s and
`rho` = 1440 kg/m³ that is **3.8 to 22.9 kPa**.

Every `sigma_star` fitted here is one to three decades below that. So during an
impact the stress-activation half of the model is pinned at `f = 1`:

| material | `sigma*` | decades below impact | `f` at impact | `phi_J − phi_m` left | vs its own `phi_m` band |
|---|---|---|---|---|---|
| Guy 3770 nm | 9.9 Pa | 2.6 | 0.997–1.000 | 0.00018 | 0.3% |
| Richards corn starch | 4.8 Pa | 2.9 | 0.998–1.000 | 0.00016 | 0.3% |
| Guy 404 nm | 514 Pa | 0.9 | 0.874–0.978 | 0.00565 | 60% |

Saturation is judged against the calibration's own resolution, not against an
arbitrary closeness to 1: what the residual `1 - f` buys is a shift in the
jamming point, `(1-f)(phi_0 - phi_m)`, and that only matters if it is comparable
to the width of the `phi_m` band the calibration returns.

**So, for micron-scale particles:**

> The half of Wyart-Cates a flow curve measures well is saturated during a drop
> impact and carries no information. The half that decides the outcome —
> `sign(phi − phi_m)` — is the half a flow curve does not measure.

**And for sub-micron particles the first half is back in play.** Bisecting the
saturation criterion puts the boundary at `sigma_star ≈ 236 Pa`; Guy's own
`sigma* ~ d^-1.9` scaling puts that at `d ≈ 370 nm` for their chemistry. Guy's
404 nm PMMA sits just on the unsaturated side, which is why its row above reads
60%. **Do not assume saturation; it depends on particle size, and the drop-impact
programme's 800 nm silica is near the boundary rather than far from it.**

### The flip, and it survives the steady restriction

`phi_m` from one material's own curves, and the range of packings over which its
own calibrations disagree about whether it jams. The steady rows use the source
lab's own admission rule (`calibration/steady.py`), parsed out of Richards'
deposited plotting script and taken from Guy et al.'s caption:

| material | subset | protocol | `phi_m` band | n | width | × the 0.002 φ tolerance |
|---|---|---|---|---|---|---|
| Guy 404 nm | all | free | [0.5537, 0.5631] | 5 | 0.0094 | 5 |
| Guy 404 nm | all | pinned | [0.5319, 0.6066] | 5 | 0.0747 | 37 |
| Guy 404 nm | steady | free | [0.5537, 0.5631] | 4 | 0.0094 | 5 |
| Guy 3770 nm | all | free | [0.5297, 0.5876] | 6 | 0.0579 | 29 |
| Guy 3770 nm | all | pinned | [0.5609, 0.5839] | 6 | 0.0230 | 11 |
| **Guy 3770 nm** | **steady** | **free** | **[0.5297, 0.5429]** | 4 | **0.0132** | **7** |
| Guy 3770 nm | steady | pinned | [0.5609, 0.5791] | 4 | 0.0182 | 9 |
| corn starch | all | free | [0.4475, 0.5070] | 7 | 0.0595 | 30 |
| **corn starch** | **steady** | **free** | **[0.4475, 0.4593]** | 4 | **0.0118** | **6** |
| corn starch | steady | pinned | [0.4516, 0.4597] | 4 | 0.0082 | 4 |

**The dramatic version of this claim does not survive the restriction and the
load-bearing version does.** Cutting to the curves each lab calls steady takes
the bands from 0.058 and 0.060 down to 0.013 and 0.012, a factor of about five.
That matches what `dst_fit`'s gate N8 found and it must be quoted that way. But
0.012 is still **six times** the tolerance to which a suspension's volume
fraction is actually set, and the verdict flips right across it. Worked, on the
steady subset, all three materials:

| material | `phi` | at the low `phi_m` | at the high `phi_m` |
|---|---|---|---|
| Guy 3770 nm | 0.5363 | jammed, `eta_r` = ∞ | flowing, `eta_r` = 1.5e4 |
| Guy 404 nm | 0.5584 | jammed, `eta_r` = 8.0e4 | flowing, `eta_r` = 5569 |
| corn starch | 0.4534 | jammed, `eta_r` = ∞ | flowing, `eta_r` = 5678 |

Same material, same paper, same figure, only the choice of which of its own
steady flow curves you calibrated on. Three for three.

Two notes on this table, both found by the units audit and both fixed rather
than caveated. The two PMMA `eta_r` values are evaluated at **n = 2.2**, the
exponent Guy et al. fitted their own `phi_m` at; the tier-0 module's
`relative_viscosity` was hardcoded at n = 2, which understated them by 2.4x. A
`phi_m` fitted at one exponent cannot be evaluated at another, and the exponent
is now a parameter with no default silently in play. And **corn starch's "x
tolerance" column is not comparable with the PMMA rows**: its `phi` is a weight
fraction, its `phi_0 = 0.546` is in weight units, and the ±0.002 yardstick is a
volume-fraction tolerance from the silica rig. The corn-starch factor of 4 to 6
is a within-material statement about that band's own width, not a like-for-like
comparison with Guy's. Converting would need a granule porosity nobody has
measured, which is why `dst_fit` ships the conversion and never calls it.

## 5. What to actually do

For anyone taking numbers from a flow curve into the Basilisk closure:

1. **Take `sigma_star`.** State the definition with it — the fitted value is
   where `f = 1/e`, which runs 3.4x to 5.0x above the value a paper reads off
   where thickening visibly begins. The ratio between two systems is reliable;
   the absolute is definition-bound.
2. **Do not take `phi_m` as a number.** Take the band, and check whether the
   experiment's `phi` sits inside it. If it does, the closure cannot tell you
   whether the material jams, and no amount of solver resolution will fix that.
3. **Do not take `phi` from the fit.** Use the measured one, and check the
   Peclet number before pinning it — pinning on a Brownian series is worse than
   not pinning.
4. **Do not fit `tau`, `gamma_0` or `f_res` to a flow curve.** `gamma_0` comes
   from reversal transients and is measured at 0.54 to 1.64, against the 0.2 the
   solver currently borrows (`dst_fit` G11).
5. **What would fix `phi_m`.** It needs a measurement that is not a steady flow
   curve: the shear-jamming boundary reached through the oscillation branch
   Richards et al. access, a confining-pressure jamming point, or a direct
   compression measurement. That is a bench item, not an analysis item.

## 6. Limits of this document

- Three datasets, two materials, one of them corn starch in weight fractions.
  The transfer map is established on two of the three; the third could not test
  it.
- **The transfer map is not established on the steady subsets alone.** Every
  table is reported both ways for exactly this reason. The steady subsets hold
  four curves each, the bootstrap intervals roughly triple, and one of four
  testable rows still excludes the identity. What survives the restriction
  cleanly is the *flip*, not the map.
- **`Guy 404 nm`'s band is identical before and after the restriction**, so its
  steady row is not an independent test. Two reasons, both worth knowing. The
  rule is applied to the *true* volume fraction, and that series' sheet labelled
  "56%" is φ = 0.5562, just under the 0.56 threshold, so only the 58% curve is
  removed. And that curve's `phi_m` = 0.5596 sits in the interior of the band,
  so dropping it moves neither end.
- `Sigma = rho k U0^2` is a scale, not a solution. The claim it supports —
  one to three decades above `sigma_star` — is robust to the whole `k` range and
  to a factor of a few in `U0`, and nothing here rests on its precision.
- The corn-starch Peclet numbers use a lower-bound granule radius of 2.5 µm
  because the deposit records no size; the paper's own first sentence calls the
  system non-Brownian, and the bound clears the threshold by one to three
  decades, so the unknown does not bite.
- `sigma_yield` was not calibrated. Doing it honestly needs a thinning term the
  model does not have.
