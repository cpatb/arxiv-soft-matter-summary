# Uroflow splashback: jet breakup, stream tracking, and a splash threshold

A self-contained analytical model of a bodily-fluid directed jet, and a sibling of the
bidet nozzle-optimisation paper. Both aim a liquid jet at a surface and control what
happens on impact; here the objective is inverted, and the goal is to minimise splashback.

**v2, 2026-08-07.** The v1 model (2026-07-08) is preserved in `_archive/` together with
its output. Five defects were found and fixed, and two of the three original
recommendations changed as a result. See `PHYSICS_NOTES.md` for the derivations and the
full defect list.

## The physics chain

1. **Exit.** Peak uroflow `Q ~ 15-25 mL/s` through an effective meatus diameter
   `d ~ 3-6 mm` gives `v0 = Q/(pi d^2/4) ~ 1-4 m/s`, `We ~ 100-600`, `Oh ~ 0.002`.

2. **The falling jet accelerates and thins.** Energy gives `v(s)`; continuity gives
   `r(s) = r0 sqrt(v0/v(s))`. Over a 0.75 m fall the diameter drops to 63% of its exit
   value, which matters twice: the local Rayleigh growth rate goes as `r^-3/2`, and the
   impact Weber number goes as `d`.

3. **Breakup by accumulated disturbance.** The local fastest-growing varicose mode is
   found in closed form from Weber's viscous long-wave dispersion relation. The jet
   fragments where `G(s) = int omega/v ds` reaches `Lambda = ln(r0/delta0)`.

4. **Trajectory.** Ballistic for the coherent column, with a floor check. After breakup
   the fragments fly with Schiller-Naumann drag, so the drag correction is computed
   rather than assumed (it is a fraction of a per cent here).

5. **Impact against published thresholds.** The target carries an explicit unit normal.
   Velocity is resolved into normal and tangential parts, and the splash criterion is
   evaluated on the normal part: Mundo et al. (1995) for a dry solid, Cossali et al.
   (1997) for a wetted one. The output is the **splash excess** `K/K_c - 1`, which is
   negative when no splash is predicted and comparable across geometries.

## What's here

| file | purpose |
|------|---------|
| `uroflow.py` | the physics: exit, falling jet, stability, trajectory, splash criteria |
| `test_uroflow.py` | 44 checks against analytic limits, published constants and conservation laws |
| `run_uroflow.py` | six figures and a quantitative report |
| `PHYSICS_NOTES.md` | derivations, the defect list, and what changed in the conclusions |
| `figures/` | output |
| `_archive/` | v1 model, v1 figures, v1 output |

Reproduce: `python -m pytest test_uroflow.py -q` then `python run_uroflow.py`.

## Results

**Breakup is set by the exit disturbance, not by ideal theory.** The parallel-jet
textbook result gives 63 cm against an observed 15-20 cm. Calibrating `Lambda` against
the observation gives `Lambda = 2.7-3.6`, implying an exit radius perturbation
`delta0/r0` of 3-7%. That is a falsifiable side-prediction: a high-speed clip of the
first 2 cm of the stream should show a ripple of a few per cent, not a smooth column.

**Wall and pool need opposite advice, and this is the cleanest result in the model.**
At a vertical wall the normal is horizontal, so the splash-driving velocity is
`v0 cos(theta)`, which contains no launch height at all. Aim is the entire lever, and
the three height curves in `fig6` coincide exactly wherever both are reachable. At a
horizontal pool the normal is vertical, so the splash-driving velocity is the fall
speed, height is the entire lever, and aim barely enters. **Conflating the two is
exactly the error v1 made**: it reasoned verbally about a pool while integrating to a
wall, and so recommended the head-on impact as the driest configuration.

**The practical number.** Standing 30 cm from a wetted vertical wall, the impact stays
below the Cossali threshold only for aim angles steeper than about 37 degrees below
horizontal. At 45 cm no reachable aim clears it. From a sitting height a steep aim lands
short, so the safe window closes at about 14 cm.

**Splash switches on discontinuously at the breakup point.** The right panel of `fig6`
shows the excess jumping from below threshold to well above it as the fall height passes
the 15-20 cm coherent length. The Splash Lab reported that a coherent stream splashes far
less than a drop train; here that is derived rather than assumed, and it is sharp rather
than gradual.

**A wetted surface splashes less.** Cossali's threshold rises with film depth, so over
the correlated range `0.1 <= h/D <= 1` a wetter target is monotonically harder to splash.
Drying the target is counterproductive. This lever was absent from v1 entirely.

## Honest limitations

- **Cossali's correlation is fitted for `0.1 <= delta <= 1`.** A deep toilet bowl is
  `delta >> 1`, where crown and cavity splash replace film splash. The code flags the
  extrapolation and does not correct it.
- **The meatus is a slit, not a circle.** The effective circular diameter absorbs stream
  splitting and asymmetry, which is the largest unmodelled geometric effect.
- **`Lambda` is calibrated against one reported observation range**, so the breakup point
  inherits that anchor's uncertainty. One measured breakup distance on a real stream
  would tighten it directly.
- **Gravitational stretching stabilises a falling jet** (Frankel & Weihs 1985) and is not
  modelled, so the *ideal* breakup length is a lower bound. The calibrated `Lambda`
  absorbs the difference, which means `Lambda` should be read as an effective disturbance
  level rather than as a literal `ln(r0/delta0)`.
- **Rayleigh breakup gives one drop size.** Real streams produce a distribution, and
  satellite drops are smaller and therefore harder to splash.

## Where to go next

1. **Calibrate against one phone-camera clip.** A single measured breakup distance and one
   measured impact angle would pin `Lambda` and check the predicted threshold crossing.
   This is an afternoon and it converts the model from anchored-to-literature to
   anchored-to-measurement.
2. **A curved target.** A real urinal wall is neither a plane nor horizontal; the local
   normal rotates along the surface, and the model already has the machinery to sweep it.
3. **Deep-pool splash.** Replace the extrapolated Cossali branch with a cavity/crown
   criterion so the bowl case stops being an extrapolation.
4. **Two-stream asymmetry**, the well-known stream splitting, which currently sits inside
   the effective diameter.

## References

- Rayleigh (1878); Weber (1931). Jet breakup and the viscous dispersion relation.
- Eggers & Villermaux, *Rep. Prog. Phys.* **71**, 036601 (2008). Breakup review.
- Frankel & Weihs (1985). Stability of a gravitationally stretched jet.
- Mundo, Sommerfeld & Tropea, *Int. J. Multiphase Flow* **21**, 151 (1995). Dry-surface
  splash threshold, `K = We^(1/2) Re^(1/4) ~ 57.7`.
- Cossali, Coghe & Marengo, *Exp. Fluids* **22**, 463 (1997). Wetted-surface threshold,
  `K_c = 2100 + 5880 delta^1.44`.
- Watson, *J. Fluid Mech.* **20**, 481 (1964). Radial wall jet; the transition radius is
  reused from the validated bidet model.
- Hurd, Truscott et al., "Urinal Dynamics", APS-DFD 2013 (Splash Lab). The coherent-length
  observation used as the calibration anchor.
