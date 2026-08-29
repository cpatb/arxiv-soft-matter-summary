# Physics notes: derivations, defects, and what changed

Written 2026-08-07 alongside the v2 rewrite. Everything asserted here is either derived
below or checked in `test_uroflow.py`.

---

## 1. The closed-form growth rate

Weber's viscous long-wave dispersion relation for varicose perturbations on a jet of
radius `r`, with `x = k r` and `nu = mu/rho`:

```
omega(x) = -1.5 nu k^2 + sqrt[ (1.5 nu k^2)^2 + (sigma / (2 rho r^3)) x^2 (1 - x^2) ]
```

Write `A = 1.5 nu / r^2`, `B = sigma / (2 rho r^3)`, `y = x^2`. Then

```
omega(y) = -A y + sqrt[ (A^2 - B) y^2 + B y ].
```

Setting `d omega / d y = 0` gives `f'(y) = 2 A sqrt(f(y))` with
`f(y) = (A^2 - B) y^2 + B y`. Squaring once and cancelling a factor of `B`:

```
4 (A^2 - B) y^2 + 4 B y - B = 0     ->     y_opt = sqrt(B) / (2 (A + sqrt(B))).
```

The inviscid limit `A -> 0` gives `y = 1/2`, so `x = 1/sqrt(2)` and
`omega_max = sqrt(sigma / (8 rho r^3)) = 0.35355 sqrt(sigma/(rho r^3))`. That is the
maximum of the long-wave relation; the exact Rayleigh (Bessel) result is
`0.3433 sqrt(sigma/(rho r^3))`, so the long-wave approximation carries a **+3.0% bias**
in the growth rate. The tests assert both the closed form and the size of that bias, so
it is a known offset rather than a hidden one.

The root was verified symbolically (sympy) and is cross-checked numerically against a
bounded optimiser across three diameters and four viscosities spanning twelve decades.

**Why it matters practically:** the closed form made the falling-jet integral about a
thousand times faster, which is what made the parameter maps in `run_uroflow.py`
feasible at all. The v1 code called a scalar optimiser once per station.

---

## 2. The falling jet

Energy along a free streamline, with air drag on the column and viscous dissipation both
negligible at `Oh ~ 2e-3`:

```
v(s)^2 = v0^2 + 2 g * fall(s)
```

Continuity, `Q = pi r^2 v = const`:

```
r(s) = r0 sqrt( v0 / v(s) )
```

For `v0 = 1.67 m/s` and a 0.75 m fall this gives `v = 4.18 m/s` and `d = 2.53 mm`
against a 4.00 mm exit, i.e. **63% of the exit diameter**. Both relations are asserted in
the tests as conservation checks rather than as formulas.

### Breakup criterion

Rather than a single `L_b` computed from exit conditions, the disturbance is accumulated
along the jet:

```
G(s) = int_0^s omega_max(r(s')) / v(s') ds' ,        breakup where G = Lambda.
```

With no vertical drop this reduces exactly to the textbook `L_b = v/omega_max * Lambda`,
which the tests check to 0.2%.

### Calibration of Lambda

`Lambda = ln(r0/delta0)` is the logarithmic amplification the initial disturbance needs
before the jet pinches. The textbook value for a quiet laboratory jet is about 12. Fitting
`Lambda` so that the model breaks at the observed 15-20 cm gives

```
Lambda = 2.70 to 3.55      ->      delta0/r0 = 2.9% to 6.7%
```

A meatus is a slit with a turbulent, wetted exit, so a few per cent initial radius ripple
is entirely plausible. This turns v1's admitted 3-4x over-prediction into a calibrated
model with a stated anchor and a testable side-prediction.

**Caveat on interpreting Lambda.** A gravitationally stretched jet is more stable than a
parallel one (Frankel & Weihs 1985), and that stabilisation is not modelled. The fitted
`Lambda` therefore absorbs both the true initial disturbance and the missing
stabilisation, so it should be read as an effective disturbance level. The predicted
breakup point is unaffected, because it is what was fitted.

---

## 3. Drop size after breakup

One wavelength of jet collapses into one drop. With `lambda = 2 pi r / x_opt`,

```
(pi/6) D^3 = pi r^2 lambda    ->    D = (12 pi / x_opt)^(1/3) r ... = 1.89 d_jet
```

for the inviscid optimum `x_opt = 1/sqrt(2)`. v1 used the **jet** diameter as the drop
diameter in the impact Weber number, under-stating `We` by a factor of 1.89 in exactly
the fragmented regime it identified as the worst case.

---

## 4. Target geometry, and the v1 inversion

This is the defect that changed the conclusions, so it is worth stating precisely.

v1 integrated the trajectory to a **vertical plane** at `x = horizontal_target`, then
scored splash with

```
angle_factor = sin( arctan2(|vz|, |vx|) )
```

For a vertical wall the outward normal is horizontal. The normal velocity component is
therefore `|vx|`, and the tangential component is `|vz|`. So
`sin(arctan2(|vz|,|vx|)) = |vz|/|v|` is the **tangential** fraction, not the normal one.
v1 multiplied its splash index by the tangential fraction, which means:

| arrival | normal component | v1 factor | v1 verdict | correct verdict |
|---|---|---|---|---|
| horizontal (head-on) | maximum | 0.00 | driest | **wettest** |
| 45 degrees | mixed | 0.71 | mid | mid |
| steep (grazing the wall) | minimum | 1.00 | wettest | **driest** |

The verbal reasoning in the v1 README ("aim shallow, grazing impact makes sin(angle)
small") is correct for a **horizontal** target, where the normal is vertical. The code
integrated to a vertical one. The two were never reconciled, and the recommendation
inherited the mismatch.

v2 makes the target an object carrying an explicit unit normal, resolves the impact
velocity against it, and evaluates every criterion on the normal component. Two
regression tests (`test_wall_normal_is_horizontal`, `test_pool_normal_is_vertical`) pin
this down.

### The consequence, which is a genuine result

At a vertical wall, `v_normal = v0 cos(theta)` and there is **no launch height in it**.
At a horizontal pool, `v_normal` is the fall speed and there is **no aim angle in it**.
So the two geometries have orthogonal levers, and advice for one is useless for the
other. That is now the headline figure.

---

## 5. Splash criteria

v1 used `coherence_penalty * sqrt(We) * sin(angle)`, a composite that is monotone in the
right variables but has no threshold: it could rank configurations and never predict
whether anything splashes.

v2 uses published criteria, evaluated on the normal component (standard practice for
oblique impact):

- **Dry solid**, Mundo, Sommerfeld & Tropea (1995):
  `K = We^(1/2) Re^(1/4)`, deposition below `K ~ 57.7`. Roughness and gas pressure move
  this by order-unity factors, so it is treated as an order-of-magnitude criterion.
- **Wetted surface**, Cossali, Coghe & Marengo (1997):
  `K = We * Oh^(-0.4) = We^0.8 Re^0.4`, threshold `K_c = 2100 + 5880 delta^1.44` with
  `delta = h_film / D`, correlated over roughly `0.1 <= delta <= 1`. Outside that range
  the code sets an `extrapolated` flag, because a silently extrapolated correlation is
  how a wrong headline gets made.

The reported quantity is the **splash excess** `K/K_c - 1`. It is zero at threshold,
negative when no splash is predicted, and comparable across geometries and film depths.

Scaling checks in the tests: `K_dry ~ v^(5/4)` and `K_film ~ v^2`, both asserted exactly.

---

## 6. Drag, computed rather than asserted

v1 asserted that drag is negligible. v2 integrates each fragment with Schiller-Naumann
sphere drag (`Cd = 24/Re (1 + 0.15 Re^0.687)` below `Re = 1000`, else 0.44) and reports
the fractional speed deficit against the drag-free flight over the same time.

For the ~6 mm Rayleigh drops this project produces, the deficit is **a few tenths of a
per cent**, so v1's assertion was correct. It is now a number in the output instead of a
claim in a docstring, and the tests confirm the deficit grows as drops shrink.

---

## 7. The floor guard

A defect introduced in the v2 rewrite and caught by its own critical-angle sweep: the
wall trajectory did not check whether the stream reaches the floor before the wall, so a
steeply-aimed stream was reported as arriving at a distant wall from below ground. Fixed
with an explicit floor crossing time and a dedicated `StreamMissesTarget` exception, so
callers can skip unreachable configurations without swallowing genuine input errors.

This is worth recording because it is the same failure class the drop-impact programme
keeps hitting: a diagnostic that returns a number for a configuration that does not
exist. The regression test is `test_steep_aim_that_lands_short_is_rejected`.

---

## 8. Summary of what changed in the conclusions

| v1 conclusion | v2 verdict |
|---|---|
| "Stand close so the stream stays coherent" | Right idea, wrong window. The coherent length is 15-20 cm, not 63 cm, so from a standing height it is unreachable. That, rather than path length in general, is why a lower launch helps. |
| "Aim shallow; grazing impact makes sin(angle) small" | **Reversed for a wall.** Shallow means head-on and is the worst case. The correct rule is to minimise the normal component, which at a wall means arriving steeply. |
| "Sit or lower the source: drier at every distance" | Only at a wall via reachability and coherence, not via the impact itself, since the wall-normal component contains no height. At a pool, height is the dominant term and lower genuinely is drier. |
| Splash index spans about 30x between best and worst | Replaced. The model now predicts a threshold crossing: below about 37 degrees of downward aim at 30 cm, no splash. |
| Standing water not modelled | Now a first-class lever, and a wetter surface splashes less over the correlated range. |
