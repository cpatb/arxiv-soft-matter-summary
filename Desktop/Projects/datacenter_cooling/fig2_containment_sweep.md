# fig2_containment_sweep.png

Generic hall: 20 racks, 200 kW IT, `seed=1`, ASHRAE inlet limit 27 °C,
`recirc_model="cfd"` (φ(η) = 0.007 + 0.228·e^(−η/0.79), calibrated in `cfd/`).
Regenerate with `python optimize_cooling.py`.

## (a) PUE vs containment effectiveness

Two curves. The gap between them is the point.

**Red, "fixed airflow"**: provisioning held at the design point *p* = 1 and the supply
temperature raised until the hottest rack sits exactly on the 27 °C limit
(`containment_roi()`). This is the containment lever as an operator buys it:

| η | 0.0 | 0.3 | 0.5 | 1.0 |
|---|---|---|---|---|
| PUE | 1.2939 | 1.2875 | 1.2849 | 1.2813 |
| supply T (°C) | 23.01 | 24.47 | 25.09 | 26.00 |

ΔPUE across the full sweep is **0.0126**.

**Grey, "free airflow"**: PUE minimized over both (T_s, *p*) at each η
(`containment_sweep()`). It is flat at 1.2355 to 13 decimal places. That is not a plotting
error and not a solver failure. At the free optimum *p* ≈ 0.684, so the mass-balance floor
on recirculation, φ ≥ 1 − *p* = 0.316, exceeds the geometric term φ(η)·(1 + under) for
**every** η in [0, 1] (it ranges 0.094 to 0.309). Once the floor binds, η drops out of the
influence matrix entirely and containment buys literally nothing. The flat line is
therefore a statement about the airflow axis, not about containment, which is why the red
curve is the one the README quotes.

## (b) Supply temperature vs containment effectiveness

The mechanism behind (a): better containment lowers the inlet rise, so the supply air can
be warmer before the hottest rack hits 27 °C, and chiller COP rises with supply
temperature. 23.0 → 26.0 °C across the sweep.

## Axes

Panel (a) is anchored at **PUE = 1.0**, the physical floor (zero cooling overhead), rather
than at the data range. The containment lever really is worth about 1% of PUE, and a
cut axis would draw it as though it filled the panel. Panel (b) is anchored at 0 °C. Both
x-axes are anchored at η = 0. No axis is log-scaled.

## Caveats that belong with this figure

- The *p* < 1 branch of the model is **provisional**, see `DESIGN.md` §8. Over-provisioning
  carries no thermal benefit in this model (`under = max(0, 1−p)`), so the airflow axis is
  one-sided, and server-fan power is not charged, so the model can improve PUE by moving
  work into the IT denominator. Both push the free optimum toward less airflow.
- The mass-balance floor φ ≥ 1 − *p* was added 2026-08-29. Before it, the free optimum sat
  on the `p_min` = 0.6 bound with a recirculated fraction of 0.23 where conservation
  demands 0.40, and the free-airflow sweep ran 1.2270 → 1.2046 (ΔPUE 0.0224). That curve
  is superseded.
- The off-diagonal rack-to-rack coupling in `B` is still uncalibrated (2D CFD gives φ(η)
  only), so absolute hotspot temperatures, and hence the exact supply temperatures in
  panel (b), are indicative rather than predictive.
