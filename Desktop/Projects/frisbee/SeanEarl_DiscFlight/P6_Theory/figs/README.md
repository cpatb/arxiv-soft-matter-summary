# The theory figures

`figs/` is built from `out/sweep_speed.json` (the wind-speed axis, 54 cells);
`figs_flight/` is the same three figures built from `out/sweep_flight.json` (the
flight-variation axis, 30 cells). Regenerate either with

```bash
python analyse_theory.py --sweep out/sweep_speed.json
python analyse_theory.py --sweep out/sweep_flight.json --figdir figs_flight
```

Nothing is written on the axes that is not a variable name or a unit; what each
figure means is here.

## `fig_wind_speed.png` — the three contributions to landing error

Three panels on **one shared, zero-based y axis**, which is the whole point: the
comparison the figure exists to make is between the panels, not within them.
Left is the thrower's own release scatter, flat by construction because it is
measured with the wind switched off. Middle is the mean-wind bias — the part a
thrower who knows the field can aim off. Right is the gust dispersion, which
they cannot.

Colour is the wind direction it blows *from*; line style is the throw type.

The reading: **the flat line in the left panel sits above almost everything in
the other two.** The gust dispersion only overtakes the release scatter around
`Uref` = 15–20 m/s. Both wind panels are close to straight lines through the
origin, which is the theory's prediction that the response is linear in wind
speed; the curvature that is visible at the top end is the disc's own
aerodynamic nonlinearity, not the site's.

Sharing the y axis was a deliberate late change. With three autoscaled panels a
0.2 m gust effect is drawn the same height as a 7.8 m release scatter, and the
figure then says the opposite of what the numbers say.

## `fig_law_test.png` — measured dispersion against both laws

Log axes on both sides, because the quantity spans more than a decade across the
sweep; the dashed line is 1:1. Red is the brief's ballistic `Φ(N)` law, blue is
the Taylor tracer law it displaced. Both are evaluated with `σ_u` and `L`
measured along the flights themselves and with each flight's own `T_f`, `V_d`
and `L_f`. **Nothing is fitted.**

The reading, on `figs/` (wind-speed axis): the two clouds bracket the 1:1 line —
ballistic high by ~3.3×, tracer low by ~0.43× — and both are *parallel* to it.
On that axis they are not distinguishable, and the analysis says so: the two
predictions separate by only 2.13× there while the fit residual is 1.39×.

The reading on `figs_flight/` (flight-variation axis) is the one that matters.
The flight now spans `L_f` 13–71 m and `V_d` 5.5–14.4 m/s, the two predictions
separate by 3.96×, and the red cloud is visibly **shallower than 1:1** — slope
0.671, CI [0.522, 0.818]. The ballistic law's shape is rejected. The blue cloud
stays parallel — slope 0.950, CI [0.640, 1.218] — so the tracer law is not.

## `fig_wind_fraction.png` — the share of the error the wind is responsible for

Percentage of the total mean-square landing error that would vanish in still
air, against wind speed. Zero-based, capped at 100.

The reading: at 5 m/s it is 4–19 %. It reaches half only between 13 and 20 m/s
depending on direction and throw. Every value on this plot scales with the
assumed release scatter, which is **not calibrated against a measured player**;
`crossover_speed()` in `analyse_theory.py` gives the same information in a form
that carries its own scaling (`U* = σ_release / c`) so a reader can substitute
their own thrower.
