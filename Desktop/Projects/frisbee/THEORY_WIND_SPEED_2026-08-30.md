# Throwability as a function of wind speed

**2026-08-30. Code in `SeanEarl_DiscFlight/P6_Theory`, results in its `out/`.**

Everything this project had computed sat at one wind speed: 5 m/s at 10 m. This
adds the wind-speed axis and asks what a general theory of disc flight in wind
would have to say. Sean Earl Field becomes the worked example rather than the
subject.

Three things came out of it. One is a gate that makes the whole axis affordable.
One is a law that works better than it has any right to. One is a result about
which nobody has been asking the right question.

---

## 1. One CFD solve covers every wind speed, to 0.15 %

A speed axis with nine points and three directions is 27 OpenFOAM solves, which
is days. The standard escape is Reynolds-number independence: around sharp-edged
bluff bodies at high Re the separation is fixed by the geometry, so the steady
solution is self-similar in the inflow speed and

```
U → s·U        k → s²·k        ω → s·ω        ε = Cμkω → s³·ε
```

Three consequences carry the rest of the argument:

| | scaling | |
|---|---|---|
| `σ_u = sqrt(2k/3)` | `~ s` | the disturbance amplitude scales |
| `T_L = k/ε` | `~ 1/s` | eddies turn over faster |
| `L = σ_u·T_L` | **`~ s⁰`** | **the integral length scale does not move** |
| `TI = σ_u/U` | `~ s⁰` | nor does the turbulence intensity |

So a site's turbulence *structure* is wind-speed independent and only its
amplitude changes. **Every wind-speed dependence of throwability comes from the
disc, not from the site.**

That is not free here. `nu = 1.5e-5` is finite, and the `buildings` patch uses
the smooth-wall `nutkWallFunction`, whose blending carries an explicit Reynolds
number; only the `ground` patch's roughness-based `atmNutkWallFunction` is
scale-free. So it was measured rather than assumed. The identical gridL1 mesh
was re-solved at `Uref = 2.5, 5, 10` — same geometry, same schemes, same
800-iteration cycle average, only the three inlet `Uref` entries changed, with
the case generator asserting it patched exactly three of them.

Dividing out the expected powers, over the playing surface, relative to each
field's own RMS:

| Uref | `|U|` median / max | `k` | `ε` |
|---|---|---|---|
| 2.5 m/s | 0.14 % / 0.78 % | 0.04 % / 0.43 % | 0.03 % / 0.37 % |
| 10 m/s | 0.06 % / 0.25 % | 0.03 % / 0.21 % | 0.02 % / 0.15 % |

and the derived quantities the theory actually consumes:

| Uref | `U_mean` | `σ_u` | `T_L` | `L` | `TI` |
|---|---|---|---|---|---|
| 2.50 | 0.473 | 0.118 | 27.04 | **3.196** | 0.250 |
| 5.00 | 0.946 | 0.236 | 13.51 | **3.194** | 0.250 |
| 10.00 | 1.892 | 0.473 | 6.753 | **3.192** | 0.250 |

Fitted exponents against predicted, over a 4× speed range, pre-registered
tolerance ±0.10: `σ_u` +1.0000 (pred +1), `T_L` −1.0007 (−1), `L` −0.0007 (0),
`U_mean` +0.9992 (+1), `TI` +0.0008 (0). **PASS on all five.** The 2.5 m/s case
is about twice as far off as the 10 m/s one, which is the direction the
wall-function argument predicts.

Two limits worth stating rather than burying. The inlet is a *neutral* Richards
& Hoxey ABL, so this is a family of neutral boundary layers at different speeds
and says nothing about stratification — a real 1 m/s day usually is not neutral.
And the test ran on gridL1 at one heading; whether the residual is mesh- or
direction-dependent was not measured.

---

## 2. The gust dispersion is linear in wind speed, and it is the smallest term

Each cell of the sweep is flown three ways on the **same release draws**: no
wind at all, the RANS mean with gusts suppressed, and the full field. The gust
and mean-wind contributions come out as paired per-throw displacements rather
than as differences of variances.

That choice was forced. Release scatter alone is `σ ≈ 7.4–7.8 m` on a 27 m throw
while the gust term is 1–2 m, so `σ_full² − σ_mean²` is a ~10 % difference of two
variances each carrying a 20 % relative error at n = 48. The first attempt
returned a **negative variance**. The paired displacement measures the same thing
with the shared release draw cancelling exactly instead of statistically.

### The law holds

`σ_gust ~ Uref^p`, six independent (direction, throw) combinations, nine speeds
each:

| wind from | throw | p | 95 % CI | r² |
|---|---|---|---|---|
| 90° easterly | backhand | 1.104 | [1.042, 1.211] | 0.996 |
| 90° easterly | forehand | 1.001 | [0.996, 1.012] | 1.000 |
| 210° SSW | backhand | 1.027 | [0.985, 1.089] | 0.998 |
| 210° SSW | forehand | 0.889 | [0.846, 0.908] | 0.999 |
| 270° westerly | backhand | 1.171 | [0.989, 1.423] | 0.977 |
| 270° westerly | forehand | 1.260 | [1.126, 1.469] | 0.990 |

Predicted exactly 1. Measured 0.89 to 1.26. The mean-wind bias behaves the same
way, `q = 0.85` to `1.12`. The scatter around 1 is the disc's own nonlinearity:
the wind changes the airspeed, hence the lift, hence the flight, and the
tabulated Potts & Crowther coefficients are strongly nonlinear in angle of
attack. The departure is largest at 270°, which is also where the flight is most
disturbed.

### One number per direction per throw

Because the response is linear, a site-direction-throw combination is fully
characterised by a single slope. Fitting `σ_wind = c·Uref` through the origin,
where `σ_wind` collects bias, mean-shear and gust terms:

| wind from | throw | `c` (m of error per m/s of reference wind) | `U*` | linear residual |
|---|---|---|---|---|
| 90° easterly | backhand | 0.379 | 20.5 | 3.5 % |
| 90° easterly | forehand | 0.279 | 26.5 | 5.9 % |
| 210° SSW | backhand | 0.498 | 15.6 | 4.0 % |
| 210° SSW | forehand | 0.580 | 12.7 | 11.4 % |
| 270° westerly | backhand | 0.581 | 13.4 | 12.8 % |
| 270° westerly | forehand | 0.850 | 8.7 | 26.3 % |

`c` is the wind-sensitivity coefficient of the site in that direction, and it is
the quantity a general theory should be reporting. **The westerly is 1.5–3×
more wind-sensitive than the easterly.** The residual column is a nonlinearity
meter: it stays under 6 % where the throw is least disturbed and reaches 26 % on
the 270° forehand, marking where the disc's own aerodynamics start to matter
more than the wind's amplitude.

`U*` is the wind speed at which the wind terms equal the thrower's own scatter.
It is reported alongside `c` because **it scales linearly with whatever release
scatter you assume**: `U* = σ_release / c`. The release distribution used here
is `trajectories.RELEASE_DISTRIBUTIONS`, which its own docstring calls
"practitioner-typical" and which is *not* calibrated against a measured player.
It produces a 7.4–7.8 m landing scatter on a 27 m throw, which is poor. A player
with half that scatter reaches `U*` at half the wind speed.

### And this is the result nobody was asking about

The error budget, at a 3 m catch radius:

| Uref | wind's share of the landing error (270° backhand) |
|---|---|
| 2 m/s | 1.8 % |
| 5 m/s | 9.5 % |
| 10 m/s | 27.6 % |
| 14 m/s | 51.7 % |
| 20 m/s | 72.2 % |

At the 5 m/s this entire project was computed at, **the wind accounts for 4–19 %
of the landing error and the thrower for the rest.** The whole apparatus — the
CFD, the gust model, the 97,200-throw sweep, GATE G1 — is resolving a correction
of order 10 %, and the corresponding statement was never made because the
pipeline had no still-air arm to compare against.

This does not make the wind work wrong. The direction contrast is real and it is
measured on paired throws, so it survives a large common floor. But any sentence
of the form "this site is hard to throw on" has to name the wind speed, and at
this site, in ordinary wind, the honest answer is that the thrower dominates.

---

## 3. Neither published law predicts the magnitude

The brief's own ballistic law and the tracer law it displaced, both evaluated
with path-averaged `σ_u` and `L` measured on the same flights and the flights'
own `T_f`, `V_d`, `L_f` — nothing fitted:

| law | predicted / measured, median | 10–90 % |
|---|---|---|
| ballistic `Φ(N)` | **3.27×** | [2.02, 5.03] |
| tracer `1/√N` | **0.43×** | [0.28, 0.55] |
| geometric mean of the two | 1.20× | |

The ballistic law over-predicts by about a factor 3 and the tracer law
under-predicts by about a factor 2.3. Their geometric mean happens to land at
1.2×, which is almost certainly a coincidence and is quoted only so nobody
mistakes it for a proposal.

A constant factor is the least interesting kind of wrong — it could be `κ`. What
matters is whether either law *tracks* the measurement as the flight and the
wind change. Regressing `log(measured)` on `log(predicted)`:

| law | slope b | 95 % CI | prefactor needed | residual scatter |
|---|---|---|---|---|
| ballistic | 1.075 | [0.990, 1.175] | 0.306× | 1.39× |
| tracer | 1.061 | [0.978, 1.162] | 2.352× | 1.36× |

Both have slope 1 — **and on this sweep that is not evidence for either.** The
two laws' predictions differ only through `T_f`, `V_d` and `L_f`, and across a
wind-speed sweep those move by 1.7×, 1.3× and 1.9× while `σ_u` moves by 42.6×.
The laws separate by only 2.13× across the whole set against a 1.39× residual
scatter. The analysis reports that ratio with the verdict and refuses to call
the test informative below 1.6×, because a discriminator that passes both
candidates has discriminated nothing.

Separating them needs the flight varied hard, which is what
`--release-speeds`/`--elevations` are for and what
`P1_AnalyticalModel/throw_catalogue.py` calibrates a range ladder for. That is
the brief's WP2 P2 and it is the obvious next step.

### A gap in one of the brief's formulas

§C5 states the ratio of the two laws' amplitudes as `(κ/2)·Fr⁻¹`, "i.e. ~6–8×".
Working it through with `T_f = L_f/V_d`,

```
σ_ballistic / σ_tracer = (κ/2) · Fr⁻¹ · sqrt(Φ(N)·N),      Fr⁻¹ = g L_f / V_d²
```

The `sqrt(Φ(N)·N)` factor is missing from the stated formula. At `N = 4` it is
1.334, taking `(κ/2)Fr⁻¹ = 6.08` to 8.11. The brief's own quoted numbers
(3.4/5.5/6.1 m against a tracer 0.4–0.8 m) imply ~8.5 — **the numbers are right
and the prose formula is incomplete.** Recovering the brief's inputs from two of
its other statements (`Fr⁻¹ = 2.45` at `L_f = 25` m fixes `V_d = 10.0` m/s; the
"6–8×" fixes κ as the isotropic 4.96) reproduces 3.13 / 5.31 / 5.58 m, within
4–9 % of the quoted values, the residual being that the brief band-averages
`σ_u` over 1–6 m where this uses the 2 m pitch mean.

`Φ`'s anchors were re-derived symbolically rather than quoted: `Φ(0) = 1`,
`N·Φ → 8/3`, the `V_d` exponent is exactly −3, and the brief's "1.7" for `L_f`
is the *local* log-slope at `N ≈ 4` — 1.6955 — not a universal exponent. It
drifts with `N` and is now recomputed at each cell's own `N`.

---

## 4. What is open

- **The discriminating test has not been run at enough flight variation.** The
  release-speed × elevation sweep is the instrument; the throw-length catalogue
  (`pass_15m` at 10.53 m/s → 15.0 m; `huck_50m` at 24.80 m/s → 50.0 m) is the
  other.
- **Why the ballistic law over-predicts by 3×.** The first candidate is that its
  derivation contains no spin: a disc is a gyroscope, and a gust's rolling
  moment precesses into a heading change a quarter cycle later, so gyroscopic
  stiffness should suppress exactly the response the law integrates. The sweep
  now carries a `--spins` axis to test it. This is question 3 of
  `FRISBEE_GENERALIZATION.md`.
- **The release-scatter floor is uncalibrated** and it sets every absolute number
  in the error budget. One afternoon of measured throws would fix it, and it is
  cheaper than any of the CFD.
- **Stratification.** The whole scaled family is neutral.

---

## 5. Reproducing

```bash
cd SeanEarl_DiscFlight
bash P3_CFD/runScaleTest.sh                      # ~40 min of CFD, once
python P6_Theory/wind_scaling.py                 # the gate
python P1_AnalyticalModel/calibrate_catalogue.py # the throw ladder
python P6_Theory/sweep_theory.py --out P6_Theory/out/sweep_speed.json
python P6_Theory/analyse_theory.py --sweep P6_Theory/out/sweep_speed.json
python -m pytest P6_Theory/test_theory.py -q     # 16 tests
```
