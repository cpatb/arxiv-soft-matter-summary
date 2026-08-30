# P6_Theory — throwability as a function of wind speed

What P1–P5 built is a study of one site in three wind directions at one wind
speed. This stage turns it into a statement about wind, of which Sean Earl Field
is then the worked example.

The question is the obvious one nobody in the pipeline could previously ask:
**how does the wind speed change how well a disc can be thrown, and what else
does it depend on?**

---

## 1. The move that makes it affordable

Every OpenFOAM case in `P3_CFD/cases` is solved at `Uref = 5.0` m/s at 10 m. One
solve per direction takes hours; a speed axis with nine points would take days.

The escape is Reynolds-number independence. For flow around sharp-edged bluff
bodies at high Re the separation points are fixed by the geometry, not by the
boundary layer, and the steady RANS solution is self-similar in the inflow
speed:

```
U(x; s·Uref) = s·U(x; Uref)      k → s²·k      ω → s·ω      ε = Cμkω → s³·ε
```

Three consequences, and they are the whole backbone:

| quantity | scaling | meaning |
|---|---|---|
| `σ_u = sqrt(2k/3)` | `~ s` | the amplitude of the disturbance scales |
| `T_L = k/ε` | `~ 1/s` | eddies turn over faster in a faster wind |
| `L = σ_u·T_L` | **`~ s⁰`** | **the integral length scale does not change** |
| `TI = σ_u/U` | `~ s⁰` | nor does the turbulence intensity |

So under this scaling a site's turbulence *structure* — its length scales, its
turbulence intensity, its spatial pattern of shelter — is independent of wind
speed. Only the amplitude moves. **Every wind-speed dependence of throwability
therefore comes from the disc, not from the site.**

That is strong enough to be worth doubting, and it is not exact here: `nu` is
finite and the `buildings` patch uses the smooth-wall `nutkWallFunction`, whose
blending carries a Reynolds number. The `ground` patch's `atmNutkWallFunction`
is roughness-based and scale-free; the buildings' is not.

### The test, and its result

`P3_CFD/runScaleTest.sh` re-solves the *identical* gridL1 mesh at
`Uref = 2.5, 5, 10` m/s — same geometry, same mesh, same schemes, same 800-iteration
cycle average, only the three inlet `Uref` entries changed
(`P3_CFD/makeScaleCase.py` writes the case and asserts it patched exactly three).
`wind_scaling.py` divides out the expected powers and reports what is left, over
the playing surface, relative to each field's own RMS.

| Uref | `|U|` | `k` | `ε` |
|---|---|---|---|
| 2.5 m/s | 0.14 % median, 0.78 % max | 0.04 %, 0.43 % | 0.03 %, 0.37 % |
| 10 m/s | 0.06 % median, 0.25 % max | 0.03 %, 0.21 % | 0.02 %, 0.15 % |

and the derived quantities the theory actually uses:

| Uref | `U_mean` | `σ_u` | `T_L` | `L` | `TI` |
|---|---|---|---|---|---|
| 2.50 | 0.473 | 0.118 | 27.04 | 3.196 | 0.250 |
| 5.00 | 0.946 | 0.236 | 13.51 | 3.194 | 0.250 |
| 10.00 | 1.892 | 0.473 | 6.753 | 3.192 | 0.250 |

Fitted exponents against predicted, over a 4× speed range:
`σ_u` +1.0000 (pred +1), `T_L` −1.0007 (pred −1), `L` −0.0007 (pred 0),
`U_mean` +0.9992, `TI` +0.0008. Pre-registered tolerance ±0.10 in the exponent.
**PASS on all five.**

The 2.5 m/s case is about twice as far off as the 10 m/s case, which is the
right direction: lower Reynolds number, more room for the wall functions to
matter. So **wind speed is a free parameter of this solution to better than
0.15 % over 2.5–10 m/s**, and one solve per direction covers the axis.

Two caveats that are not in the numbers. The inlet is a *neutral* Richards &
Hoxey ABL, so the scaled family is a family of neutral boundary layers; a real
1 m/s day is usually not neutral, and this says nothing about stratification.
And the test was run on gridL1 at one heading — whether the residual is
mesh-dependent or direction-dependent is not measured.

---

## 2. The dispersion law, and what it predicts about wind speed

Not new here. `CLAUDE_CODE_FRISBEE_HARDENING.md` §C5 derives it and
`P4_Coupling/regrade_gate_g1.py` implements `Φ`; `dispersion_law.py` re-exports
that implementation rather than making a second one.

```
σ_landing = (κ · g · σ_u · T_f²) / (2 V_d) · sqrt(Φ(N)),        N = L_f / L
Φ(N) = 8[ 1/(3N) − 1/(2N²) − e^−N/N³ + (1 − e^−N)/N⁴ ]
```

The disc is a high-Stokes body — ballistic length 51 m against `L ~ 3–5 m`, so
`St ~ 8.5` — and does **not** follow eddies. It responds through force, the
displacement is a *double* time integral, and the Taylor tracer law
`σ ~ σ_u T_f/√N` does not apply. `dispersion_law.py` implements the tracer law
too, because a test that only ever evaluates the favoured law is not a test.

`check_exponents_symbolically` re-derives the anchors with sympy instead of
quoting them: `Φ(0) = 1`, `N·Φ → 8/3`, the `V_d` exponent is exactly −3, and the
brief's "1.7" for `L_f` is the **local** log-slope at `N ≈ 4` (1.6955) — it is
not a universal exponent and drifts as `N` moves, so it is recomputed at each
cell's own `N`.

Substituting the scaling of §1 at a **fixed flight**: `σ_u ~ Uref` and `N ~ Uref⁰`,
so

> **`σ_landing ~ Uref`, exactly linear.** No threshold, no exponent to fit.

That is the prediction worth attacking, because a real throw is not a fixed
flight: the mean wind changes airspeed, hence lift, hence `T_f`, `V_d` and `L_f`,
and the disc's coefficients are strongly nonlinear in angle of attack. So the
law predicts a linear **floor**, and any departure from it is disc aerodynamics
rather than site meteorology.

### A gap in one of the brief's prose formulas

§C5 says the ratio of the ballistic to the tracer amplitude is `(κ/2)·Fr⁻¹`,
"i.e. ~6–8×". Working it through, with `T_f = L_f/V_d`,

```
σ_ballistic / σ_tracer = (κ/2) · Fr⁻¹ · sqrt(Φ(N)·N),     Fr⁻¹ = g L_f / V_d²
```

The `sqrt(Φ(N)·N)` factor is missing from the stated formula. At `N = 4` it is
1.334, so `(κ/2)Fr⁻¹ = 6.08` becomes 8.11. The brief's own quoted numbers
(3.4/5.5/6.1 m ballistic against 0.4–0.8 m tracer) imply ~8.5, i.e. the numbers
are right and the prose formula is incomplete. Reproducing the brief's inputs
from two of its other statements (`Fr⁻¹ = 2.45` at `L_f = 25` m fixes
`V_d = 10.0` m/s; the "6–8×" fixes κ as the isotropic 4.96) recovers
3.13 / 5.31 / 5.58 m, within 4–9 % of the quoted values, the residual being that
the brief band-averages `σ_u` over 1–6 m where this uses the 2 m pitch mean.

---

## 3. How the sweep measures it

`sweep_theory.py` flies each cell three ways on the **same release draws**:

```
still   no wind at all               → release scatter alone, the floor
mean    CFD mean field, no gusts     → floor + mean-wind bias and shear
full    mean field + sweeping gusts  → everything
```

and takes **paired per-throw displacements**, not differences of variances:

```
d_gust,i = landing(full,i) − landing(mean,i)
d_wind,i = landing(mean,i) − landing(still,i)
```

with `bias = mean(d)` (correctable — a thrower who knows the field aims off it)
and `σ = rms(d − mean(d))` (not correctable).

The variance-difference route was tried first and does not work at this site.
Release scatter alone is `σ ≈ 7 m` against a gust term of order 1–2 m, so
`σ_full² − σ_mean²` is a ~10 % difference of two variances each carrying a 20 %
relative error at n = 48 — the first smoke test returned a *negative* variance.
The paired displacement measures the same quantity with the shared release draw
cancelling exactly rather than statistically.

Common random numbers run across the whole speed axis: throw *i*'s release draw
and gust phases come from a seed sequence that excludes wind speed, heading,
bearing and arm, so the wind-speed response is a paired comparison.

Field statistics are taken **along the flown path**, arc-length weighted, not at
the release point — the release point is the most sheltered part of the flight
by construction, being at head height in the lee of the same buildings.

---

## 4. Files

| file | what it is |
|---|---|
| `wind_scaling.py` | `ScaledField`; the Reynolds-independence test and its verdict |
| `dispersion_law.py` | the ballistic law and the tracer alternative; symbolic exponent derivation |
| `sweep_theory.py` | the three-arm Monte Carlo over speed, direction, throw and spin |
| `analyse_theory.py` | the four questions, the error budget, the figures |
| `test_theory.py` | 16 tests, anchored on sympy limits and hand-built fields |
| `out/scale_invariance.json` | the §1 result |
| `out/sweep_speed.json` | the sweep |
| `out/theory_results.json` | the graded answers |

Run order:

```bash
bash ../P3_CFD/runScaleTest.sh            # once; ~40 min of CFD
python wind_scaling.py                    # the gate
python ../P1_AnalyticalModel/calibrate_catalogue.py
python sweep_theory.py --out out/sweep_speed.json
python analyse_theory.py --sweep out/sweep_speed.json
python -m pytest test_theory.py -q
```
