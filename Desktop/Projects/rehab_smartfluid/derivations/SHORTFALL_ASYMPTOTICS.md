# Closed forms for the quasi-steady shortfall

**Date:** 2026-09-01 · **Script:** `derivations/shortfall_asymptotics.py`
All five sections below regenerate with `python derivations/shortfall_asymptotics.py`.
Every symbolic identity is asserted inside the script, so a broken one raises rather
than printing.

`SCALING_2026-08-31.md` §4 gave the shortfall as a table of integrated numbers. A
table cannot be checked and does not explain itself. This replaces it with closed
forms, and then measures where those closed forms stop being accurate instead of
assuming they hold.

---

## 1. Setup

Writing `D = phi_0 - phi_m` for the width of the jamming interval, the jamming
fraction is linear in the frictional fraction,

    phi_J(f) = (1 - f) phi_0 + f phi_m = phi_0 - f D

and the viscosity is `eta = eta_s g^n` with `g = phi_J / (phi_J - phi)`.

The amplification of a change in frictional fraction follows by differentiating:

    d ln eta / d f  =  n D phi / [ (phi_0 - f D)(phi_0 - f D - phi) ]

verified symbolically against SymPy's own derivative. **This is the quantity that
makes the shortfall worse near jamming**: both factors in the denominator shrink as
`phi` rises toward `phi_J`, so the same change in `f` moves the viscosity much
further.

## 2. The approximation, stated once

The exact problem is nonlinear: `f` grows, so viscosity grows, so stress grows, so
`f_ss` grows. Freeze that feedback by evaluating `f_ss` once at the **steady**
stress. The strain clock then integrates exactly,

    f(A) = f_s (1 - exp(-A)),          A = gamma_dot tau / gamma_0

and everything below is closed form.

The true stress starts at the lubricated value and climbs, so the true `f_ss` is
smaller than the frozen one throughout. **The frozen result is therefore an upper
bound on the delivered stress**, and §4 checks that the bound is never violated
rather than trusting the argument.

## 3. Results

**The A -> 0 plateau is exactly `eta(0)/eta(f_s)`**, equivalently
`sigma_lubricated / sigma_steady`. Both sides come from steady quantities with no
integration at all, which makes this an anchor outside the code under test.

| phi | integrated | closed form | rel. diff |
|---|---|---|---|
| 0.40 | 0.752675 | 0.752674 | 3.1e-07 |
| 0.46 | 0.405846 | 0.405845 | 5.5e-07 |
| 0.50 | 0.140566 | 0.140565 | 9.9e-07 |
| 0.52 | 0.043145 | 0.043145 | 1.3e-06 |

**The initial slope** is the plateau times the amplification at `f = 0`:

    d(ratio)/dA |_0  =  [eta(0)/eta(f_s)] * n D phi f_s / [phi_0 (phi_0 - phi)]

matching a finite difference of the closed form to 1e-7 at all four packings
(0.1847, 0.2580, 0.1496, 0.0578).

**The crossover** `A*` at which the ratio reaches `1 - eps` inverts in closed form.
Solving `g(f) = (1-eps)^(1/n) g(f_s)` for `f` gives

    f* = [ phi_0 - G phi / (G - 1) ] / D,     G = (1-eps)^(1/n) g(f_s)
    A* = -ln(1 - f*/f_s)

| phi | A* closed (bound) | A* integrated |
|---|---|---|
| 0.40 | 3.485 | **4.585** |
| 0.46 | 4.860 | **6.908** |
| 0.50 | 5.963 | **7.491** |
| 0.52 | 6.797 | **7.733** |

The closed form sits below the integrated value at every packing, which is what an
upper bound on the ratio must do. **This supersedes the "safe above A ~ 10" in
`SCALING` §4**: at the 1% level the crossover is 4.6 to 7.7 and it rises with
packing fraction.

## 4. Where the closed form stops being good

Gap between the frozen bound and the full integration, in ratio units:

| phi \ A | 0.1 | 0.3 | 1.0 | 3.0 | 10.0 |
|---|---|---|---|---|---|
| 0.40 | 0.004 | 0.011 | 0.024 | 0.016 | 0.0002 |
| 0.46 | 0.007 | 0.023 | 0.071 | 0.077 | 0.001 |
| 0.50 | 0.004 | 0.014 | 0.063 | **0.123** | 0.001 |
| 0.52 | 0.001 | 0.005 | 0.031 | **0.121** | 0.001 |

The gap is not monotone in `A`. It vanishes at both ends, because at `A -> 0`
nothing has happened yet and at `A -> 10` both have converged, and peaks near
`A ~ 1` to `3` where the stress feedback has time to act but has not finished.
Worst case 0.1225 at `phi = 0.50, A = 3`.

**Read this as a measurement, not a defect.** The gap *is* the strength of the
stress feedback. Where it is small the shortfall is a kinematic effect of the
strain clock alone; where it is 0.12 the nonlinearity is doing real work, and a
device operating there cannot be designed with any closed form.

## 5. What this changes

1. `SCALING` §4's crossover claim is replaced by the table in §3 above.
2. The plateau closed form becomes a regression anchor in the test suite: it
   checks the integrator against steady-state algebra rather than against itself.
3. The amplification expression in §1 is the analytic statement of the numerical
   finding that the shortfall worsens with `phi`. It is now a formula rather than
   four numbers, so it extrapolates.

## 6. What this does not do

- It assumes the constant-rate step of `SCALING` §7.3. A gait waveform ramps.
- It inherits every limitation in `SCALING` §7, in particular that Wyart-Cates has
  no microstructure variable and cannot describe flow reversal.
- `n = 2` and `beta = 1` throughout. The amplification expression carries `n`
  symbolically, so its sensitivity is readable, but nothing here was tested at
  another exponent.
