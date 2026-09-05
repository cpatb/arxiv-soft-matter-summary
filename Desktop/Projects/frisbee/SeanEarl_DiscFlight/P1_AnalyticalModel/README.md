# P1 — Analytical 6-DOF Disc Model

Pure-Python implementation of the 13-state quaternion-based 6-DOF rigid-body
model derived in `Project1_AnalyticalTreatment.tex`.

## Files

- `coefficients.py` — Hummel (2003), Hummel-long, Potts–Crowther coefficient sets and disc geometry.
- `disc_6dof.py` — core integrator (right-hand side, quaternion kinematics, aero force/moment, simulate wrapper).
- `trajectories.py` — canonical throws (backhand, forehand, hammer) + sampling distributions for Monte Carlo.
- `validation.py` — four-case scenario runner that prints summary statistics.
- `test_disc_6dof.py` — pytest sanity tests (range, drift, spin decay, rotation correctness).

## Quick start

```bash
pip install numpy scipy pytest matplotlib
python validation.py
pytest test_disc_6dof.py -v
```

## State vector

13 components: `r(3), v(3), q(4), omega_B(3)`.
- `r`, `v` in inertial frame (z up).
- `q = (qw, qx, qy, qz)`, body-to-inertial rotation quaternion.
- `omega_B` body-frame angular velocity; `omega_B[2]` is the spin rate.

## Wind interface

`wind(r, t) -> np.ndarray(3)` returns the inertial-frame wind vector at the
disc's position. P4 will replace `still_air` / `uniform_wind` with a
CFD-field interpolator from OpenFOAM output. No change to the integrator
required — just a different `wind` callable.

## What's tested

- Quaternion → rotation-matrix correctness (identity + 90° z).
- Still-air range falls in physical ballpark.
- Headwind reduces range.
- Crosswind drift matches `U_c * T` to within ~50%.
- Spin decay over a flight is < 30%.

## Validation status (2026-07-25)

`validation_hummel.py` compares the integrator to two independently published,
fully-specified results in Hummel (2003) — the same source that supplied the
coefficients, so agreement tests the *integrator and frame conventions* rather
than the coefficient fit. Both are asserted in `test_disc_6dof.py`.

**Case A — "common throw", Hummel Sec. 2.4.3 / Fig. 2-9.** Release conditions
are given exactly (v=14 m/s, spin=50 rad/s, θ=11°, β=6°, α=5°, z₀=1 m):

| quantity | Hummel (published) | this model | error |
|---|---|---|---|
| range | ~25 m | 26.99 m | +8% |
| apex height | 2.55 m | 2.45 m | −4% |
| flight time | 3.5 s | 3.48 s | −1% |
| terminal speed | ~v₀/3 | 4.83 m/s | +4% |
| lateral | right, reverses ~3 s | right 1.21 m, reverses 3.34 s | ✓ |

**Case B — steady glide, Hummel Sec. 2.4.2** (analytic, no integration):

| quantity | Hummel | this model |
|---|---|---|
| trim α (CM=0) | 9° | 9.00° |
| glide angle γ | −19.3° | −19.32° |
| C_L at trim | 0.67 | 0.672 |
| trim speed | 9.1 m/s | 8.28 m/s (−9%) |

Only the **tabulated** Potts–Crowther set reproduces these; the linear sets miss
range/apex/time by 30–57% and never fade back left. It is therefore the default.

Figure: `validation_hummel.png`.

## Corrections made during validation

Five substantive discrepancies against the primary source were found and fixed.
Ground truth is Hummel's own MATLAB (`discfltEOM.m`, thesis Appendix B), which
is the code the Table 3-2 coefficients were *fitted with* — the values are only
meaningful inside those exact functional forms.

1. **`CL` and `CM` used `(α − α₀)`** where Hummel uses `α` itself. α₀ appears
   **only** in the drag parabola. This shifted the zero-pitching-moment angle
   from 8.8° to 4.8° — the single number that sets where turn flips to fade.
2. **α₀ was hardcoded to −4°.** It is a derived quantity, `α₀ = −CL0/CLa`,
   which is −9.9° for the Hummel short-flight set (−3.9° for Potts–Crowther,
   which is where the −4° came from).
3. **Rate-damping terms were non-dimensionalised by `d/(2V)`.** Hummel's
   `CMq`, `CRp` multiply *raw* rates in rad/s, and `CRr` carries an explicit
   `sqrt(d/g)` time scale. The normalisation under-damped pitch and roll by
   ~2 orders of magnitude.
4. **`CRr` value mismatch.** Table 3-2 lists `1.4e-2`, but Appendix A carries
   `CRr = 0.00171` with the note *"for nondimensionalization = sqrt(d/g)"*.
   Only the latter is consistent with the implemented form.
5. **Moments of inertia used the thin-uniform-disc idealisation**
   (Is=1.58e-3, Ip=7.9e-4) instead of Hummel's measured Table 3-1 values
   (Is=2.35e-3, Ip=1.22e-3) — low by 33% and 35%. Since precession goes as
   `M/(Is·ω_z)`, this overstated turn and fade by ~50%.

Also fixed: `disc_6dof.py` had a duplicated fragment after the final `return`
that made the module a syntax error, so nothing in P1 imported.

`HUMMEL_LONG` was additionally corrected: Hummel Sec. 3.3.4 states that for the
long flights only `CMq`, `CRp`, `CNr` were re-fitted and the other seven were
taken from Potts & Crowther. The previous set paired the short-flight force
coefficients with long-flight damping, a hybrid matching no published result.

## Known model limitations (carry into the paper)

- **No side-force coefficient `CS(α)`** for rim asymmetry. Consequence: the
  model has no yaw stiffness, so launch sideslip is shed only slowly, and predicted
  crosswind drift comes out at about a quarter of the naive `U_c·T` full-advection
  estimate at 5 m/s. Do not quote crosswind drift as well-validated. (The
  Galilean-invariance test confirms the wind coupling itself is exact, so this is a
  physics-model gap, not a numerical bug.)
  - ⚠ **RETRACTED 2026-09-03 (WP4b item 5): this used to say crosswind drift is a
    "lower bound", and that conflates two different omissions.**
  - ⚠ **AMENDED 2026-09-04:** the 09-03 retraction itself said "missing yaw stiffness
    is one-sided — with `CS(α)` the disc would weathercock and shed sideslip faster."
    **That is withdrawn too.** Yaw stiffness means the yawing moment depends on
    *sideslip angle*; for a body of revolution the geometry is identical at every
    azimuth, so no such dependence exists. An α-dependent side force is not a
    β-dependent yawing moment — the model does not *omit* yaw stiffness, an
    axisymmetric disc **has none**. And a disc spinning at 90–110 rad/s
    (`Is = 2.35e-3 kg m²`, `L ≈ 0.21–0.26 kg m²/s`) responds to a moment by
    **precessing perpendicular to it**, not by rotating into the wind, so it cannot
    weathercock at any coefficient. Raised in `CLAUDE_CODE_FRISBEE_HARDENING.md` §(a).
    Net effect: **both** missing terms are unsigned, which makes the conclusion below
    follow from symmetry rather than from a race between a one-sided and a two-sided
    term.
  - The model also omits the **Robins–Magnus** side force, which goes as
    `ω × v` and therefore **reverses with the spin direction**: a backhand and a
    forehand spin opposite ways, so the same crosswind gives a lateral force of
    opposite sign for the two throws. The net effect is handedness- and
    throw-type-dependent and can add to the modelled drift or subtract from it.
    **Lateral predictions are uncertain in sign, not bounded below.** The same
    retraction is applied at two places in the manuscript.
- **Angle-of-attack range.** The tabulated curves span −10°…+30°; a realistic
  flight reaches ~33° in its final second, so samples in the landing phase are
  clamped to the table edge. The linear sets are far worse (44% out of range, α to 47°).
  - ⚠ **CORRECTED 2026-09-04.** This used to say "~7% of samples are clamped". That
    figure does not reproduce on either reading, and the true exposure is larger.
    Measured by `check_alpha_coverage.py` (seed 20260904, n = 200 draws per throw):

    | throw | nominal | ensemble mean | ensemble α_max | draws never leaving the table |
    |---|---|---|---|---|
    | `right_backhand` | **17.6%** | **12.0%** | 30.0° [5.1, 52.6] | 52% |
    | `right_forehand` | 0.0% | 8.8% | 25.5° [2.1, 51.9] | 65% |
    | `hammer` | 0.0% | 0.1% | 6.2° [−3.0, 27.1] | **100%** |

    Two things to carry: the backhand's exposure is 12–18%, not 7%; and **the hammer is
    the *cleanest* throw by this measure, not the worst** — it never leaves the table in
    200 draws. The manuscript said the opposite until 2026-09-04 (see
    `wp4b_item6_retract.py`). The `hammer()` roll of +80° is what keeps its α small.
- No Robins–Magnus side force; Hummel notes its omission may bias other
  coefficients.

## Remaining (not blocking the paper)

1. Digitise Potts–Crowther's own figures directly rather than via Hummel's
   transcription, and extend the tables past 30° to cover the landing phase.
2. Add advance-ratio dependence to `CL`/`CD` (currently only `CRr` carries it).
3. Add `CS(α)` and re-check crosswind drift against field measurements.
