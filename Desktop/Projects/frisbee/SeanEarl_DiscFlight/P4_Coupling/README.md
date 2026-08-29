# P4 — CFD-to-Trajectory Coupling and Monte Carlo

Wraps the P1 6-DOF integrator with:
1. A CFD wind-field interpolator that reads P3's `flow_fields.npz` and
   provides `U(r, t)` to the integrator.
2. A synthetic-turbulence gust generator that adds correlated
   `u'(t)` on top of the RANS mean, with local intensity from `k`.
3. A Monte Carlo driver that sweeps 10⁴ release-parameter realizations
   per (heading, launch position, throw type), records landing
   statistics, and feeds P5.

## Files

- `cfd_field.py` — trilinear interpolation of `flow_fields.npz` to arbitrary `(x, y, z)`.
- `synthetic_turbulence.py` — Mann (1998)-style isotropic Fourier sampler scaled to local TKE and integral length scale.
- `monte_carlo.py` — 10⁴-throw sweep driver; parallelizes over launch points × throw types.
- `landing_stats.py` — landing-distribution analyzer: 90% ellipse axes, success rate at canonical targets.

## Quick test (no CFD needed)

```bash
cd P4_Coupling
python -c "
import numpy as np
import sys; sys.path.insert(0, '../P1_AnalyticalModel')
from monte_carlo import sweep_one_launch_point
from cfd_field import StubField
field = StubField(uniform_wind=np.array([3.0, 0.0, 0.0]))
stats = sweep_one_launch_point(
    launch=(0, 0, 1.5), heading_deg=0,
    throw_type='right_backhand', n_throws=200, field=field, seed=0)
print(stats)
"
```

## Workflow with real CFD data

```bash
# after P3 produces flow_fields.npz:
python monte_carlo.py \
    --field ../P3_CFD/flow_fields.npz \
    --launch-grid -25 25 5 -25 25 5 \
    --headings 0 15 30 45 60 75 90 105 120 135 150 165 180 \
    --throws right_backhand right_forehand hammer \
    --throw-bearings 0 30 60 90 120 150 180 210 240 270 300 330 \
    --n-throws 1000 \
    --out landing_stats.npz

# landing_stats.npz is the input to P5
```

## Output schema (landing_stats.npz)

All arrays are indexed `[i_launch_x, i_launch_y, i_heading, i_throw_type,
i_bearing]` with any component axes appended. The authoritative, machine-
readable version of this table is written into every artifact as the
`schema_frames` entry; if the two ever disagree, `schema_frames` is right.

**Reference frames — these are not interchangeable.**

| entry | frame | population |
|---|---|---|
| `landing_mean[..., 3]` | **world** (absolute site x, y, z in m) | delivered |
| `landing_bias[..., 2]` | **throw** (downrange−25 m, crossrange), m | delivered |
| `target_xy[..., 2]` | **world**; `launch_xy + R·(cos b, sin b)` | — |
| `landing_cov[..., 2, 2]` | world-aligned axes; translation-invariant | delivered |
| `ellipse_axes[..., 2]` | 90% containment semi-major/minor, m | delivered |
| `success_rate[...]` | fraction within `target_radius_m` of `target_xy` | valid |

`landing_mean` is **not** an offset from the target. It is an absolute
position. An earlier version of this file described it as "(Δx, Δy) from
target", and a downstream reader who believed that reconstructed the target as
a fixed point at bearing 0 and measured 25–50 m of bias that does not exist.
Use `landing_bias`, or subtract the saved `target_xy`.

Populations: `valid` = reached the ground inside `t_max`; `delivered` = valid
and not landing inside a building footprint. Dispersion statistics use
`delivered`, because the 6-DOF model has no collision detection and the landing
point of a blocked throw describes no physical trajectory. `success_rate` uses
`valid`, because a throw into a building is a failed delivery rather than a
missing observation. `n_valid`, `n_delivered`, `n_airborne` and
`blocked_fraction` are all saved so the two denominators are reconstructable.

Scalars: `launch_x`, `launch_y`, `launch_z`, `headings`, `throws`, `bearings`,
`n_throws`, `seed`, `target_range_m`, `target_radius_m`, `field_source`,
`gusts_enabled`, `coeff_set`, plus the `prov_*` provenance block (git SHA,
library versions, SHA-256 of the input field, integrator tolerances, and the
realised gust-model parameters).

## Status (2026-07-26) — validated end-to-end on a stub field

`pytest test_coupling.py` — 18 tests passing. The whole P4→P5 chain was run
start-to-finish on a synthetic field before any CFD existed
(`stub_flow_fields.py` writes a `flow_fields.npz` in P3's exact schema).

### Three substantive fixes found by that smoke test

1. **The white-noise gust generator was unusable.** It drew a fresh random
   number inside every RHS evaluation, so the ODE right-hand side was
   non-deterministic. An adaptive integrator estimates local error by comparing
   stage evaluations; with a random RHS it measures the noise instead, rejects
   the step, and shrinks forever. Measured: **>100 s for a single throw**
   versus ~0.1 s without gusts — a 10⁴-throw sweep would never have finished.
   Worse, the result is not a solution of any well-posed equation: the landing
   scatter is a step-size artefact.
   Replaced with `CorrelatedGusts`, a Shinozuka–Deodatis spectral
   representation (sum of cosines with random phase) matching an exponential
   autocorrelation `exp(-Δt/T_L)` with `T_L = k/ε` from the RANS. It is
   analytic in `t`, so the integrator steps freely: gusts now cost ~0 extra
   (193 solver steps with gusts vs 222 without). Autocorrelation is asserted
   against the target to within 0.06.
2. **`success_rate` was in the documented schema but never computed.** P5 was
   falling back to a proxy, `min(1, (r_target/a_mean)²)`, which is not a
   probability — it ignores the offset between the landing centroid and the
   target, so a tight but badly-aimed group scored as a success. Now measured
   directly: the fraction of throws landing within 3 m of a 25 m target.
3. **Field interpolation dominated the sweep.** Five `RegularGridInterpolator`
   calls per RHS evaluation. Replaced with a hand-rolled trilinear that locates
   the cell once and shares the eight-corner weights across U, k and ε
   (`U_and_k` fuses the two the gust path needs). Asserted bit-equal to scipy.

Also: seeding now goes through `np.random.SeedSequence` spawning, so results are
reproducible and independent of worker count; the sweep is parallelised over
cells (~21 throws/s on 20 workers); out-of-domain queries are clamped rather
than extrapolated and the clamped fraction is reported.

### Measured throughput

~21 trajectories/s on 20 workers (0.95 s/throw/worker, varying by a factor ~4
between quiet and high-TKE cells). Budget sweeps accordingly.

## Open implementation tasks

1. **Synthetic turbulence**. The Mann-sampler stub generates white-in-time Gaussian noise with magnitude √(2k/3); to get real correlated gusts we need the von Kármán / Mann spectrum, which `synthetic_turbulence.py` only sketches. For first-pass Monte Carlo the white-noise approximation is fine.
2. **Vertical interp at z<1 m**. P3 samples at 1, 2, 3 m AGL. Below 1 m we either extrapolate via log-law or zero out; the integrator should never go below 0.3 m AGL anyway (release ≥ 1.5 m, landing at z=0 stops integration).
3. **Wind-direction at point**. The CFD heading at the inlet is θ_wind; at a given launch position, the local wind direction can differ substantially (wake, jet). The interpolator gives `U(r)` directly so no special handling is needed.
