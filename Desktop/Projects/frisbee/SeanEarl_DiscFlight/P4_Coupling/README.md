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

Keyed array entries:
- `success_rate[i_launch, i_heading, i_throw_type, i_bearing]`  — 0–1
- `landing_mean[..., 2]`        — (Δx, Δy) from target, m
- `landing_cov[..., 2, 2]`      — covariance of landing distribution, m²
- `ellipse_axes[..., 2]`        — 90% containment semi-major/minor, m
- `launch_grid[..., 3]`         — launch point coords
- `headings_deg[...]`
- `throw_types[...]`
- `bearings_deg[...]`

## Open implementation tasks

1. **Synthetic turbulence**. The Mann-sampler stub generates white-in-time Gaussian noise with magnitude √(2k/3); to get real correlated gusts we need the von Kármán / Mann spectrum, which `synthetic_turbulence.py` only sketches. For first-pass Monte Carlo the white-noise approximation is fine.
2. **Vertical interp at z<1 m**. P3 samples at 1, 2, 3 m AGL. Below 1 m we either extrapolate via log-law or zero out; the integrator should never go below 0.3 m AGL anyway (release ≥ 1.5 m, landing at z=0 stops integration).
3. **Wind-direction at point**. The CFD heading at the inlet is θ_wind; at a given launch position, the local wind direction can differ substantially (wake, jet). The interpolator gives `U(r)` directly so no special handling is needed.
