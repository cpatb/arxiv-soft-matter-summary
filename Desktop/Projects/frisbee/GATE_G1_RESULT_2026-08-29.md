# GATE G1: the gust-model correction is real, and its predicted size is wrong

**2026-08-29.** WP1's instrument (`P4_Coupling/gust_model_comparison.py`) had been written
and never run. There was no `gust_model_comparison.json` on disk. It is run now, on the
real CFD field, and the gate does not pass.

## The result

Paired arms, common random numbers, 108 cells × 40 paired throws = 8,640 trajectories per
model. Every throw is flown twice from an identical release draw, once with gusts, once on
the RANS mean alone, so the estimator measures turbulence-induced displacement rather than
thrower scatter, which is ~98% of total landing variance.

```
sigma_turb(temporal)  = 1.7154 m per axis
sigma_turb(sweeping)  = 0.9886 m per axis
RATIO                 = 1.735
```

**Pre-registered acceptance: a fall of 2.0-2.6×. Measured 1.735×.**

The band is not merely missed, it is **excluded**. Cluster bootstrap over the 108 cells
(the correct resampling unit, throws within a cell share a launch point, heading, bearing
and field, with the same cell indices drawn for both arms, so the paired design is
preserved), 20,000 resamples:

| quantity | value |
|---|---|
| ratio | 1.735 |
| 95% CI | [1.638, 1.838] |
| bootstrap sd | 0.051 |
| P(ratio ≥ 2.0) | **0.0000** |
| P(ratio < 1.2) | 0.0000 |

## What it means

The brief derived the band from `√N_eddy` with N ≈ 4-6 independent eddies encountered
along the flight, and two independent panels put it at 2.3-2.6×. Inverting the measurement:

**N_eddy = ratio² = 3.01, 95% CI [2.68, 3.38].**

The disc integrates about **three** independent eddies over its flight, not four to six.
That is the one-line finding.

Both of the brief's failure branches are also excluded, so this lands in a case the brief
did not enumerate:

- It is **not** the "< 20% change, the defect was harmless" branch, P = 0.0000. The
  correlation-model defect was real and large: the old temporal model overstated
  turbulence-induced landing dispersion by 74%.
- It is **not** the "≫ 3×, the disc is not a linear filter" branch.

So the correction stands and the mechanism that motivated it survives, a sweeping,
arclength-correlated gust field really does disperse the disc far less than an
Eulerian-decay-time field. What fails is the *quantitative* prediction of how much. Per the
brief's own instruction not to proceed on a theory whose central prediction just failed,
**do not start WP2 on the 4-6 eddy count.** Re-derive N first; the measurement now
constrains it to 2.7-3.4.

## Robustness

Re-run with `--isotropic` (isotropic σ in both arms, so the comparison isolates the
correlation-model change with no anisotropy in it):

```
sigma_turb(temporal)  = 2.4338 m per axis
sigma_turb(sweeping)  = 1.3026 m per axis
RATIO                 = 1.868      95% CI [1.789, 1.955]   P(ratio >= 2.0) = 0.0018
```

Also below the band, and still excluding it at alpha = 0.05, though by a narrower margin
than the anisotropic run (P = 0.0018 against 0.0000). Implied N_eddy = 3.49, 95% CI
[3.20, 3.82].

**The verdict does not depend on the sigma choice.** Under both, the measured fall is
between 1.7x and 1.9x, the 2.0-2.6 band is excluded, and N_eddy comes out near 3 rather
than 4 to 6. Anisotropy accounts for some of the gap but not the sign of it: switching the
gust field from isotropic to sigma_v/sigma_u = 0.75, sigma_w/sigma_u = 0.5 moves the ratio
from 1.87 to 1.74, which is the wrong direction and too small to reach 2.0.

One caveat worth recording: the isotropic temporal arm dropped 7 of 4,320 pairs (a throw in
one arm did not reach the ground within t_max). The anisotropic run dropped none. 7 pairs
in 4,320 cannot move a ratio of this size, but it is the only asymmetry between the arms.

Run it with:

```
python gust_model_comparison.py --field ../P3_CFD/flow_fields.npz \n    --n-throws 40 --jobs 20 --isotropic --out gust_model_comparison_isotropic.json
```

## Reproduce

```
cd P4_Coupling
python gust_model_comparison.py --field ../P3_CFD/flow_fields.npz \
    --n-throws 40 --jobs 20 --out gust_model_comparison.json
```
~16 min on 20 workers. No cells dropped, no pairs lost, no scale clipping in either arm.

## The separate blocker this run exposed

`P4_Coupling/landing_stats.npz` is dated 2026-07-26 03:19; `monte_carlo.py` is dated
14:19 the same day: **eleven hours later**. The artifact therefore predates the frame fix
and carries none of the schema the current `P4_Coupling/README.md` documents: no
`landing_bias`, no `target_xy`, no `n_delivered`, no `schema_frames`, and no `prov_*` block.
Its keys are:

```
alpha_landing, bearings, blocked_fraction, coeff_set, ellipse_axes, field_source,
headings, landing_cov, landing_mean, launch_x, launch_y, launch_z, n_throws,
n_valid, seed, success_rate, target_radius_m, target_range_m, throws
```

Consequences:

1. `landing_stats.py`'s new guard will refuse to summarise it, by design.
2. **There is no record of which gust model built it**, `gust_model` is not among the saved
   scalars, so it cannot be established from the artifact whether the manuscript's Q rests
   on `temporal` or `sweeping`. Given the file predates the rewrite, `temporal` is the
   likely answer, and G1 has just measured that as a 74% overstatement.
3. `P5_Synthesis/Q.npz` derives from it, and every manuscript macro derives from Q.

**The fix is one command and about 1.5 hours of compute**, and it regenerates the whole
chain consistently:

```
bash SeanEarl_DiscFlight/runDownstream.sh 20 150
```

97,200 trajectories at the ~19.5 traj/s measured here. It was not run in this session
because it moves every downstream number, and the paired-significance verdict
(`worst_supported: 270`) could move with them, that is a result the author should be
present for.
