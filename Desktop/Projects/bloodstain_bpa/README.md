# Bloodstain Pattern Analysis — ballistic region-of-origin

A self-contained analytical + simulation project on **the drop impact of a bodily
fluid**. It is the forensic inverse problem of your drop-impact work: instead of
predicting a splat from an impact, it reconstructs the impact (and its 3-D origin)
from the splat.

## The physics, in one paragraph

A bloodstain is an ellipse. Three measurements decode the impact:

1. **Impact angle** `α = arcsin(W/L)` from the ellipse aspect ratio (Balthazard 1939).
2. **Direction of travel** from the major-axis azimuth + tail.
3. **Impact speed** from the stain *size*, via a **drop-impact spread factor**
   `ξ = D_stain/D₀ = f(We_n)` (Pasandideh-Fard et al. 1996). This is the bridge to
   your drop-impact simulations. Because `ξ(We_n)` is monotonic, measuring the stain
   size *measures the normal impact speed*.

Together these give the full **impact velocity vector** at each stain. Integrating
the equations of motion (gravity + Schiller–Naumann drag) **backward in time** traces
each drop's flight path; the point where all paths converge is the region of origin.

The classical forensic **straight-line "stringing"** method is the special case that
ignores gravity, and is known to *overestimate* the height of origin. This project
quantifies that bias and shows the ballistic method removes it.

## What's here

| file | purpose |
|------|---------|
| `bpa.py` | core module: spread model, ballistic flight (`solve_ivp`), forward synthetic-stain generator, straight-line and ballistic reconstructions |
| `test_bpa.py` | 6 correctness tests (spread invertibility, free-fall vs analytic, forward/backward reversibility, round-trip origin recovery) |
| `run_experiments.py` | 4 figures + numeric report |
| `uncertainty.py` | **new 2026-08-07**: error propagation, bias attribution, posterior over the region of origin |
| `test_uncertainty.py` | 17 checks on the error layer |
| `run_uncertainty.py` | 4 further figures + the error-rate report |
| `figures/` | output |

Reproduce: `python -m pytest test_bpa.py test_uncertainty.py -q`, then
`python run_experiments.py` and `python run_uncertainty.py` (the latter takes a few
minutes; every posterior draw is a full back-projection).

## Key findings

**The ballistic method removes the classical bias.** Straight-line stringing overestimates
the height of origin by ~1.3–2.2× for this impact-spatter-onto-floor geometry, and the bias
grows with launch speed, because faster drops fly flatter and accumulate more gravitational
curvature over the flight. With matched physics and no measurement error, the ballistic
reconstruction recovers the known origin to a few centimetres.

**But the stain-size route to impact speed is a ~3.5× error amplifier**, and this was not
stated anywhere before. The spread correlation is roughly `ξ ~ v^(1/3.5)` at BPA speeds, so

```
d ln(v_n) / d ln(ξ) = 3.1 to 3.9   over v_n = 1 to 16 m/s
```

A 5% error in the stain minor axis is an ~18% error in impact speed; combined with 10%
uncertainty in the parent drop diameter it is ~38%. This is the dominant methodological
fact about the method and it belongs in any report that uses it.

**The amplifier makes the height estimate biased, not merely noisy.** Because `v ~ ξ^k`
with `k >> 1`, symmetric measurement noise on `ξ` gives `E[v] > v_true`, and over-estimated
speeds back-project to flatter, longer trajectories that converge lower. Attributing the
bias one noise source at a time (`figU2`):

| noise source | height bias | scatter |
|---|---|---|
| stain width, 5% | +3.9 cm | 15.4 cm |
| stain length, 5% | +1.4 cm | 4.0 cm |
| azimuth, 3° | +3.1 cm | 0.3 cm |
| **parent drop D₀, 10%** | **−17.2 cm** | 23.1 cm |
| all combined | −22.1 cm | 26.2 cm |

**Parent-drop-size uncertainty carries essentially the whole bias**, which is awkward
because `D₀` is the least well known quantity in casework.

⚠ **An earlier version of this README claimed the ballistic estimate "stays near-unbiased"
under noise while the run output showed −18 cm on a 1.15 m height.** The output was right.

**How well must `D₀` be known?** (`figU3`) The estimator is unbiased only while
`σ_D₀ ≲ 4%`; at 10% it under-states the height by ~13 cm and at 20% by ~35 cm.

**The deliverable is an interval, not a point.** For the reference scene with the full
noise model, the 95% credible interval on the height of origin is **[0.46, 1.33] m** against
a true 1.15 m. The truth is inside, so the posterior is honest, but the interval is 87 cm
wide. That is the error rate the 2009 NAS and 2016 PCAST reports asked BPA for, and it is
the number a report should quote.

**A correction that was tried and rejected.** A first-order Jensen shrinkage of `ξ` makes
the bias worse (−22 cm → −33 cm), because `D₀` enters the inversion in three places at once
and also sets the drag, so a correction derived from the `ξ` pathway alone has the wrong
magnitude. It is kept in the code, off by default, so the negative result stays
reproducible.

## Honest limitations / where to go next

- **The exact bias magnitude is geometry-dependent** (emission cone, floor vs. wall,
  source height). The *direction and mechanism* of the bias are robust; the percentage
  is illustrative of this scenario, not a universal constant.
- **Parent drop diameter D₀** is needed to convert stain size → speed. Here it's known
  (synthetic) or assumed; in casework D₀ is uncertain. The spine-count estimator
  `n ≈ 1.14·We^{1/2}` is included as a size-independent cross-check but is noisy.
- **Spread model** is the smooth-surface maximum-spread correlation; real surfaces
  (fabric, textured drywall) deviate. A surface-dependent `ξ(We_n, surface)` calibration
  is the obvious extension, and a natural place to *reuse your actual drop-impact
  simulation data* to build a blood-specific spread law.
- **Next steps**: (1) wall + floor combined patterns and true 3-D convergence clouds;
  (2) Monte-Carlo posterior over the region of origin instead of a point estimate;
  (3) validate against a public BPA benchmark image set; (4) couple to a Basilisk blood-
  drop impact sim to replace the analytic spread law with a simulated one.

## References (method anchors)

- Balthazard et al. (1939) — impact-angle/ellipse relation.
- Pasandideh-Fard, Qiao, Chandra & Mostaghimi, *Phys. Fluids* **8**, 650 (1996) — max spread.
- Schiller & Naumann (1935) — sphere drag correlation.
- Attinger et al., *Forensic Sci. Int.* (2013) — fluid dynamics of BPA (review).
