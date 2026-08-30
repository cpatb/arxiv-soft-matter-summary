"""Solve for the release speed that gives each catalogue throw its target range.

The range of a thrown disc is not a monotone-simple function of release speed
-- a disc released faster at fixed nose angle generates more lift, climbs, and
can turn over -- so this brackets first and bisects on the MEDIAN still-air
range of an ensemble, and reports the bracket it used. If the target range is
not bracketed inside the allowed speed range the entry FAILS rather than
returning the closest endpoint, because a huck silently pinned at the speed
ceiling would look like a calibrated throw and would not be one.

"Still air" here means literally zero wind, not the site's mean field: the
catalogue must be a property of the throw, transferable to any site.

Writes throw_catalogue.json. Deterministic: fixed seed, and the SAME release
draws are reused at every candidate speed (common random numbers), so the
range-vs-speed curve is smooth and the bisection is not chasing sampling noise.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import disc_6dof as d6
from coefficients import DEFAULT_COEFFS, ULTRASTAR
from trajectories import RELEASE_DISTRIBUTIONS
from throw_catalogue import _CATALOGUE_SPEC

# 60 draws, not 200. The bisection runs on COMMON RANDOM NUMBERS, so the
# range-vs-speed curve is a deterministic smooth function of speed and the
# bisection is not chasing sampling noise -- more draws buy a better estimate of
# the population median, not a better root. Measured cost at the top of the
# bracket is 0.7 s per integration, so 200 draws x 22 evaluations x 3 throws was
# 2.5 hours; this is ~25 minutes. The achieved median and its p10/p90 are
# reported from the same 60 draws, so the residual sampling error is visible.
N_DRAWS = 60
SEED = 20260830
SPEED_LO, SPEED_HI = 8.0, 34.0      # m/s, the bracket searched
RTOL, ATOL = 1e-6, 1e-8
NO_WIND = lambda r, t: np.zeros(3)   # noqa: E731


def median_range(speed, elevation, spin, base, draws):
    """Median still-air range over the shared release draws at this speed."""
    d = RELEASE_DISTRIBUTIONS[base]
    out = []
    for z_speed, z_elev, z_alpha, z_roll, z_spin in draws:
        r0, v0, q0, w0 = d6.build_initial_conditions(
            launch_pos=(0.0, 0.0, 1.5),
            speed=speed + z_speed * d["speed_std"],
            heading_deg=0.0,
            elevation_deg=elevation + z_elev * d["elevation_std"],
            release_alpha_deg=d["release_alpha_mean"] + z_alpha * d["release_alpha_std"],
            roll_deg=d["roll_mean"] + z_roll * d["roll_std"],
            spin_rad_s=max(20.0, spin + z_spin * d["spin_std"]),
            handedness=d["handedness"])
        res = d6.simulate(r0, v0, q0, w0, coeffs=DEFAULT_COEFFS, geom=ULTRASTAR,
                          wind=NO_WIND, t_max=12.0, rtol=RTOL, atol=ATOL,
                          compute_diag=False)
        r = np.asarray(res.r)
        out.append(float(np.hypot(r[-1, 0] - r[0, 0], r[-1, 1] - r[0, 1])))
    return float(np.median(out)), np.asarray(out)


def solve(name, spec, draws):
    f = lambda v: median_range(v, spec["elevation_mean"], spec["spin_mean"],  # noqa: E731
                               spec["base"], draws)[0]
    target = spec["target_range_m"]
    lo, hi = SPEED_LO, SPEED_HI
    flo, fhi = f(lo), f(hi)
    print(f"  bracket: {lo:.1f} m/s -> {flo:.1f} m   {hi:.1f} m/s -> {fhi:.1f} m "
          f"  (target {target:.0f} m)")
    if not (min(flo, fhi) <= target <= max(flo, fhi)):
        raise SystemExit(
            f"FAIL: {name}: target {target} m is not bracketed by "
            f"[{flo:.1f}, {fhi:.1f}] m over speeds [{lo}, {hi}] m/s at "
            f"elevation {spec['elevation_mean']} deg. Raise the elevation or "
            f"widen the bracket -- do NOT clip to an endpoint.")
    for _ in range(12):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if (fm - target) * (flo - target) <= 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
        # 0.05 m/s is ~3% of the release-speed spread the draws already
        # carry, so bisecting finer measures nothing.
        if hi - lo < 0.05:
            break
    speed = 0.5 * (lo + hi)
    med, allr = median_range(speed, spec["elevation_mean"], spec["spin_mean"],
                             spec["base"], draws)
    return dict(speed_mean=round(speed, 4),
                elevation_mean=spec["elevation_mean"],
                spin_mean=spec["spin_mean"],
                target_range_m=target,
                achieved_median_range_m=round(med, 3),
                range_p10_m=round(float(np.percentile(allr, 10)), 2),
                range_p90_m=round(float(np.percentile(allr, 90)), 2),
                n_draws=len(allr), seed=SEED,
                wind="none (still air)")


def main():
    rng = np.random.default_rng(SEED)
    # Common random numbers: standard normals reused at every candidate speed.
    draws = rng.standard_normal((N_DRAWS, 5))
    out = {}
    for name, spec in _CATALOGUE_SPEC.items():
        print(f"\n[{name}] target {spec['target_range_m']:.0f} m at "
              f"elevation {spec['elevation_mean']:.0f} deg")
        out[name] = solve(name, spec, draws)
        r = out[name]
        print(f"  -> speed_mean = {r['speed_mean']:.2f} m/s, "
              f"median range {r['achieved_median_range_m']:.1f} m "
              f"[p10 {r['range_p10_m']:.1f}, p90 {r['range_p90_m']:.1f}]")
    p = Path(__file__).with_name("throw_catalogue.json")
    p.write_text(json.dumps(out, indent=2))
    print(f"\nwrote {p}")


if __name__ == "__main__":
    main()
