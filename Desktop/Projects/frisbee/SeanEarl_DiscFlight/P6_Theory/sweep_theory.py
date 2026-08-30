"""Monte Carlo sweep over WIND SPEED, wind direction, throw direction and throw
type -- the data the general theory is tested against.

WHAT THE SWEEP IS FOR
---------------------
The production sweep answers "how throwable is each part of this field, in the
three wind directions we solved, at 5 m/s". This one answers the question the
theory needs: how does throwability depend on WIND SPEED, and on the flight
itself, and is that dependence the one `dispersion_law.py` predicts?

Wind speed is a free parameter because `wind_scaling.py` measured that it is:
the RANS solution at 10 m/s reproduces the 5 m/s solution to a median 0.06% of
RMS over the pitch after the expected powers are divided out, and the three
derived quantities the theory uses (sigma_u ~ Uref, T_L ~ 1/Uref, L ~ Uref^0)
hold to within 0.0005 in the exponent. So one solve per direction covers the
whole speed axis.

THREE ARMS, AND WHY
-------------------
Each cell is flown three ways on the SAME release draws:

    still   no wind at all              -> release scatter alone, the floor
    mean    CFD mean field, no gusts    -> floor + mean-wind bias and shear
    full    mean field + sweeping gusts -> everything

The gust and mean-wind contributions are then measured as PAIRED PER-THROW
DISPLACEMENTS, not by subtracting variances:

    d_gust,i = landing(full, i) - landing(mean, i)
    d_wind,i = landing(mean, i) - landing(still, i)

with `sigma_gust = rms(d_gust - mean(d_gust))` and the mean of `d_wind` being
the correctable aiming bias.

Subtracting variances instead would not work here, and the reason is worth
recording. At this site the release scatter alone is sigma = 6.96 m while the
gust contribution is predicted at 1.6-2.8 m, so recovering the gust term from
sigma_full^2 - sigma_mean^2 means taking a ~10% difference of two variances,
each estimated with a relative error of sqrt(2/n) = 20% at n = 48. The first
smoke test did exactly that and returned a NEGATIVE variance -- the estimator
was 2x noisier than the signal. The paired displacement measures the same
quantity directly, with the shared release draw cancelling exactly rather than
statistically.

COMMON RANDOM NUMBERS ACROSS WIND SPEED
---------------------------------------
The RNG for throw i is derived from a seed sequence that EXCLUDES the wind
speed, the heading and the arm. So the same release and the same gust phases
are used at every speed, and the wind-speed response is measured on paired
throws rather than across independent samples. This is the same design the
production sweep uses across headings, and for the same reason: ~98% of
per-throw variance is the shared release draw, so an unpaired speed sweep would
be swamped by it.

WHAT A THROW MEASURES
---------------------
Per trajectory: landing point, flight time T_f, path length L_f (the actual arc,
not the chord -- the law's L_f is the distance flown through the turbulence),
mean airspeed V_d = L_f/T_f, apex height, and whether it ended inside a
building. Everything the law needs is measured on the same flights whose scatter
it is asked to predict, so no quantity is carried over from another calculation.
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "P1_AnalyticalModel"))
sys.path.insert(0, str(ROOT / "P4_Coupling"))
sys.path.insert(0, str(ROOT / "P5_Synthesis"))

DEFAULT_FIELD = ROOT / "P3_CFD" / "vertical_fields.npz"

# Wind speeds at the 10 m reference. 0 is the still-air control; 20 m/s at 10 m
# is a strong gale (Beaufort 8), well past the point anyone plays in.
U_REFS = (0.0, 1.0, 2.0, 3.5, 5.0, 7.0, 10.0, 14.0, 20.0)
HEADINGS = (90, 210, 270)
THROWS = ("pass_15m", "right_backhand", "huck_50m", "pull_65m")
BEARINGS = (0.0,)                 # due east, the production sweep's bearing
N_THROWS = 48
T_MAX = 15.0
SEED = 20260830

_FIELD_CACHE = {}


def _base_field(field_path, heading):
    key = (str(field_path), int(heading))
    if key not in _FIELD_CACHE:
        from cfd_field import CFDField
        _FIELD_CACHE[key] = CFDField(str(field_path), heading=int(heading))
    return _FIELD_CACHE[key]


def _path_field_stats(base, u_ref, r, seg, stride=4):
    """Arc-length-weighted mean sigma_u and integral length along a trajectory.

    `base` is the unscaled field; the scaling is applied here rather than
    through a `ScaledField` wrapper because sigma_u ~ s and L ~ s^0 exactly, so
    only one multiply is needed and the per-point wrapper call is not.

    `stride` subsamples the integrator's own points. The field is trilinear on
    a 1-2 m grid, so sampling every 4th point of a ~250-point trajectory still
    resolves it, and this runs inside the innermost loop of the whole sweep.
    """
    s = u_ref / 5.0
    idx = np.arange(0, len(r) - 1, stride)
    if len(idx) == 0:
        return float("nan"), float("nan")
    w = np.add.reduceat(seg, idx)
    sig = np.empty(len(idx))
    Li = np.empty(len(idx))
    for j, i0 in enumerate(idx):
        p = r[i0]
        k = max(float(base.k(p, 0.0)), 0.0) * s ** 2
        eps = max(float(base.eps(p, 0.0)) * s ** 3, 1e-12)
        su = np.sqrt(2.0 / 3.0 * k)
        sig[j] = su
        Li[j] = su * (k / eps)
    W = w.sum()
    if W <= 0:
        return float(sig.mean()), float(Li.mean())
    return float((sig * w).sum() / W), float((Li * w).sum() / W)


def _release_rng(i):
    """RNG for throw i, independent of wind speed, heading, bearing and arm."""
    return np.random.default_rng(
        np.random.SeedSequence(entropy=SEED, spawn_key=(int(i),)))


def fly_cell(heading, u_ref, throw, bearing, launch, field_path,
             n_throws=N_THROWS, arms=("still", "mean", "full"), spin=None,
             release_speed=None, elevation=None):
    """Fly one cell and return per-arm per-throw measurements.

    `spin` overrides the catalogue entry's mean spin rate, in rad/s, leaving
    everything else alone. A disc is a gyroscope: a gust's rolling moment
    precesses into a heading change a quarter cycle later, so a faster-spinning
    disc should be stiffer against exactly the perturbation the dispersion law
    describes. Whether that is true is question 3 of
    `FRISBEE_GENERALIZATION.md` and, since the ballistic law as derived contains
    no spin at all, it is also the first candidate explanation for any
    disagreement between that law and the measured dispersion.
    """
    import disc_6dof as d6
    import monte_carlo as mc
    import synthetic_turbulence as st
    from throw_catalogue import build_catalogue
    from wind_scaling import ScaledField

    dist = dict(build_catalogue()[throw])
    if spin is not None:
        dist["spin_mean"] = float(spin)
    # Release speed and elevation are overridable so the flight LENGTH and
    # AIRSPEED can be varied without a calibrated catalogue entry. The
    # dispersion law is stated in terms of the flight that actually happened
    # (measured T_f, V_d, L_f), so a throw does not need a target range
    # attached to it to be useful here -- it only needs to be a different
    # flight. Varying speed and elevation on separate axes also breaks the
    # collinearity between L_f and V_d that a speed-only sweep would have.
    if release_speed is not None:
        dist["speed_mean"] = float(release_speed)
    if elevation is not None:
        dist["elevation_mean"] = float(elevation)
    base = _base_field(field_path, heading)
    launch = np.asarray(launch, dtype=float)

    def make_wind(arm, rng):
        if arm == "still" or u_ref == 0.0:
            return (lambda r, t: np.zeros(3))
        f = ScaledField(base, u_ref)
        # The gust generator is built for BOTH wind arms, even though the mean
        # arm suppresses its output: `CorrelatedGusts` draws its phases from
        # `rng`, and building it only on one arm would leave the two arms with
        # different downstream RNG state. `wind_mean_only` takes it and ignores
        # it, which is exactly what the pairing needs.
        g = st.CorrelatedGusts(f, rng, t_max=T_MAX, ref_pos=launch,
                               model="sweeping", sigma_ratios=st.SIGMA_RATIOS)
        return (st.wind_mean_only(f, g) if arm == "mean"
                else st.wind_with_gusts(f, g))

    out = {a: dict(land=[], T_f=[], L_f=[], V_d=[], apex=[], blocked=[],
                   path_sigma_u=[], path_L_int=[])
           for a in arms}
    for i in range(n_throws):
        rng = _release_rng(i)
        r0, v0, q0, w0 = d6.build_initial_conditions(
            launch_pos=launch,
            speed=rng.normal(dist["speed_mean"], dist["speed_std"]),
            heading_deg=bearing,
            elevation_deg=rng.normal(dist["elevation_mean"], dist["elevation_std"]),
            release_alpha_deg=rng.normal(dist["release_alpha_mean"],
                                         dist["release_alpha_std"]),
            roll_deg=rng.normal(dist["roll_mean"], dist["roll_std"]),
            spin_rad_s=max(20.0, rng.normal(dist["spin_mean"], dist["spin_std"])),
            handedness=dist["handedness"])
        # The gust phases are drawn from a SEPARATE stream seeded the same way
        # for every arm and speed, so the release draw above is not consumed
        # differently between arms.
        for arm in arms:
            grng = np.random.default_rng(
                np.random.SeedSequence(entropy=SEED + 1, spawn_key=(int(i),)))
            res = d6.simulate(r0, v0, q0, w0, coeffs=mc.DEFAULT_COEFFS,
                              geom=mc.ULTRASTAR, wind=make_wind(arm, grng),
                              t_max=T_MAX, rtol=mc.MC_RTOL, atol=mc.MC_ATOL,
                              compute_diag=False)
            r = np.asarray(res.r)
            seg = np.linalg.norm(np.diff(r, axis=0), axis=1)
            T_f = float(res.t[-1])
            L_f = float(seg.sum())
            o = out[arm]
            # Field statistics along the PATH the disc actually flew, not at
            # the release point. The law's sigma_u and L are properties of the
            # turbulence the disc passes through; a 30 m throw from a sheltered
            # corner crosses several times the variation in k that a
            # release-point probe would ever see, and the release point is the
            # most sheltered part of the flight by construction (it is at head
            # height in the lee of the same buildings). Arc-length weighted, so
            # the integrator's variable step does not weight the slow parts of
            # the flight more heavily than the fast ones.
            if arm == "mean" and u_ref > 0.0:
                sig_p, L_p = _path_field_stats(base, u_ref, r, seg)
                o["path_sigma_u"].append(sig_p)
                o["path_L_int"].append(L_p)
            else:
                o["path_sigma_u"].append(np.nan)
                o["path_L_int"].append(np.nan)
            o["land"].append(r[-1])
            o["T_f"].append(T_f)
            o["L_f"].append(L_f)
            o["V_d"].append(L_f / max(T_f, 1e-9))
            o["apex"].append(float(r[:, 2].max()))
            o["blocked"].append(bool(base.is_solid(r[-1, 0], r[-1, 1])))
    return {a: {k2: np.asarray(v) for k2, v in d.items()}
            for a, d in out.items()}


def field_stats_at(heading, u_ref, launch, field_path, z_probe=None):
    """sigma_u, T_L and L at the launch point, at this wind speed."""
    from wind_scaling import ScaledField
    base = _base_field(field_path, heading)
    f = ScaledField(base, u_ref) if u_ref > 0 else None
    p = np.asarray(launch, dtype=float).copy()
    if z_probe is not None:
        p[2] = z_probe
    if f is None:
        return dict(sigma_u=0.0, T_L=float("nan"), L_int=float("nan"),
                    U_mean=0.0)
    k = float(f.k(p, 0.0))
    eps = max(float(f.eps(p, 0.0)), 1e-12)
    sig = float(np.sqrt(max(2.0 / 3.0 * k, 0.0)))
    T_L = k / eps
    return dict(sigma_u=sig, T_L=T_L, L_int=sig * T_L,
                U_mean=float(np.linalg.norm(f.U(p, 0.0)[:2])))


def throw_frame(bearing):
    """Unit vectors (downrange, crossrange) for a bearing, deg CCW from +x."""
    b = np.radians(bearing)
    return np.array([np.cos(b), np.sin(b)]), np.array([-np.sin(b), np.cos(b)])


def arm_stats(o, launch_xy, e_down, e_cross):
    """Absolute statistics of one arm."""
    L = np.asarray(o["land"])[:, :2] - launch_xy
    down, cross = L @ e_down, L @ e_cross
    mu = np.array([down.mean(), cross.mean()])
    # Radial dispersion about the group's own centroid: the spread a thrower
    # cannot aim off. The bias is reported separately because it is in
    # principle correctable.
    sig = float(np.sqrt(((down - mu[0]) ** 2 + (cross - mu[1]) ** 2).mean()))
    return dict(down_mean=float(mu[0]), cross_mean=float(mu[1]), sigma=sig,
                sigma_down=float(down.std()), sigma_cross=float(cross.std()),
                T_f=float(np.mean(o["T_f"])), L_f=float(np.mean(o["L_f"])),
                V_d=float(np.mean(o["V_d"])), apex=float(np.mean(o["apex"])),
                blocked=float(np.mean(o["blocked"])),
                path_sigma_u=float(np.nanmean(o["path_sigma_u"]))
                if np.isfinite(o["path_sigma_u"]).any() else float("nan"),
                path_L_int=float(np.nanmean(o["path_L_int"]))
                if np.isfinite(o["path_L_int"]).any() else float("nan"))


def paired_effect(o_hi, o_lo, e_down, e_cross):
    """Per-throw displacement between two arms, in the throw frame.

    `bias` is the mean displacement -- the part a thrower who knows the field
    can aim off. `sigma` is the rms of the demeaned displacement -- the part
    they cannot. Returned separately because they have completely different
    consequences and the literature routinely conflates them.
    """
    d = (np.asarray(o_hi["land"])[:, :2] - np.asarray(o_lo["land"])[:, :2])
    down, cross = d @ e_down, d @ e_cross
    mu = np.array([down.mean(), cross.mean()])
    sig = float(np.sqrt(((down - mu[0]) ** 2 + (cross - mu[1]) ** 2).mean()))
    n = len(down)
    return dict(bias_down=float(mu[0]), bias_cross=float(mu[1]),
                bias=float(np.linalg.norm(mu)), sigma=sig,
                sigma_down=float((down - mu[0]).std()),
                sigma_cross=float((cross - mu[1]).std()),
                # Standard error of an rms over n paired samples, to 1/sqrt(2n).
                sigma_se=float(sig / np.sqrt(2.0 * max(n - 1, 1))))


def run(args):
    from concurrent.futures import ProcessPoolExecutor
    import site_view as sv

    launch = sv.LAUNCH
    spins = args.spins if args.spins else [None]
    rspeeds = args.release_speeds if args.release_speeds else [None]
    elevs = args.elevations if args.elevations else [None]
    if args.launches:
        if len(args.launches) % 2:
            raise SystemExit("--launches takes x y pairs")
        launches = [(args.launches[i], args.launches[i + 1], args.launch_z)
                    for i in range(0, len(args.launches), 2)]
    else:
        launches = [tuple(launch)]
    for lp in launches:
        sv.assert_launch_clear(args.field, np.asarray(lp, float), args.headings)
    cells = list(itertools.product(args.headings, args.u_refs, args.throws,
                                   args.bearings, spins, rspeeds, elevs,
                                   launches))
    print(f"{len(cells)} cells x {args.n_throws} throws x 3 arms = "
          f"{len(cells) * args.n_throws * 3} integrations")

    t0 = time.time()
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        from concurrent.futures import as_completed
        futs = {ex.submit(_worker, (h, u, th, br, np.asarray(lp, float),
                                    args.field, args.n_throws, sp, rs, el)):
                (h, u, th, br, sp, rs, el, lp)
                for (h, u, th, br, sp, rs, el, lp) in cells}
        for n, f in enumerate(as_completed(futs), 1):
            rows.append(f.result())
            if n % 10 == 0 or n == len(cells):
                print(f"  {n}/{len(cells)} cells  "
                      f"({time.time() - t0:.0f} s)")
    rows.sort(key=lambda r: (r["launch"], r["heading"], r["throw"],
                             r["bearing"],
                             r["spin"] if r["spin"] is not None else -1.0,
                             r["u_ref"]))
    return rows


def _worker(spec):
    (heading, u_ref, throw, bearing, launch, field, n_throws, spin,
     rspeed, elev) = spec
    cell = fly_cell(heading, u_ref, throw, bearing, launch, field, n_throws,
                    spin=spin, release_speed=rspeed, elevation=elev)
    e_down, e_cross = throw_frame(bearing)
    launch_xy = np.asarray(launch, dtype=float)[:2]
    row = dict(heading=int(heading), u_ref=float(u_ref), throw=throw,
               bearing=float(bearing), n_throws=int(n_throws),
               spin=(float(spin) if spin is not None else None),
               release_speed=(float(rspeed) if rspeed is not None else None),
               elevation=(float(elev) if elev is not None else None),
               launch=[round(float(v), 3) for v in np.asarray(launch, float)])
    for arm, o in cell.items():
        for k2, v in arm_stats(o, launch_xy, e_down, e_cross).items():
            row[f"{arm}__{k2}"] = v
    row.update({f"field__{k2}": v for k2, v in
                field_stats_at(heading, u_ref, launch, field).items()})

    # Paired per-throw decomposition (see the module docstring for why this is
    # not a difference of variances).
    for tag, hi, lo in (("gust", "full", "mean"), ("wind", "mean", "still"),
                        ("total", "full", "still")):
        for k2, v in paired_effect(cell[hi], cell[lo], e_down, e_cross).items():
            row[f"{tag}__{k2}"] = v
    row["sigma_gust"] = row["gust__sigma"]
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--field", default=str(DEFAULT_FIELD))
    ap.add_argument("--u-refs", nargs="+", type=float, default=list(U_REFS))
    ap.add_argument("--headings", nargs="+", type=int, default=list(HEADINGS))
    ap.add_argument("--throws", nargs="+", default=list(THROWS))
    ap.add_argument("--bearings", nargs="+", type=float, default=list(BEARINGS))
    ap.add_argument("--n-throws", type=int, default=N_THROWS)
    ap.add_argument("--launches", nargs="+", type=float, default=None,
                    metavar="X Y",
                    help="release points as flat x y pairs, m in the site "
                         "frame; omit to use site_view.LAUNCH. Every point is "
                         "asserted clear of buildings in every heading before "
                         "any flying happens.")
    ap.add_argument("--launch-z", type=float, default=1.5)
    ap.add_argument("--release-speeds", nargs="+", type=float, default=None,
                    help="override mean release speed, m/s")
    ap.add_argument("--elevations", nargs="+", type=float, default=None,
                    help="override mean release elevation, deg")
    ap.add_argument("--spins", nargs="+", type=float, default=None,
                    help="override mean release spin, rad/s; omit to use the "
                         "catalogue entry's own value")
    ap.add_argument("--workers", type=int,
                    default=max(1, (os.cpu_count() or 4) // 2))
    ap.add_argument("--out", default=str(HERE / "out" / "sweep_wind_speed.json"))
    args = ap.parse_args()

    rows = run(args)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(
        meta=dict(field=args.field, u_refs=args.u_refs, headings=args.headings,
                  throws=args.throws, bearings=args.bearings,
                  spins=args.spins, release_speeds=args.release_speeds,
                  elevations=args.elevations, launches=args.launches,
                  launch_z=args.launch_z,
                  n_throws=args.n_throws, seed=SEED, t_max=T_MAX,
                  gust_model="sweeping"),
        rows=rows), indent=2))
    print(f"wrote {out}  ({len(rows)} cells)")


if __name__ == "__main__":
    main()
