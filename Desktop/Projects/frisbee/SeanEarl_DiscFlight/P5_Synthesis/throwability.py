"""
Compute the throwability metric Q(x, y, theta_wind) from P3's CFD field and
P4's Monte Carlo landing statistics.

METRIC
------
Following the manuscript,

    Q(x, theta) = (1 - |U_2m(x)| / U_ref) * (1 - I_2m(x)) * S_bar(x, theta)

with every factor clipped to [0, 1]:

  * wind-magnitude penalty, from the RANS mean speed at 2 m AGL;
  * turbulence penalty, from I = sqrt(2k/3) / |U| at the same height;
  * simulated success rate, the fraction of Monte Carlo throws delivered
    within 3 m of a 25 m target, averaged over throw types and bearings.

Q = 1 is calm-air-equivalent, Q = 0 unplayable.

DEVIATION FROM THE DRAFT EQUATION. The draft writes the first factor with a
field-averaged speed, <|U_2m|>_field, which is constant in x and therefore
cannot contribute any spatial structure. Since the paper's central claim is a
spatial Q gradient driven by the local flow, we use the LOCAL speed |U_2m(x)|.
This is recorded in the manuscript's methods; `--field-average-wind` restores
the literal draft form for comparison.

The success-rate factor is the MEASURED fraction of Monte Carlo throws that hit
the target, computed in P4 from the actual landing points. An earlier version of
this module estimated it from the 90% ellipse axes as
min(1, (r_target/a_mean)^2), which is not a probability -- it ignores the offset
between the landing centroid and the target, so a tightly-grouped but badly
aimed throw scored as a success. The measured rate is used instead.
"""

from __future__ import annotations

import argparse

import numpy as np


def wind_factor(U_2m, U_ref):
    """Wind-magnitude penalty, clipped to [0, 1]."""
    mag = np.linalg.norm(U_2m, axis=-1)
    return np.clip(1.0 - mag / U_ref, 0.0, 1.0)


def turbulence_intensity(k_2m, U_2m, I_cap=1.0):
    """I = sqrt(2k/3) / |U|, with a defined value where |U| -> 0."""
    sigma = np.sqrt(2.0 * np.maximum(k_2m, 0.0) / 3.0)
    mag = np.linalg.norm(U_2m, axis=-1)
    return np.where(mag > 1e-3, sigma / np.maximum(mag, 1e-3), I_cap)


def turbulence_factor(k_2m, U_2m, I_cap=1.0):
    """Turbulence-intensity penalty (1 - I), clipped to [0, 1].

    SATURATION WARNING. The draft defines this factor as (1 - I), which goes
    negative whenever I > 1 -- i.e. wherever the RMS fluctuation exceeds the
    local mean speed. That is not exotic: it is the normal state of a
    recirculating building wake, exactly the sheltered region the paper is
    interested in. Clipping at zero is the only sane choice (a calm sheltered
    spot must not score worse than a gale), but it means Q collapses to
    identically zero across the whole saturated region and carries no spatial
    information there.

    `compute_Q` reports the saturated fraction per heading so this can be
    stated rather than silently absorbed. If a large fraction saturates, the
    honest reading is that (1 - I) is the wrong functional form for sheltered
    sites and a bounded penalty such as 1/(1 + I) should be used instead; that
    is flagged as a methodological finding rather than changed unilaterally,
    since the (1 - I) form is what the manuscript defines.
    """
    I = turbulence_intensity(k_2m, U_2m, I_cap)
    return np.clip(1.0 - I, 0.0, 1.0)


def compute_Q(fields_npz, landings_npz, U_ref, field_average_wind=False,
              allow_stub=False):
    """Return (Q, launch_x, launch_y, headings, parts).

    `parts` carries the three factors separately so the paper can report which
    term drives the spatial structure.
    """
    from scipy.interpolate import RegularGridInterpolator

    F = np.load(fields_npz, allow_pickle=False)
    L = np.load(landings_npz, allow_pickle=False)

    if bool(F["is_stub"]) if "is_stub" in F.files else False:
        if not allow_stub:
            raise SystemExit(
                f"{fields_npz} is a STUB field generated for pipeline testing, "
                "not CFD output. Refusing to produce results from it. Pass "
                "--allow-stub if you are deliberately smoke-testing."
            )
        print("  ** WARNING: using a STUB wind field. Not physical. **")

    if "success_rate" not in L.files:
        raise SystemExit(
            f"{landings_npz} has no 'success_rate' entry. Re-run "
            "P4/monte_carlo.py, which now measures it directly from the "
            "landing points."
        )

    headings = L["headings"]
    launch_x, launch_y = L["launch_x"], L["launch_y"]
    nx, ny, nh = len(launch_x), len(launch_y), len(headings)

    Q = np.full((nx, ny, nh), np.nan)
    wf_all = np.full((nx, ny, nh), np.nan)
    tf_all = np.full((nx, ny, nh), np.nan)
    sr_all = np.full((nx, ny, nh), np.nan)
    I_all = np.full((nx, ny, nh), np.nan)

    zk = int(np.argmin(np.abs(F["z"] - 2.0)))
    if not np.isclose(F["z"][zk], 2.0, atol=0.51):
        print(f"  note: nearest sampled altitude to 2 m is {F['z'][zk]:.2f} m")

    XX, YY = np.meshgrid(launch_x, launch_y, indexing="ij")
    pts = np.column_stack([XX.ravel(), YY.ravel()])

    for ih, h in enumerate(headings):
        key = f"dir_{int(h):03d}"
        if f"{key}__U" not in F.files:
            print(f"  heading {h}: no CFD field, skipping")
            continue

        U_2m = F[f"{key}__U"][:, :, zk, :]
        k_2m = F[f"{key}__k"][:, :, zk]

        Ui = [RegularGridInterpolator((F["x"], F["y"]), U_2m[..., i],
                                      bounds_error=False, fill_value=None)
              for i in range(3)]
        ki = RegularGridInterpolator((F["x"], F["y"]), k_2m,
                                     bounds_error=False, fill_value=None)

        Ulaunch = np.stack([ui(pts) for ui in Ui], axis=-1).reshape(nx, ny, 3)
        klaunch = ki(pts).reshape(nx, ny)

        if field_average_wind:
            mean_speed = np.linalg.norm(U_2m, axis=-1).mean()
            wf = np.full((nx, ny), np.clip(1.0 - mean_speed / U_ref, 0.0, 1.0))
        else:
            wf = wind_factor(Ulaunch, U_ref)
        tf = turbulence_factor(klaunch, Ulaunch)

        # Measured success rate, averaged over throw types and bearings.
        # (nx, ny, nt, nb) -> (nx, ny)
        sr = np.nanmean(L["success_rate"][:, :, ih, ...], axis=(-1, -2))

        Q[:, :, ih] = wf * tf * sr
        wf_all[:, :, ih], tf_all[:, :, ih], sr_all[:, :, ih] = wf, tf, sr
        I_all[:, :, ih] = turbulence_intensity(klaunch, Ulaunch)

    sat = np.nanmean(I_all >= 1.0, axis=(0, 1))
    for ih, h in enumerate(headings):
        if sat[ih] > 0.05:
            print(f"  ** heading {int(h)}: turbulence factor saturates at 0 "
                  f"over {100*sat[ih]:.0f}% of the grid (I >= 1); Q carries no "
                  f"spatial information there **")

    parts = dict(wind_factor=wf_all, turbulence_factor=tf_all,
                 success_rate=sr_all, turbulence_intensity=I_all,
                 saturated_fraction=sat)
    return Q, launch_x, launch_y, headings, parts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fields", required=True)
    ap.add_argument("--landings", required=True)
    ap.add_argument("--u-ref", type=float, default=5.0)
    ap.add_argument("--field-average-wind", action="store_true",
                    help="use the draft's spatially-averaged wind penalty")
    ap.add_argument("--allow-stub", action="store_true",
                    help="permit a stub (non-CFD) wind field")
    ap.add_argument("--out", default="Q.npz")
    args = ap.parse_args()

    Q, x, y, h, parts = compute_Q(
        args.fields, args.landings, args.u_ref,
        field_average_wind=args.field_average_wind,
        allow_stub=args.allow_stub)

    np.savez_compressed(args.out, Q=Q, x=x, y=y, headings=h,
                        u_ref=args.u_ref, **parts)
    print(f"wrote {args.out}: Q shape {Q.shape}, "
          f"range [{np.nanmin(Q):.4f}, {np.nanmax(Q):.4f}]")
    print(f"  {'heading':>8s} {'Q_mean':>8s} {'Q_min':>8s} {'Q_max':>8s} "
          f"{'wind':>7s} {'turb':>7s} {'success':>8s}")
    for ih, hh in enumerate(h):
        print(f"  {int(hh):>7d}d {np.nanmean(Q[..., ih]):>8.4f} "
              f"{np.nanmin(Q[..., ih]):>8.4f} {np.nanmax(Q[..., ih]):>8.4f} "
              f"{np.nanmean(parts['wind_factor'][..., ih]):>7.3f} "
              f"{np.nanmean(parts['turbulence_factor'][..., ih]):>7.3f} "
              f"{np.nanmean(parts['success_rate'][..., ih]):>8.3f}")


if __name__ == "__main__":
    main()
