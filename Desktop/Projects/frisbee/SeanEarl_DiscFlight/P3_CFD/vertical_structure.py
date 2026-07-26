"""
Vertical structure of the disc-altitude flow over Sean Earl Field.

WHY THIS EXISTS
---------------
The RANS solution is fully three-dimensional, but the manuscript sampled it
only on six horizontal planes at z = 1..6 m and reported the 2 m plane. The
vertical structure -- which is what decides how high a disc can safely be
thrown -- was computed and then discarded at the sampling step. This script
recovers it from `vertical_fields.npz` (the same solved fields, re-sampled to
z = 0.5..100 m).

WHAT IT REPORTS, AND WHY EACH QUANTITY MATTERS TO A DISC
-------------------------------------------------------
  |U|(z)      mean speed. Sets the steady aerodynamic trim; a disc can be
              aimed off to compensate, so mean wind is the *correctable* part.
  TI(z)       turbulence intensity, sqrt(2k/3)/|U|. The UNcorrectable part:
              the thrower cannot aim off a fluctuation they cannot see.
  sigma_u(z)  sqrt(2k/3), the velocity fluctuation in absolute terms. TI alone
              is misleading in deep shelter, where |U| -> 0 drives TI -> inf
              while the absolute perturbation is small. Both are reported.
  |w|(z)      mean vertical velocity. Updraught/downdraught off a building
              face adds directly to the disc's angle of attack, so it changes
              lift without changing heading -- it makes throws land short or
              long rather than left or right.
  dU/dz(z)    vertical shear. A disc is a gyroscope: a rolling moment from
              shear precesses into a HEADING change 90 deg later. Shear is
              therefore the disc-specific hazard that a ballistic point-mass
              model would miss entirely.

METHOD NOTES / ASSUMPTIONS
--------------------------
* Statistics are taken over the playing surface only (the same
  x in [-31, 47], y in [-10, 15] box the manuscript's macros use), and are
  masked with the 3-D fluid mask `__valid3d`. Points inside buildings were
  nearest-neighbour filled by the sampler; including them would read that fill
  as if it were flow. Masked cells are excluded, not zero-filled.
* k is the MODELLED turbulent kinetic energy of a steady k-omega SST solution.
  It is the variance of the unresolved fluctuations, not a resolved gust
  field. Anything here about gust TIMING is therefore a statement about the
  statistics only; see the LES caveat in the manuscript.
* dU/dz is a second-order central difference on a non-uniform z grid via
  np.gradient, which handles the varying spacing correctly.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

# Playing surface, matching P5_Synthesis/make_results_tex.py exactly.
PITCH_X = (-31.0, 47.0)
PITCH_Y = (-10.0, 15.0)

# Band a thrown disc actually occupies: release ~1.5 m, apex of a normal
# backhand ~4-6 m. Used to report band-averaged quantities.
DISC_BAND = (1.0, 6.0)

# Canonical 25 m backhand, used only to form dimensionless ratios against the
# turbulence scales. Flight time from the P1 6-DOF model's still-air solution.
T_FLIGHT = 2.5     # s
L_FLIGHT = 25.0    # m


def pitch_mask(x, y):
    fx = (x >= PITCH_X[0]) & (x <= PITCH_X[1])
    fy = (y >= PITCH_Y[0]) & (y <= PITCH_Y[1])
    return fx, fy


def profiles(F, theta):
    """Return dict of pitch-averaged vertical profiles for one heading.

    Averages are over fluid cells of the playing surface at each height. A
    height whose pitch cells are entirely inside buildings yields NaN rather
    than a fabricated value.
    """
    x, y, z = F["x"], F["y"], F["z"]
    fx, fy = pitch_mask(x, y)
    U = F[f"dir_{theta:03d}__U"][np.ix_(fx, fy)]          # (nx,ny,nz,3)
    k = F[f"dir_{theta:03d}__k"][np.ix_(fx, fy)]          # (nx,ny,nz)
    eps = F[f"dir_{theta:03d}__eps"][np.ix_(fx, fy)]      # (nx,ny,nz)
    ok = F[f"dir_{theta:03d}__valid3d"][np.ix_(fx, fy)]   # (nx,ny,nz)

    mag = np.linalg.norm(U, axis=-1)
    w = U[..., 2]
    sigma = np.sqrt(2.0 * np.maximum(k, 0.0) / 3.0)

    # Turbulence integral scales. T_L = k/eps is the eddy turnover (Lagrangian
    # integral) time; L = sigma_u * T_L the corresponding length. These set the
    # regime: comparing T_L with the ~2.5 s flight time of a 25 m throw says
    # whether the disc flies through a frozen field or a fluctuating one, and
    # comparing L with the 25 m flight length says how many independent eddies
    # it samples -- i.e. how much the disturbance averages out along the path.
    T_L = k / np.maximum(eps, 1e-12)
    L_int = sigma * T_L

    # Mask building interiors out of every statistic.
    m = np.where(ok, 1.0, np.nan)
    magm, wm, sigm = mag * m, w * m, sigma * m

    with np.errstate(invalid="ignore"):
        Uz = np.nanmean(magm, axis=(0, 1))
        sz = np.nanmean(sigm, axis=(0, 1))
        wz = np.nanmean(np.abs(wm), axis=(0, 1))
        # TI from the ratio of the means, not the mean of the ratio: the
        # pointwise ratio blows up wherever |U| -> 0 in the wake core and the
        # average would then be dominated by near-stagnant cells.
        TIz = sz / np.maximum(Uz, 1e-6)
        # Shear of the pitch-mean profile.
        dUdz = np.gradient(Uz, z)
        TLz = np.nanmean(T_L * m, axis=(0, 1))
        Lz = np.nanmean(L_int * m, axis=(0, 1))
        nfluid = np.nansum(ok, axis=(0, 1))

    return dict(z=z, U=Uz, sigma=sz, TI=TIz, w_abs=wz, dUdz=dUdz,
                T_L=TLz, L_int=Lz, n_fluid_cells=nfluid)


def band_average(z, v, lo, hi):
    """Thickness-weighted mean of v(z) over [lo, hi] by trapezoid."""
    sel = (z >= lo) & (z <= hi)
    if sel.sum() < 2:
        return float("nan")
    return float(np.trapezoid(v[sel], z[sel]) / (z[sel][-1] - z[sel][0]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fields", default="vertical_fields.npz")
    ap.add_argument("--control", default="abl_control.npz",
                    help="empty-fetch case sampled on the same z stack; the "
                         "no-buildings counterfactual")
    ap.add_argument("--out-json", default="vertical_structure.json")
    ap.add_argument("--outdir", default="../P5_Synthesis/figs")
    ap.add_argument("--u-ref", type=float, default=5.0)
    args = ap.parse_args()

    F = np.load(args.fields, allow_pickle=False)
    headings = [int(h) for h in F["headings"]]
    z = F["z"]

    print("=" * 74)
    print("VERTICAL STRUCTURE over the playing surface")
    print(f"  pitch x {PITCH_X} m, y {PITCH_Y} m; reference {args.u_ref} m/s "
          "at 10 m")
    print("=" * 74)

    out = {}
    P = {}
    for h in headings:
        p = P[h] = profiles(F, h)
        print(f"\n--- theta_wind = {h} deg " + "-" * 45)
        print(f"{'z [m]':>7} {'|U|':>7} {'sigma_u':>8} {'TI':>7} "
              f"{'|w|':>7} {'dU/dz':>8}  {'fluid cells':>11}")
        for i, zz in enumerate(z):
            print(f"{zz:7.1f} {p['U'][i]:7.2f} {p['sigma'][i]:8.2f} "
                  f"{p['TI'][i]:7.2f} {p['w_abs'][i]:7.3f} "
                  f"{p['dUdz'][i]:8.3f}  {int(p['n_fluid_cells'][i]):11d}")

        band = {q: band_average(z, p[q], *DISC_BAND)
                for q in ("U", "sigma", "TI", "w_abs", "dUdz",
                          "T_L", "L_int")}
        out[str(h)] = dict(
            z=z.tolist(),
            **{q: np.asarray(p[q]).tolist()
               for q in ("U", "sigma", "TI", "w_abs", "dUdz",
                         "T_L", "L_int")},
            disc_band=band,
        )
        print(f"  disc band {DISC_BAND[0]}-{DISC_BAND[1]} m: "
              f"|U|={band['U']:.2f} m/s, sigma_u={band['sigma']:.2f} m/s, "
              f"TI={band['TI']:.2f}, dU/dz={band['dUdz']:.3f} 1/s")
        print(f"    integral scales: T_L={band['T_L']:.1f} s "
              f"(flight {T_FLIGHT} s -> T_L/T_f={band['T_L']/T_FLIGHT:.1f}), "
              f"L={band['L_int']:.1f} m "
              f"(flight {L_FLIGHT} m -> N_eddy={L_FLIGHT/band['L_int']:.1f})")

    # ---- the question a thrower actually asks: does it get worse with height?
    print("\n" + "=" * 74)
    print("HEIGHT DEPENDENCE  (what changes as you throw higher)")
    print("=" * 74)
    for h in headings:
        p = P[h]
        i2 = int(np.argmin(np.abs(z - 2.0)))
        for target in (6.0, 10.0, 20.0):
            it = int(np.argmin(np.abs(z - target)))
            print(f"  {h:3d} deg  z={target:5.1f} m vs 2 m: "
                  f"|U| x{p['U'][it]/p['U'][i2]:4.1f}, "
                  f"sigma_u x{p['sigma'][it]/p['sigma'][i2]:4.1f}, "
                  f"TI {p['TI'][i2]:.2f} -> {p['TI'][it]:.2f}")

    # ---- the counterfactual: what would the same wind do with no buildings?
    # This is the scientifically load-bearing comparison. Mean wind is the
    # CORRECTABLE part of the disturbance (a thrower aims off it); sigma_u is
    # the UNCORRECTABLE part. If shelter removes mean wind but leaves -- or
    # raises -- sigma_u, then the site trades an error the thrower can cancel
    # for one they cannot, and "sheltered" is not the same as "easy".
    try:
        C = np.load(args.control, allow_pickle=False)
    except FileNotFoundError:
        C = None
        print(f"\n!! control {args.control} not found; skipping counterfactual")
    if C is not None:
        ctheta = int(C["headings"][0])
        cp = profiles(C, ctheta)
        out["control"] = {q: np.asarray(cp[q]).tolist()
                          for q in ("z", "U", "sigma", "TI", "w_abs", "dUdz",
                                    "T_L", "L_int")}
        print("\n" + "=" * 74)
        print("COUNTERFACTUAL: sheltered pitch vs the SAME ABL over empty "
              "fetch")
        print("=" * 74)
        print("  correctable = mean wind (aim off it); "
              "uncorrectable = sigma_u (cannot)")
        i2 = int(np.argmin(np.abs(z - 2.0)))
        print(f"\n  empty fetch at 2 m: |U|={cp['U'][i2]:.2f} m/s, "
              f"sigma_u={cp['sigma'][i2]:.2f} m/s, TI={cp['TI'][i2]:.2f}")
        print(f"\n{'theta':>6} {'|U| ratio':>10} {'sigma ratio':>12} "
              f"{'TI pitch':>9} {'TI empty':>9}")
        ratios = {}
        for h in headings:
            p = P[h]
            ru = p["U"][i2] / cp["U"][i2]
            rs = p["sigma"][i2] / cp["sigma"][i2]
            ratios[str(h)] = dict(U_ratio=float(ru), sigma_ratio=float(rs),
                                  TI_pitch=float(p["TI"][i2]),
                                  TI_empty=float(cp["TI"][i2]))
            print(f"{h:6d} {ru:10.2f} {rs:12.2f} {p['TI'][i2]:9.2f} "
                  f"{cp['TI'][i2]:9.2f}")
        out["counterfactual_2m"] = ratios
        print("\n  Reading: a ratio < 1 means the buildings REMOVED that "
              "component;\n  > 1 means they ADDED it.")

    with open(args.out_json, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {args.out_json}")


if __name__ == "__main__":
    main()
