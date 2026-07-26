"""
How much does turbulence cost you, as a function of release height?

Reads the paired sweep written by `release_height_sweep.sh` and separates the
landing scatter into the part a thrower could in principle correct for and the
part they could not.

THE DECOMPOSITION
-----------------
Both arms share a seed, hence identical release-parameter draws, and differ
only in whether the synthetic gust field is applied. Writing the total landing
variance of the gust arm as sigma_g^2 and of the mean-only arm as sigma_m^2:

    sigma_m       release scatter + mean-flow advection + longer flight time
    sigma_excess  = sqrt(max(sigma_g^2 - sigma_m^2, 0))

sigma_excess is the dispersion attributable to turbulence. The subtraction in
quadrature assumes the two contributions are uncorrelated, which holds here
because the gust realization is drawn independently of the release draw.
Where the difference is negative (sampling noise on a small effect) it is
clipped to zero and flagged, rather than silently producing a NaN.

sigma is reported as the total 2-D RMS landing displacement,
sqrt(var_x + var_y) = sqrt(trace(cov)), which is rotation-invariant and so
does not depend on the arbitrary choice of throw bearing.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re

import numpy as np


def sigma_total(npz):
    """Total 2-D RMS landing displacement per cell, from the covariance."""
    cov = npz["landing_cov"]                      # (...,2,2)
    var = cov[..., 0, 0] + cov[..., 1, 1]
    return np.sqrt(np.maximum(var, 0.0))


def load_arm(path):
    d = np.load(path, allow_pickle=False)
    return d, sigma_total(d)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default="height_sweep")
    ap.add_argument("--out-json", default="height_sweep.json")
    ap.add_argument("--figure", default="../P5_Synthesis/figs/"
                                        "fig7_release_height.png")
    args = ap.parse_args()

    heights = sorted({float(re.match(r"h([\d.]+)_", os.path.basename(p)).group(1))
                      for p in glob.glob(f"{args.dir}/h*_gust.npz")})
    if not heights:
        raise SystemExit(f"no sweep output in {args.dir}/")

    print("=" * 78)
    print("RELEASE HEIGHT vs LANDING SCATTER  (paired gust / mean-only)")
    print("=" * 78)
    print(f"{'z_rel':>6} {'sigma_gust':>11} {'sigma_mean':>11} "
          f"{'sigma_excess':>13} {'excess/total':>13} {'success_g':>10}")

    rows = []
    per_heading = {}
    for z in heights:
        gp, mp = f"{args.dir}/h{z}_gust.npz", f"{args.dir}/h{z}_mean.npz"
        if not (os.path.exists(gp) and os.path.exists(mp)):
            print(f"{z:6.1f}   incomplete pair, skipped")
            continue
        g, sg = load_arm(gp)
        m, sm = load_arm(mp)

        # Only compare cells where BOTH arms produced valid landings.
        ok = (g["n_valid"] > 0) & (m["n_valid"] > 0)
        if not ok.any():
            print(f"{z:6.1f}   no valid cells")
            continue
        sg_m = float(np.sqrt(np.mean(sg[ok] ** 2)))
        sm_m = float(np.sqrt(np.mean(sm[ok] ** 2)))
        diff = sg_m ** 2 - sm_m ** 2
        clipped = diff < 0
        sx = float(np.sqrt(max(diff, 0.0)))
        succ = float(np.mean(g["success_rate"][ok]))
        frac = sx / sg_m if sg_m > 0 else float("nan")
        print(f"{z:6.1f} {sg_m:11.2f} {sm_m:11.2f} {sx:13.2f} "
              f"{frac:12.1%} {succ:10.3f}"
              + ("   [excess clipped at 0]" if clipped else ""))
        rows.append(dict(z=z, sigma_gust=sg_m, sigma_mean=sm_m,
                         sigma_excess=sx, excess_fraction=frac,
                         success_gust=succ, clipped=bool(clipped)))

        # Per-heading breakdown: headings axis is index 2.
        hh = [int(v) for v in g["headings"]]
        per_heading[z] = {}
        for i, h in enumerate(hh):
            okh = ok[:, :, i]
            if not okh.any():
                continue
            a = float(np.sqrt(np.mean(sg[:, :, i][okh] ** 2)))
            b = float(np.sqrt(np.mean(sm[:, :, i][okh] ** 2)))
            per_heading[z][h] = dict(
                sigma_gust=a, sigma_mean=b,
                sigma_excess=float(np.sqrt(max(a * a - b * b, 0.0))))

    if not rows:
        raise SystemExit("no complete height pairs to report")

    z = np.array([r["z"] for r in rows])
    sx = np.array([r["sigma_excess"] for r in rows])
    print("\nTurbulence-induced scatter grows by a factor "
          f"{sx[-1]/max(sx[0], 1e-9):.2f} from z = {z[0]:g} m to "
          f"{z[-1]:g} m.")

    with open(args.out_json, "w") as fh:
        json.dump(dict(rows=rows, per_heading=per_heading), fh, indent=2)
    print(f"wrote {args.out_json}")

    # ---------------- figure ------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))
    ax = axes[0]
    ax.plot(z, [r["sigma_gust"] for r in rows], "-o", color="#c1272d",
            label=r"with turbulence  $\sigma_g$")
    ax.plot(z, [r["sigma_mean"] for r in rows], "-s", color="#1f4e9c",
            label=r"mean flow only  $\sigma_m$")
    ax.plot(z, sx, "-^", color="#0b7a3b",
            label=r"turbulence-induced  $\sigma_{\rm excess}$")
    ax.set_xlabel("release height  [m]")
    ax.set_ylabel("landing scatter  [m]")
    ax.set_title("Scatter vs release height", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    for h, c in ((90, "#c1272d"), (210, "#0b7a3b"), (270, "#1f4e9c")):
        zz = [k for k in sorted(per_heading) if h in per_heading[k]]
        if not zz:
            continue
        ax.plot(zz, [per_heading[k][h]["sigma_excess"] for k in zz],
                "-o", color=c, label=rf"$\theta_{{\rm wind}}={h}^\circ$")
    ax.set_xlabel("release height  [m]")
    ax.set_ylabel(r"turbulence-induced scatter  $\sigma_{\rm excess}$  [m]")
    ax.set_title("By wind direction", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    fig.suptitle("Cost of throwing higher: the uncorrectable part of the "
                 "landing error", fontsize=11)
    fig.tight_layout()
    fig.savefig(args.figure, dpi=160)
    print(f"wrote {args.figure}")


if __name__ == "__main__":
    main()
