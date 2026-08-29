"""
Convenience reader / aggregator for `landing_stats.npz`.

Computes derived quantities:
    - success rate at canonical 25 m target
    - directional curl (turn + fade) per throw
    - mean drift vs heading
"""

from __future__ import annotations
import argparse
import numpy as np


def success_rate_at_target(landings, target_xy, radius_m=3.0):
    dx = landings[:, 0] - target_xy[0]
    dy = landings[:, 1] - target_xy[1]
    hit = (dx ** 2 + dy ** 2) < radius_m ** 2
    return float(hit.mean())


def summarize(npz_path: str):
    d = np.load(npz_path)
    print("launch grid:",
          d["launch_x"].shape, "x", d["launch_y"].shape)
    print("headings:", d["headings"])
    print("throws:", d["throws"])
    print("bearings:", d["bearings"])

    print()
    print("Mean landing bias, throw frame, averaged over launch points and")
    print("bearings (downrange = along the throw, cross = thrower's left):")
    # MUST use landing_bias, not landing_mean. `landing_mean` is an ABSOLUTE
    # world-frame position, so averaging it over a 65 m launch grid and six
    # bearings returns the centroid of the launch grid plus a rotation-averaged
    # range vector -- a number with no interpretation at all, which is what
    # this function used to print under the label "mean drift".
    if "landing_bias" not in d.files:
        print("  landing_bias not present: this artifact predates the frame")
        print("  fix (see README 'Output schema'). Re-run monte_carlo.py.")
        return
    bias = d["landing_bias"]                       # (nx, ny, nh, nt, nb, 2)
    for ih, h in enumerate(d["headings"]):
        for it, throw in enumerate(d["throws"]):
            b = bias[:, :, ih, it, :, :]
            ave = np.nanmean(b, axis=(0, 1, 2))
            sd = np.nanstd(b, axis=(0, 1, 2))
            print(f"  heading {h}°  throw {str(throw):20s}  "
                  f"downrange {ave[0]:+.2f} +- {sd[0]:.2f} m   "
                  f"cross {ave[1]:+.2f} +- {sd[1]:.2f} m")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    args = ap.parse_args()
    summarize(args.npz)


if __name__ == "__main__":
    main()
