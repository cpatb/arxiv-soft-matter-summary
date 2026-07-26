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
    print("Mean drift across all launch points, per heading, per throw:")
    drift = d["landing_mean"]  # (nx, ny, nh, nt, nb, 3)
    for ih, h in enumerate(d["headings"]):
        for it, throw in enumerate(d["throws"]):
            ave = drift[:, :, ih, it, :, :2].mean(axis=(0, 1, 2))
            print(f"  heading {h}°  throw {str(throw):20s}  "
                  f"mean landing dx={ave[0]:.2f} dy={ave[1]:.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    args = ap.parse_args()
    summarize(args.npz)


if __name__ == "__main__":
    main()
