"""
Plot Q(x, y, theta) maps — one contour figure per heading + a summary figure
showing best/worst heading at each launch point.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("Q_npz")
    ap.add_argument("--outdir", default="figs")
    args = ap.parse_args()

    Q_data = np.load(args.Q_npz)
    Q = Q_data["Q"]
    x = Q_data["x"]
    y = Q_data["y"]
    headings = Q_data["headings"]

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Per-heading panel
    n = len(headings)
    cols = min(4, n)
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(3 * cols, 3 * rows),
                             squeeze=False)
    XX, YY = np.meshgrid(x, y, indexing="ij")
    for k, h in enumerate(headings):
        ax = axes[k // cols][k % cols]
        cf = ax.contourf(XX, YY, Q[..., k], levels=np.linspace(0, 1, 11),
                          cmap="viridis", vmin=0, vmax=1)
        ax.set_title(f"$\\theta_{{wind}} = {int(h)}^\\circ$")
        ax.set_aspect("equal")
    for k in range(len(headings), rows * cols):
        axes[k // cols][k % cols].axis("off")
    fig.colorbar(cf, ax=axes.ravel().tolist(), shrink=0.7, label="Q")
    fig.suptitle("Throwability $Q(\\mathbf{x}, \\theta_{wind})$")
    fig.savefig(outdir / "Q_panel_by_heading.png", dpi=150, bbox_inches="tight")
    print(f"wrote {outdir / 'Q_panel_by_heading.png'}")

    # Best heading per launch point
    best_h = headings[np.argmax(Q, axis=-1)]
    best_Q = Q.max(axis=-1)

    fig2, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    cf1 = a1.contourf(XX, YY, best_Q, levels=np.linspace(0, 1, 11),
                       cmap="viridis", vmin=0, vmax=1)
    a1.set_title("Max throwability across all headings")
    a1.set_aspect("equal"); plt.colorbar(cf1, ax=a1, label="$Q_{max}$")

    cf2 = a2.contourf(XX, YY, best_h, cmap="twilight",
                       levels=np.linspace(0, 180, 13))
    a2.set_title("Best wind heading at each launch point")
    a2.set_aspect("equal"); plt.colorbar(cf2, ax=a2, label="$\\theta_{wind}$ (deg)")

    fig2.savefig(outdir / "Q_summary.png", dpi=150, bbox_inches="tight")
    print(f"wrote {outdir / 'Q_summary.png'}")


if __name__ == "__main__":
    main()
