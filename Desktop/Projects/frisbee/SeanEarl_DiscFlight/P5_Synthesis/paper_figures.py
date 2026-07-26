"""
Build the four paper figures the manuscript draft references.
This is a layout stub — final styling and panel ordering should be tuned
after the actual data is available.
"""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt


def fig1_site_map(outdir: Path):
    """Site overview: Sean Earl Field, buildings, lake shoreline, IBL schematic."""
    fig, ax = plt.subplots(figsize=(6, 5))
    # Lake shoreline (east, schematic)
    ax.fill_between([60, 250], -250, 250, alpha=0.3, color="steelblue",
                    label="Lake Michigan")
    # Sean Earl Field rectangle
    ax.add_patch(plt.Rectangle((-30, -45), 60, 90, fill=False, edgecolor="green",
                                linewidth=2, label="Sean Earl Field"))
    # Mertz (W of field)
    ax.add_patch(plt.Rectangle((-80, -30), 30, 60, fill=True, color="gray",
                                alpha=0.6, label="Mertz Hall"))
    # Norville Center (S of field)
    ax.add_patch(plt.Rectangle((-50, -90), 100, 30, fill=True, color="gray",
                                alpha=0.6))
    # IBL schematic
    xs = np.linspace(60, 250, 50)
    delta = 0.10 * (xs - 60) ** 0.8
    ax.plot(xs, delta + 30, "k--", alpha=0.7, label="$\\delta_{IBL}(x)$")
    ax.set_xlim(-150, 260); ax.set_ylim(-150, 150)
    ax.set_aspect("equal")
    ax.set_xlabel("Local x (East), m"); ax.set_ylabel("Local y (North), m")
    ax.legend(loc="upper left", fontsize=8)
    ax.set_title("Sean Earl Field — site schematic")
    fig.savefig(outdir / "fig1_site_map.png", dpi=200, bbox_inches="tight")


def fig2_trajectories(outdir: Path):
    """3-DOF trim glide example trajectories under (still, head, tail) wind."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "P1_AnalyticalModel"))
    import disc_6dof as d6
    from trajectories import right_handed_backhand

    fig, ax = plt.subplots(figsize=(7, 4))
    for label, U in [("still", [0, 0, 0]),
                     ("headwind 5 m/s", [-5, 0, 0]),
                     ("tailwind 5 m/s", [5, 0, 0])]:
        r0, v0, q0, w0 = right_handed_backhand(release_alpha_deg=8.0, roll_deg=0.0)
        res = d6.simulate(r0, v0, q0, w0, wind=d6.uniform_wind(U), t_max=8.0)
        ax.plot(res.r[:, 0], res.r[:, 2], label=label)
    ax.set_xlabel("x (m)"); ax.set_ylabel("z (m)")
    ax.set_title("6-DOF trajectories under uniform wind (Hummel coeffs)")
    ax.legend()
    ax.set_aspect("equal")
    fig.savefig(outdir / "fig2_trajectories.png", dpi=200, bbox_inches="tight")


def fig3_cfd_snapshots(outdir: Path, fields_npz: str | None):
    """CFD mean-flow snapshots at 2 m AGL for 3 representative headings."""
    if fields_npz is None or not Path(fields_npz).exists():
        print("  [fig3] no CFD field provided; skipping")
        return
    F = np.load(fields_npz)
    headings_to_show = [0, 90, 180]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    zk = int(np.argmin(np.abs(F["z"] - 2.0)))
    XX, YY = np.meshgrid(F["x"], F["y"], indexing="ij")
    for ax, h in zip(axes, headings_to_show):
        key = f"dir_{int(h):03d}"
        if f"{key}__U" not in F.files:
            ax.text(0.5, 0.5, f"no data\n{h}°", transform=ax.transAxes,
                    ha="center"); ax.axis("off"); continue
        U = F[f"{key}__U"][:, :, zk, :]
        mag = np.linalg.norm(U, axis=-1)
        cf = ax.contourf(XX, YY, mag, levels=20, cmap="viridis")
        ax.streamplot(F["x"], F["y"], U[..., 0].T, U[..., 1].T,
                      color="white", density=0.7, linewidth=0.6)
        ax.set_title(f"$\\theta_{{wind}} = {h}^\\circ$, |U| at 2 m AGL")
        ax.set_aspect("equal")
        plt.colorbar(cf, ax=ax, label="|U| (m/s)")
    fig.savefig(outdir / "fig3_cfd_snapshots.png", dpi=200, bbox_inches="tight")


def fig4_throwability(outdir: Path, Q_npz: str | None):
    """Throwability heatmap at best vs worst heading."""
    if Q_npz is None or not Path(Q_npz).exists():
        print("  [fig4] no Q.npz provided; skipping")
        return
    d = np.load(Q_npz)
    Q = d["Q"]; x = d["x"]; y = d["y"]; h = d["headings"]
    best_h = h[np.argmax(Q.mean(axis=(0, 1)))]
    worst_h = h[np.argmin(Q.mean(axis=(0, 1)))]
    XX, YY = np.meshgrid(x, y, indexing="ij")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    for ax, hh in zip((a1, a2), (best_h, worst_h)):
        ih = int(np.argmin(np.abs(h - hh)))
        cf = ax.contourf(XX, YY, Q[..., ih], levels=np.linspace(0, 1, 11),
                          cmap="viridis", vmin=0, vmax=1)
        ax.set_title(f"$\\theta_{{wind}} = {int(hh)}^\\circ$")
        ax.set_aspect("equal")
        plt.colorbar(cf, ax=ax, label="Q")
    fig.suptitle("Throwability metric at best (left) vs. worst (right) wind heading")
    fig.savefig(outdir / "fig4_throwability.png", dpi=200, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Q", default=None)
    ap.add_argument("--fields", default=None)
    ap.add_argument("--outdir", default="figs")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    fig1_site_map(outdir)
    fig2_trajectories(outdir)
    fig3_cfd_snapshots(outdir, args.fields)
    fig4_throwability(outdir, args.Q)
    print(f"figures saved under {outdir}/")


if __name__ == "__main__":
    main()
