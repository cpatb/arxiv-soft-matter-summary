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


def fig1_site_map(outdir: Path, geometry_json=None, half=260.0):
    """Site overview from the REAL OSM footprints, height-shaded.

    Falls back to the original schematic only if the geometry file is missing,
    since a schematic that contradicts the built model would be worse than no
    figure at all.
    """
    gj = Path(geometry_json) if geometry_json else Path(
        __file__).resolve().parent.parent / "P2_Geometry" / "site_geometry.json"
    if not gj.exists():
        print(f"  [fig1] {gj} missing; falling back to schematic")
        return _fig1_schematic(outdir)

    import json
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection

    G = json.loads(gj.read_text())
    field = np.asarray(G["field"]["ring"], dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 7))
    patches, heights = [], []
    for b in G["buildings"]:
        ring = np.asarray(b["ring"])
        if abs(ring[:, 0].mean()) > half or abs(ring[:, 1].mean()) > half:
            continue
        patches.append(MplPolygon(ring, closed=True))
        heights.append(b["height_m"])

    pc = PatchCollection(patches, cmap="inferno", edgecolor="k",
                         linewidth=0.3, alpha=0.92)
    pc.set_array(np.array(heights))
    pc.set_clim(0, max(heights))
    ax.add_collection(pc)
    fig.colorbar(pc, ax=ax, shrink=0.75, label="building height (m)")

    ax.add_patch(MplPolygon(field, closed=True, facecolor="limegreen",
                            edgecolor="darkgreen", lw=2.0, alpha=0.7, zorder=5))
    ax.annotate("Sean Earl Field", (field[:, 0].mean(), field[:, 1].mean()),
                fontsize=9, ha="center", va="center", zorder=8, weight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", alpha=0.85,
                          lw=0))

    # Label the two buildings that dominate the near-field flow.
    for key, note in (("Mertz Hall", "61 m"), ("Campus Towers", "45 m")):
        hit = [b for b in G["buildings"] if b.get("name") == key]
        if not hit:
            continue
        r = np.asarray(hit[0]["ring"])
        ax.annotate(f"{key}\n{note}", (r[:, 0].mean(), r[:, 1].mean()),
                    fontsize=8, ha="center", va="center", zorder=7,
                    color="white",
                    bbox=dict(boxstyle="round,pad=0.25", fc="black",
                              alpha=0.7, lw=0))

    ax.annotate("Lake Michigan\n(east)", (half * 0.8, 0), fontsize=9,
                ha="center", color="steelblue", weight="bold")
    ax.set_xlim(-half, half); ax.set_ylim(-half, half)
    ax.set_aspect("equal")
    ax.set_xlabel("x (East), m"); ax.set_ylabel("y (North), m")
    ax.set_title("Sean Earl Field: OSM building footprints, height-shaded\n"
                 f"origin {G['origin_lat']:.5f}, {G['origin_lon']:.5f}")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(outdir / "fig1_site_map.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def _fig1_schematic(outdir: Path):
    """Original hand-drawn schematic; retained as a fallback only."""
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
    ax.set_title("Sean Earl Field — site schematic (FALLBACK)")
    fig.savefig(outdir / "fig1_site_map.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


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
    plt.close(fig)


def fig3_cfd_snapshots(outdir: Path, fields_npz: str | None):
    """CFD mean-flow snapshots at 2 m AGL for 3 representative headings."""
    if fields_npz is None or not Path(fields_npz).exists():
        print("  [fig3] no CFD field provided; skipping")
        return
    F = np.load(fields_npz)
    # Show whichever directions were actually computed, not a fixed triple.
    if "headings" in F.files:
        headings_to_show = [int(h) for h in F["headings"]]
    else:
        headings_to_show = sorted(int(n[4:7]) for n in F.files
                                  if n.endswith("__U"))
    n = max(len(headings_to_show), 1)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4.4), squeeze=False)
    axes = axes[0]
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
                      color="white", density=0.8, linewidth=0.6)
        # Outline the playing surface so the reader can see where the field
        # sits relative to the wake structure.
        ax.add_patch(plt.Rectangle((-31.4, -10.0), 78.3, 25.1, fill=False,
                                   edgecolor="limegreen", lw=2.0, zorder=6))
        ax.set_title(f"$\\theta_{{wind}} = {h}^\\circ$, |U| at 2 m AGL")
        ax.set_aspect("equal")
        ax.set_xlabel("x (East), m")
        ax.set_ylabel("y (North), m")
        plt.colorbar(cf, ax=ax, label="|U| (m/s)")
    fig.suptitle("CFD mean flow at 2 m AGL; green outline = Sean Earl Field")
    fig.tight_layout()
    fig.savefig(outdir / "fig3_cfd_snapshots.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def fig4_throwability(outdir: Path, Q_npz: str | None):
    """Throwability heatmap at best vs worst heading."""
    if Q_npz is None or not Path(Q_npz).exists():
        print("  [fig4] no Q.npz provided; skipping")
        return
    d = np.load(Q_npz)
    Q = d["Q"]; x = d["x"]; y = d["y"]; h = d["headings"]
    mean_by_heading = np.nanmean(Q, axis=(0, 1))
    best_h = h[np.nanargmax(mean_by_heading)]
    worst_h = h[np.nanargmin(mean_by_heading)]
    XX, YY = np.meshgrid(x, y, indexing="ij")

    # Scale the colour range to the data. Q is a product of three penalties, so
    # realistic values sit well below 1; a fixed 0-1 scale renders every map as
    # a uniform dark patch and hides the spatial gradient that is the result.
    vmax = float(np.nanmax(Q))
    vmax = 1.0 if not np.isfinite(vmax) or vmax <= 0 else vmax

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    for ax, hh, lbl in ((a1, best_h, "best"), (a2, worst_h, "worst")):
        ih = int(np.argmin(np.abs(h - hh)))
        cf = ax.contourf(XX, YY, Q[..., ih], levels=np.linspace(0, vmax, 11),
                         cmap="viridis", vmin=0, vmax=vmax)
        ax.set_title(f"{lbl}: $\\theta_{{wind}} = {int(hh)}^\\circ$ "
                     f"($\\bar{{Q}}$ = {mean_by_heading[ih]:.3f})")
        ax.set_aspect("equal")
        ax.set_xlabel("x (East), m")
        ax.set_ylabel("y (North), m")
        plt.colorbar(cf, ax=ax, label="Q")
    fig.suptitle("Throwability metric at best (left) vs. worst (right) wind heading")
    fig.tight_layout()
    fig.savefig(outdir / "fig4_throwability.png", dpi=200, bbox_inches="tight")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--Q", default=None)
    ap.add_argument("--fields", default=None)
    ap.add_argument("--outdir", default="figs")
    ap.add_argument("--geometry", default=None,
                    help="P2 site_geometry.json for the real site map")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    fig1_site_map(outdir, args.geometry)
    fig2_trajectories(outdir)
    fig3_cfd_snapshots(outdir, args.fields)
    fig4_throwability(outdir, args.Q)
    print(f"figures saved under {outdir}/")


if __name__ == "__main__":
    main()
