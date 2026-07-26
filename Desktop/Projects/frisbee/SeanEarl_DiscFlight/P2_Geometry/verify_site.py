"""
Sanity-check the fetched site geometry against the Loyola campus map.

Produces `site_verification.png`: a plan view of the building footprints
shaded by height, with the named campus buildings labelled and Sean Earl Field
outlined, so the layout can be compared by eye against
`luc.edu/media/lucedu/pdfs-campusmaps/2025-map-lsc-8.5x11.pdf`.

Also prints the checks that matter for the CFD setup:
  * are the buildings the paper names actually present, and where;
  * the tallest building H, which sets the COST 732 / Tominaga (2008) domain
    guidelines;
  * whether the domain proposed in the manuscript satisfies those guidelines.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

# Buildings the manuscript names, with the position it claims relative to the
# field. Checked against what OSM actually reports.
EXPECTED = {
    "Mertz": "immediately west of the field",
    "Norville": "immediately south of the field",
    "Damen": "~150 m west",
    "Cudahy": "~200 m WNW",
    "Information Commons": "~250 m west",
}


def bearing_name(cx, cy):
    ang = np.degrees(np.arctan2(cy, cx)) % 360.0
    names = ["E", "ENE", "NE", "NNE", "N", "NNW", "NW", "WNW",
             "W", "WSW", "SW", "SSW", "S", "SSE", "SE", "ESE"]
    return names[int((ang + 11.25) // 22.5) % 16]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("geometry_json", nargs="?", default="site_geometry.json")
    ap.add_argument("--out", default="site_verification.png")
    ap.add_argument("--half-size", type=float, default=320.0)
    args = ap.parse_args()

    with open(args.geometry_json) as fh:
        G = json.load(fh)
    buildings = G["buildings"]
    field = np.asarray(G["field"]["ring"], dtype=float)

    print("=" * 74)
    print("SITE VERIFICATION vs the manuscript's description")
    print("=" * 74)
    print(f"origin (Sean Earl Field centroid): "
          f"{G['origin_lat']:.6f}, {G['origin_lon']:.6f}")

    fx, fy = field[:, 0], field[:, 1]
    w, h = fx.max() - fx.min(), fy.max() - fy.min()
    print(f"\nfield polygon: {w:.0f} m (E-W) x {h:.0f} m (N-S)")
    axis = "E-W" if w > h else "N-S"
    print(f"  long axis runs {axis}")
    if axis != "N-S":
        print("  ** MANUSCRIPT SAYS 'the field's long axis runs roughly N-S' "
              "-- OSM disagrees **")

    print("\nnamed buildings the manuscript relies on:")
    for key, claim in EXPECTED.items():
        hits = [b for b in buildings
                if b.get("name") and key.lower() in b["name"].lower()]
        if not hits:
            print(f"  {key:22s} NOT FOUND in OSM   (manuscript: {claim})")
            continue
        b = max(hits, key=lambda b: b["height_m"])
        ring = np.asarray(b["ring"])
        cx, cy = ring[:, 0].mean(), ring[:, 1].mean()
        d = float(np.hypot(cx, cy))
        print(f"  {key:22s} h={b['height_m']:5.1f} m  {d:5.0f} m "
              f"{bearing_name(cx, cy):3s}  ({b['height_source']})")
        print(f"  {'':22s} manuscript: {claim}")

    H = max(b["height_m"] for b in buildings)
    tallest = max(buildings, key=lambda b: b["height_m"])
    print(f"\ntallest building in the box: {tallest.get('name') or '(unnamed)'}"
          f" at H = {H:.1f} m")

    print("\nCOST 732 / Tominaga (2008) domain guidelines for this H:")
    print(f"  5H upstream / lateral / top : {5*H:6.0f} m")
    print(f"  15H downstream              : {15*H:6.0f} m")
    print("  manuscript proposes         : 500 x 500 x 80 m")
    if 5 * H > 80:
        print(f"  ** the 80 m domain top is only {80/H:.1f}H above the "
              f"tallest building; the guideline is 5H = {5*H:.0f} m **")
    # Frontal blockage of the tallest building in a 500 x 80 m cross-section.
    ring = np.asarray(tallest["ring"])
    width = ring[:, 1].max() - ring[:, 1].min()
    blockage = (width * H) / (500.0 * 80.0)
    print(f"  ** blockage of the tallest building alone in a 500x80 m "
          f"section: {100*blockage:.1f}% (guideline < 3%) **")

    # ---- figure ----------------------------------------------------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection

    fig, ax = plt.subplots(figsize=(9.5, 9))
    patches, heights = [], []
    for b in buildings:
        ring = np.asarray(b["ring"])
        if (abs(ring[:, 0].mean()) > args.half_size
                or abs(ring[:, 1].mean()) > args.half_size):
            continue
        patches.append(MplPolygon(ring, closed=True))
        heights.append(b["height_m"])

    pc = PatchCollection(patches, cmap="inferno", edgecolor="k", linewidth=0.3,
                         alpha=0.9)
    pc.set_array(np.array(heights))
    pc.set_clim(0, max(heights))
    ax.add_collection(pc)
    cb = fig.colorbar(pc, ax=ax, shrink=0.75)
    cb.set_label("building height (m)")

    ax.add_patch(MplPolygon(field, closed=True, facecolor="limegreen",
                            edgecolor="darkgreen", linewidth=2.0, alpha=0.65,
                            zorder=5))
    ax.plot(0, 0, "k+", ms=14, mew=2, zorder=6)

    for key in EXPECTED:
        hits = [b for b in buildings
                if b.get("name") and key.lower() in b["name"].lower()]
        if not hits:
            continue
        b = max(hits, key=lambda b: b["height_m"])
        ring = np.asarray(b["ring"])
        cx, cy = ring[:, 0].mean(), ring[:, 1].mean()
        if abs(cx) > args.half_size or abs(cy) > args.half_size:
            continue
        ax.annotate(f"{b['name']}\n{b['height_m']:.0f} m", (cx, cy),
                    fontsize=8, ha="center", va="center", zorder=7,
                    color="white",
                    bbox=dict(boxstyle="round,pad=0.25", fc="black",
                              alpha=0.65, lw=0))

    ax.annotate("Sean Earl Field", (fx.mean(), fy.mean()), fontsize=9,
                ha="center", va="center", zorder=8, weight="bold",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", alpha=0.85,
                          lw=0))
    ax.annotate("Lake Michigan\n(east)", (args.half_size * 0.82, 0),
                fontsize=9, ha="center", color="steelblue", weight="bold")

    ax.set_xlim(-args.half_size, args.half_size)
    ax.set_ylim(-args.half_size, args.half_size)
    ax.set_aspect("equal")
    ax.set_xlabel("x (East), m")
    ax.set_ylabel("y (North), m")
    ax.set_title("Sean Earl Field site model — OSM footprints, height-shaded\n"
                 "origin at field centroid "
                 f"({G['origin_lat']:.5f}, {G['origin_lon']:.5f})")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.out, dpi=160)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
