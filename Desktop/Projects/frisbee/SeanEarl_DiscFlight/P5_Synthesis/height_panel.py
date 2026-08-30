"""The inset panel that makes the vertical structure of a throw legible.

WHY A PANEL AT ALL
------------------
From behind the thrower -- the view this study is rendered in -- height is the
one thing you cannot read. Altitude and downrange distance both project onto
the vertical axis of the frame, so a disc that climbs and a disc that is simply
further away look the same. Everything else in the picture (lateral drift,
landing scatter, which throws are blocked) survives that projection; z does not.

So z gets its own axes, drawn to the side of the 3-D view:

  LEFT   z against downrange, one line per trajectory in the ensemble, growing
         with the animation. This is the flight profile.
  RIGHT  sigma_u(z), the turbulence the site actually puts at each height,
         measured over the playing surface from the CFD. Shares the z axis with
         the left panel, so the height a throw reaches can be read straight
         across into the disturbance it meets there.

That pairing is the point. `FRISBEE_GENERALIZATION.md` sec 3 argues that
throwing low wins twice -- sigma_u grows with height AND the integral scale
grows with height, so a lofted disc gets a stronger disturbance with less
averaging of it. The left panel says how high the throw went; the right says
what that cost.

sigma_u is taken as sqrt(2k/3) from the RANS TKE, which is the isotropic
partition. The gust model the trajectories actually fly through uses the
anisotropic ratios (1.0, 0.75, 0.5) about the local mean-wind direction, so the
streamwise component is larger than this curve and the two lateral components
smaller. The curve is therefore the ISOTROPIC-EQUIVALENT magnitude, labelled as
such, not the streamwise sigma_u the flight sees.

Axes start at zero (both of them). Nothing here is drawn on a suppressed axis.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

BG = "#0b0f14"
FG = "#cfd8e3"
DIM = "#5d6b7d"


def sigma_u_profile(heading, vertical_npz=None, pitch_x=(-31.0, 47.0),
                    pitch_y=(-10.0, 15.0), z_max=60.0):
    """Isotropic-equivalent sigma_u = sqrt(2k/3) vs height over the pitch.

    Averaged over the playing surface only. The full sampling box extends into
    courtyards and wakes off the pitch, where the turbulence is a different
    quantity; the status card's standing warning about the box mean applies
    here too.
    """
    p = Path(vertical_npz) if vertical_npz else ROOT / "P3_CFD" / "vertical_fields.npz"
    F = np.load(p, allow_pickle=False)
    x, y, z = F["x"], F["y"], F["z"]
    key = f"dir_{int(heading):03d}__k"
    if key not in F.files:
        raise KeyError(f"{p.name} has no {key}; has "
                       f"{[f[4:7] for f in F.files if f.endswith('__k')]}")
    k = F[key]
    solid = F[f"dir_{int(heading):03d}__solid"]
    ix = (x >= pitch_x[0]) & (x <= pitch_x[1])
    iy = (y >= pitch_y[0]) & (y <= pitch_y[1])
    sub = k[np.ix_(ix, iy)]
    msk = ~solid[np.ix_(ix, iy)]
    prof = np.array([np.nanmean(sub[:, :, kz][msk]) for kz in range(len(z))])
    keep = z <= z_max
    return z[keep], np.sqrt(np.maximum(2.0 / 3.0 * prof[keep], 0.0))


def render_panel(tracks, frac, launch, axis, heading, colour, px=(660, 380),
                 dpi=110, z_lim=None, x_lim=None, sigma_cache={}):
    """Return an (H, W, 4) uint8 RGBA image of the two-panel inset.

    `frac` is the same flight fraction the 3-D view is drawn at, so the profile
    grows with the trajectories rather than being revealed all at once.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if heading not in sigma_cache:
        sigma_cache[heading] = sigma_u_profile(heading)
    zp, sig = sigma_cache[heading]

    w, h = px
    fig = plt.figure(figsize=(w / dpi, h / dpi), dpi=dpi, facecolor=BG)
    gs = fig.add_gridspec(1, 2, width_ratios=[3.0, 1.0], wspace=0.06,
                          left=0.11, right=0.97, top=0.92, bottom=0.17)
    axL = fig.add_subplot(gs[0, 0], facecolor=BG)
    axR = fig.add_subplot(gs[0, 1], facecolor=BG, sharey=axL)

    launch = np.asarray(launch, dtype=float)
    zmax_seen = 1.0
    xmax_seen = 1.0
    for r in tracks:
        n = max(2, int(len(r) * frac))
        seg = np.asarray(r[:n], dtype=float)
        down = (seg[:, 0] - launch[0]) * axis[0] + (seg[:, 1] - launch[1]) * axis[1]
        axL.plot(down, seg[:, 2], color=colour, lw=0.9, alpha=0.45,
                 solid_capstyle="round")
        zmax_seen = max(zmax_seen, float(np.nanmax(seg[:, 2])))
        xmax_seen = max(xmax_seen, float(np.nanmax(down)))

    z_hi = z_lim if z_lim else max(6.0, 1.25 * zmax_seen)
    x_hi = x_lim if x_lim else max(20.0, 1.10 * xmax_seen)
    axL.set_xlim(0, x_hi)
    axL.set_ylim(0, z_hi)                      # zero-based, never data-ranged
    axL.set_xlabel("downrange  (m)", color=FG, fontsize=9)
    axL.set_ylabel("height  z  (m)", color=FG, fontsize=9)

    axR.plot(sig, zp, color="#e8b03a", lw=1.8)
    axR.set_xlim(0, max(0.2, 1.15 * float(np.nanmax(sig[zp <= z_hi]))))
    axR.set_xlabel(r"$\sigma_u$  (m/s)", color=FG, fontsize=9)
    axR.tick_params(labelleft=False)

    # The band the ensemble actually occupied, read across both panels.
    for ax in (axL, axR):
        ax.axhspan(0, zmax_seen, color=colour, alpha=0.07, lw=0)
        ax.axhline(zmax_seen, color=colour, lw=0.8, ls=(0, (4, 3)), alpha=0.8)
        for s in ax.spines.values():
            s.set_color(DIM)
        ax.tick_params(colors=FG, labelsize=8, length=3)
        ax.grid(True, color=DIM, alpha=0.22, lw=0.5)

    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)
    return buf


def composite(frame_png, panel_rgba, margin=26, corner="lower_right"):
    """Alpha-composite the panel onto a rendered frame, in place on disk."""
    from PIL import Image
    base = Image.open(frame_png).convert("RGBA")
    pan = Image.fromarray(panel_rgba, mode="RGBA")
    W, H = base.size
    w, h = pan.size
    if corner == "lower_right":
        box = (W - w - margin, H - h - margin)
    elif corner == "upper_right":
        box = (W - w - margin, margin)
    else:
        box = (margin, H - h - margin)
    base.alpha_composite(pan, dest=box)
    base.convert("RGB").save(frame_png)
