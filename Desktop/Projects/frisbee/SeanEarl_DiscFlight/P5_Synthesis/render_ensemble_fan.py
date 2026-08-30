"""Render Monte Carlo ensemble fans over the real site geometry, one per wind condition.

The camera is FIXED behind the thrower, on the mean flight axis measured from the
ensemble actually flown, so the whole flight stays in frame from release to landing and
nothing that moves in the picture is the camera. Rendering one heading per clip lets the
three computed wind conditions be compared like for like; passing several headings
overlays them in one clip instead.

NOTE ON THE DIRECTIONAL CLAIM. These clips are rendered on the production field, which
is what the manuscript reports. That field resolves 270 deg as the worst of the three
headings (p = 0.023 and 0.0052), but that verdict does NOT survive a mesh refinement --
see GRID_INDEPENDENCE_2026-08-30.md. Read the clips as "what the reported field looks
like", not as evidence that the headings differ.

WHAT IS AND IS NOT DATA HERE
----------------------------
Everything with a physical position is computed, not drawn:

  * the buildings are `P2_Geometry/seanEarl.stl` -- 116 OSM footprints, 6,324
    triangles, extruded to tagged or type-default heights;
  * every trajectory is a fresh integration of the same 6-DOF model the paper uses
    (`P1_AnalyticalModel/disc_6dof.py`, `potts_tabulated` coefficients), flown through
    the same time-averaged RANS field (`P3_CFD/flow_fields.npz`) with the same
    `sweeping` gust model and the same release-scatter distribution
    (`monte_carlo.sample_release`) as the 97,200-throw sweep;
  * the landing scatter is those trajectories' own endpoints, and the 95% ellipse is
    fitted to them;
  * the ground plane is coloured by the disc-altitude horizontal wind speed at z = 2 m
    from the same field.

The only choices made for looks are camera, colour and opacity. No trajectory is
smoothed, resampled for shape, or hand-placed. Throws that land inside a building are
drawn in a separate colour rather than deleted, because "blocked" is a real outcome of
the model (10.2% of the sweep) and hiding them would flatter the picture.

Seeds are fixed, so the same fan renders every time.

Usage
-----
    python render_ensemble_fan.py --n-throws 60 --fps 30 --seconds 14
    python render_ensemble_fan.py --preview          # one still, no video

Output: ensemble_fan.mp4 (+ frames/ and a still) under --outdir.
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "P1_AnalyticalModel"))
sys.path.insert(0, str(ROOT / "P4_Coupling"))

# Default clip: the two headings the paper compares against each other. Override with
# --headings to render a single wind condition, or all three.
HEADINGS = (210, 270)
HEADING_LABEL = {90: "90 deg  easterly, over-lake",
                 210: "210 deg  SSW, prevailing",
                 270: "270 deg  westerly, urban"}

# A launch point near the middle of the pitch, on the sweep's own grid.
LAUNCH = np.array([-12.0, -5.0, 1.5])
BEARING_DEG = 60.0
THROW = "right_backhand"
SEED = 20260829
TRACK_SAMPLES = 72      # points kept per trajectory for rendering, see fly_ensemble

# Colour per heading, kept fixed across every clip so the same wind condition is the
# same colour whether it is rendered alone or overlaid on another.
DEFAULT_COLOUR = {90: (0.55, 0.82, 0.60),     # easterly, over-lake
                  210: (0.30, 0.69, 0.94),    # SSW, prevailing
                  270: (0.97, 0.45, 0.33)}    # westerly, urban
COLOUR = dict(DEFAULT_COLOUR)
BLOCKED_COLOUR = (0.55, 0.16, 0.20)

PITCH_X = (-31.0, 47.0)
PITCH_Y = (-10.0, 15.0)


def fly_ensemble(heading, n_throws, field_path, seed=SEED):
    """Fly `n_throws` from the same launch point, returning full trajectories.

    Uses monte_carlo's own release sampler, gust model and integrator settings, so
    this ensemble is drawn from the same distribution as the production sweep rather
    than from a lookalike.
    """
    import disc_6dof as d6
    import monte_carlo as mc
    import synthetic_turbulence as st

    field = mc._get_field(str(field_path), heading)
    rng = np.random.default_rng(seed)

    tracks, landings, blocked = [], [], []
    for _ in range(n_throws):
        gusts = st.CorrelatedGusts(field, rng, t_max=8.0, ref_pos=LAUNCH,
                                   model="sweeping",
                                   sigma_ratios=st.SIGMA_RATIOS)
        wind = mc.wind_with_gusts(field, gusts)
        r0, v0, q0, w0 = mc.sample_release(
            throw_type=THROW, heading_deg=BEARING_DEG,
            launch_pos=LAUNCH, rng=rng)
        res = d6.simulate(r0, v0, q0, w0,
                          coeffs=mc.DEFAULT_COEFFS, geom=mc.ULTRASTAR,
                          wind=wind, t_max=8.0,
                          rtol=mc.MC_RTOL, atol=mc.MC_ATOL,
                          compute_diag=False)
        r = np.asarray(res.r)
        # A tube through ~270 integrator samples carries several times the geometry a
        # smooth curve needs, and rebuilding it dominates the per-frame render cost.
        # Subsample to a fixed count, always keeping the first and last point so the
        # launch and the landing stay exact.
        if len(r) > TRACK_SAMPLES:
            idx = np.unique(np.linspace(0, len(r) - 1, TRACK_SAMPLES).astype(int))
            r = r[idx]
        tracks.append(r)
        land = r[-1]
        landings.append(land)
        is_blocked = bool(field.is_solid(land[0], land[1])) \
            if hasattr(field, "is_solid") else False
        blocked.append(is_blocked)

    return tracks, np.asarray(landings), np.asarray(blocked)


def confidence_ellipse(points, n_sigma=2.4477):
    """95% ellipse of a 2-D point cloud (2.4477 sigma for 2 dof)."""
    pts = np.asarray(points)[:, :2]
    pts = pts[np.isfinite(pts).all(axis=1)]
    if len(pts) < 3:
        return None
    mu = pts.mean(axis=0)
    cov = np.cov(pts.T)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    t = np.linspace(0, 2 * np.pi, 200)
    circle = np.column_stack([np.cos(t), np.sin(t)])
    ell = circle @ np.diag(n_sigma * np.sqrt(np.maximum(vals, 0))) @ vecs.T + mu
    return ell


def ground_speed_texture(field_path, heading):
    """Horizontal wind speed at z = 2 m, as a structured grid for the ground plane."""
    import pyvista as pv
    sys.path.insert(0, str(HERE))
    from conventions import disc_altitude_index, horizontal_speed

    F = np.load(field_path, allow_pickle=False)
    x, y = F["x"], F["y"]
    zk, _, _ = disc_altitude_index(F["z"], 2.0)
    spd = horizontal_speed(F[f"dir_{int(heading):03d}__U"][:, :, zk, :])
    solid = F[f"dir_{int(heading):03d}__solid"]

    grid = pv.StructuredGrid()
    XX, YY = np.meshgrid(x, y, indexing="ij")
    ZZ = np.full_like(XX, 0.02)          # just above the ground plane
    grid.points = np.c_[XX.ravel(order="F"), YY.ravel(order="F"), ZZ.ravel(order="F")]
    grid.dimensions = (len(x), len(y), 1)
    s = spd.copy()
    s[solid.astype(bool)] = np.nan       # do not paint speed inside buildings
    grid["speed"] = s.ravel(order="F")
    return grid


def build_scene(args):
    import pyvista as pv

    pv.global_theme.allow_empty_mesh = True
    off = not args.interactive
    pl = pv.Plotter(off_screen=off, window_size=[args.width, args.height],
                    lighting="three lights")
    pl.set_background("#0b0f14", top="#1c2732")

    # --- site geometry -------------------------------------------------------
    stl = ROOT / "P2_Geometry" / "seanEarl.stl"
    if stl.exists():
        buildings = pv.read(str(stl))
        pl.add_mesh(buildings, color="#4a5361", specular=0.25, specular_power=12,
                    smooth_shading=False, ambient=0.28, diffuse=0.75,
                    name="buildings")
    else:
        print(f"  !! {stl} missing -- rendering without buildings")

    # --- ground, coloured by disc-altitude wind speed ------------------------
    field_path = Path(args.field) if args.field else ROOT / "P3_CFD" / "flow_fields.npz"
    ground = ground_speed_texture(field_path, args.ground_heading)
    pl.add_mesh(ground, scalars="speed", cmap="cividis", opacity=0.92,
                nan_opacity=0.0, show_scalar_bar=True,
                scalar_bar_args=dict(title="2 m wind speed  (m/s)",
                                     color="#cfd8e3", n_labels=4,
                                     position_x=0.66, position_y=0.035,
                                     width=0.27, height=0.045,
                                     title_font_size=13, label_font_size=11))

    # --- the pitch outline ---------------------------------------------------
    px, py = PITCH_X, PITCH_Y
    outline = pv.MultipleLines(points=np.array([
        [px[0], py[0], 0.05], [px[1], py[0], 0.05],
        [px[1], py[1], 0.05], [px[0], py[1], 0.05], [px[0], py[0], 0.05]]))
    pl.add_mesh(outline, color="#e8eef6", line_width=2.5, name="pitch")

    # --- the ensembles -------------------------------------------------------
    ensembles = {}
    for h in args.headings:
        print(f"  flying {args.n_throws} throws at heading {h} deg ...")
        tracks, landings, blocked = fly_ensemble(h, args.n_throws, field_path)
        ensembles[h] = (tracks, landings, blocked)
        print(f"    {int((~blocked).sum())} delivered, {int(blocked.sum())} blocked")

    return pl, ensembles


_TUBE_CACHE = {}


def add_tracks(pl, ensembles, frac, args):
    """Draw each ensemble's trajectories up to fraction `frac` of their flight.

    Once the flight is fully drawn the geometry stops changing, so the finished tubes
    are built once and reused for every remaining frame. Rebuilding them each frame
    dominated the render cost.
    """
    import pyvista as pv

    full = frac >= 1.0
    for h, (tracks, landings, blocked) in ensembles.items():
        col = COLOUR[h]
        segs = []
        for r, blk in zip(tracks, blocked):
            n = max(2, int(len(r) * frac))
            seg = r[:n]
            if len(seg) < 2:
                continue
            segs.append((seg, blk))
        if not segs:
            continue

        def _tube(which, want_blocked):
            key = (which, h)
            if full and key in _TUBE_CACHE:
                return _TUBE_CACHE[key]
            parts = [pv.MultipleLines(points=s) for s, b in segs
                     if bool(b) == want_blocked]
            if not parts:
                return None
            merged = parts[0].merge(parts[1:]) if len(parts) > 1 else parts[0]
            tube = merged.tube(radius=args.tube_radius, n_sides=6)
            if full:
                _TUBE_CACHE[key] = tube
            return tube

        tube = _tube("open", False)
        if tube is not None:
            pl.add_mesh(tube, color=col, opacity=0.55, specular=0.5,
                        ambient=0.35, diffuse=0.8,
                        name=f"tracks_{h}")
        tube = _tube("blk", True)
        if tube is not None:
            pl.add_mesh(tube, color=BLOCKED_COLOUR, opacity=0.5,
                        name=f"tracks_blocked_{h}")

        # Leading discs: where each throw is right now.
        heads = np.array([r[max(1, int(len(r) * frac)) - 1] for r in tracks])
        pl.add_mesh(pv.PolyData(heads), color=col, point_size=7.0,
                    render_points_as_spheres=True, name=f"heads_{h}")


def add_landings(pl, ensembles, alpha, args):
    import pyvista as pv

    for h, (tracks, landings, blocked) in ensembles.items():
        good = landings[~blocked]
        if len(good) == 0:
            continue
        pts = good.copy()
        pts[:, 2] = 0.08
        pl.add_mesh(pv.PolyData(pts), color=COLOUR[h], point_size=13.0 * alpha + 1,
                    render_points_as_spheres=True, opacity=min(1.0, alpha),
                    name=f"land_{h}")
        if len(blocked) and blocked.any():
            bp = landings[blocked].copy()
            bp[:, 2] = 0.08
            pl.add_mesh(pv.PolyData(bp), color=BLOCKED_COLOUR,
                        point_size=11.0 * alpha + 1, opacity=min(1.0, alpha),
                        render_points_as_spheres=True, name=f"landblk_{h}")
        ell = confidence_ellipse(good)
        if ell is not None and alpha > 0.55:
            e3 = np.column_stack([ell, np.full(len(ell), 0.12)])
            pl.add_mesh(pv.MultipleLines(points=e3), color=COLOUR[h],
                        line_width=3.5, opacity=min(1.0, (alpha - 0.55) / 0.45),
                        name=f"ell_{h}")


def caption(pl, ensembles, args):
    """Caption built from the ensemble actually rendered, not typed in."""
    lines = []
    for h in ensembles:
        tracks, landings, blocked = ensembles[h]
        good = landings[~blocked]
        if len(good) >= 3:
            cov = np.cov(good[:, :2].T)
            sig = float(np.sqrt(0.5 * np.trace(cov)))
        else:
            sig = float("nan")
        lines.append(f"{HEADING_LABEL[h]}   sigma = {sig:.2f} m   "
                     f"blocked {100 * blocked.mean():.0f}%")
    pl.add_text("\n".join(lines), position="upper_left", font_size=11,
                color="#dfe7f0", name="caption")
    src = Path(args.field).name if args.field else "flow_fields.npz"
    pl.add_text(f"{args.n_throws} throws per condition  |  6-DOF x RANS x sweeping "
                f"gusts  |  {src}",
                position="lower_left", font_size=9, color="#8fa0b3", name="sub")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-throws", type=int, default=60)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--seconds", type=float, default=12.0)
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--height", type=int, default=1080)
    ap.add_argument("--tube-radius", type=float, default=0.32)
    ap.add_argument("--ground-heading", type=int, default=None,
                    help="which heading's wind field paints the ground "
                         "(default: the first heading rendered)")
    ap.add_argument("--headings", nargs="+", type=int, default=list(HEADINGS),
                    help="wind conditions to render in this clip. One heading "
                         "gives a single-condition video; two or more overlay them "
                         "for comparison.")
    ap.add_argument("--field", default=None,
                    help="flow field npz (default: the production field the "
                         "manuscript reports)")
    ap.add_argument("--name", default="ensemble_fan",
                    help="output basename, so several wind conditions can be "
                         "rendered side by side without overwriting each other")
    # Fixed behind-the-thrower camera.
    ap.add_argument("--cam-back", type=float, default=38.0,
                    help="metres behind the launch point, along the flight axis")
    ap.add_argument("--cam-height", type=float, default=25.0,
                    help="camera height in metres")
    ap.add_argument("--focus-frac", type=float, default=0.55,
                    help="focus this fraction of the way down the flight")
    ap.add_argument("--focus-z", type=float, default=1.5)
    ap.add_argument("--outdir", default="video")
    ap.add_argument("--preview", action="store_true", help="one still, no video")
    ap.add_argument("--interactive", action="store_true")
    args = ap.parse_args()

    for h in args.headings:
        COLOUR.setdefault(h, (0.65, 0.78, 0.55))
    if args.ground_heading is None:
        args.ground_heading = args.headings[0]

    outdir = HERE / args.outdir
    frames = outdir / "frames"
    frames.mkdir(parents=True, exist_ok=True)

    pl, ensembles = build_scene(args)

    def thrower_camera():
        """Fixed camera behind the thrower, looking down the mean flight axis.

        The axis is measured from the ensemble that was actually flown -- launch
        point to mean landing -- not assumed from the nominal bearing, because the
        wind deflects the flight by several degrees and the amount differs per
        heading. So the camera sits behind the thrower for the throw as the wind
        actually shapes it, which is the whole point of the view.

        Fixed for the entire clip: no orbit. The full flight stays in frame from
        release to landing, and nothing that moves in the picture is the camera.
        """
        launch = LAUNCH.astype(float)
        # Pool the landings of every heading being rendered, so a multi-heading
        # comparison gets one shared camera rather than one that favours a heading.
        lands = np.vstack([np.asarray(ens[1])[~np.asarray(ens[2])][:, :2]
                           for ens in ensembles.values()
                           if (~np.asarray(ens[2])).any()])
        mean_land = lands.mean(axis=0)

        axis = mean_land - launch[:2]
        span = float(np.linalg.norm(axis))
        axis = axis / max(span, 1e-9)

        # Focus part-way down the flight so both ends sit inside the frame.
        focus = np.array([launch[0] + axis[0] * span * args.focus_frac,
                          launch[1] + axis[1] * span * args.focus_frac,
                          args.focus_z])
        # Camera behind the launch point, along the reversed flight axis.
        pos = np.array([launch[0] - axis[0] * args.cam_back,
                        launch[1] - axis[1] * args.cam_back,
                        args.cam_height])
        pl.camera_position = [tuple(pos), tuple(focus), (0, 0, 1)]
        return launch, mean_land, span

    launch, mean_land, span = thrower_camera()
    print(f"  camera: fixed behind the thrower, {args.cam_back:.0f} m back and "
          f"{args.cam_height:.0f} m up,")
    print(f"          looking down a {span:.1f} m mean flight axis toward "
          f"({mean_land[0]:+.1f}, {mean_land[1]:+.1f})")

    if args.preview:
        add_tracks(pl, ensembles, 1.0, args)
        add_landings(pl, ensembles, 1.0, args)
        caption(pl, ensembles, args)
        still = outdir / f"{args.name}_still.png"
        pl.screenshot(str(still))
        print(f"wrote {still}")
        return

    n_frames = int(args.fps * args.seconds)
    print(f"  rendering {n_frames} frames at {args.width}x{args.height} ...")
    for i in range(n_frames):
        u = i / max(1, n_frames - 1)
        # Flight grows over the first 55%, landings settle over the next 25%.
        # The camera does not move.
        frac = min(1.0, u / 0.55)
        alpha = 0.0 if u < 0.55 else min(1.0, (u - 0.55) / 0.25)
        add_tracks(pl, ensembles, frac, args)
        if alpha > 0:
            add_landings(pl, ensembles, alpha, args)
        caption(pl, ensembles, args)
        pl.screenshot(str(frames / f"f{i:04d}.png"))
        if i % 25 == 0:
            print(f"    frame {i}/{n_frames}")

    mp4 = outdir / f"{args.name}.mp4"
    cmd = ["ffmpeg", "-y", "-framerate", str(args.fps),
           "-i", str(frames / "f%04d.png"),
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           "-movflags", "+faststart", str(mp4)]
    print("  " + " ".join(cmd))
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stderr[-3000:])
        raise SystemExit("ffmpeg failed")
    print(f"wrote {mp4}")


if __name__ == "__main__":
    main()
