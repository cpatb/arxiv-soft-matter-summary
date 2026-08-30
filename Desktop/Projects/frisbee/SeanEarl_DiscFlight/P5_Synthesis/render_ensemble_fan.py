"""Render Monte Carlo ensemble fans over the real site geometry, one per wind
condition, from behind the thrower at the west end of the field looking east.

THE VIEW
--------
The thrower stands on the pitch centreline 3 m inside the WESTERN goal line and
throws EAST, the length of the field, toward Campus Towers and the lake beyond
it. The camera sits behind them on the flight axis and does not move. That is
the geometry the site's headline result is about: a 45 m block nine metres off
the eastern edge, not the lakefront boundary layer, is what the disc flies into.
The frame, the bearing and the site frame are all fixed in `site_view.py`, with
the sign convention checked against the CFD rather than assumed.

Throw bearing 0 deg (= due east) is also the bearing the production Monte Carlo
sweep threw at, so these clips show the throw the paper's numbers describe.

HEIGHT
------
From behind the thrower, altitude and downrange distance project onto the same
axis of the frame, so z is the one thing the main view cannot show. It gets its
own inset (`height_panel.py`): the flight profile z(downrange) beside the
site's own sigma_u(z), on a shared height axis, so how high a throw went can be
read straight across into the disturbance it met up there. `--camera sideline`
renders the same ensemble from the south touchline, where height is legible
directly.

WHAT IS AND IS NOT DATA HERE
----------------------------
Everything with a physical position is computed, not drawn:

  * the buildings are `P2_Geometry/seanEarl.stl` -- 116 OSM footprints, 6,324
    triangles, extruded to tagged or type-default heights;
  * every trajectory is a fresh integration of the same 6-DOF model the paper
    uses (`P1_AnalyticalModel/disc_6dof.py`, `potts_tabulated` coefficients),
    flown through the same time-averaged RANS field with the same `sweeping`
    gust model and the same release-scatter distribution as the 97,200-throw
    sweep;
  * the landing scatter is those trajectories' own endpoints, and the 95%
    ellipse is fitted to them;
  * the ground plane is coloured by the disc-altitude horizontal wind speed at
    z = 2 m from the same field, and the sigma_u(z) curve in the inset is
    sqrt(2k/3) from the same solution.

The only choices made for looks are camera, colour and opacity. No trajectory
is smoothed, resampled for shape, or hand-placed. Throws that land inside a
building are drawn in a separate colour rather than deleted, because "blocked"
is a real outcome of the model and hiding it would flatter the picture.

THROW LENGTH. `--throw` selects from `P1_AnalyticalModel/throw_catalogue.py`.
The default `right_backhand` is the published distribution the sweep used and
carries 27 m. The longer entries are calibrated to a target still-air RANGE,
with the release speed solved for -- they are "a throw that carries 50 m in
this model", not a measured release. Read `throw_catalogue.py`'s docstring
before quoting any of their release parameters.

FIELD RESOLUTION IN Z. The production `flow_fields.npz` is sampled at
z = 1..6 m only, so a disc above 6 m flies through wind frozen at its 6 m
value. That is fine for the 27 m backhand (apex 3.8 m) and NOT fine for a
lofted throw. Pass `--field ../P3_CFD/vertical_fields.npz` for those: same CFD
solution, sampled to z = 100 m at 2 m horizontal spacing instead of 1 m. The
renderer refuses to fly a throw whose ensemble apex leaves the field's z range
unless `--allow-clamped` is given.

NOTE ON THE DIRECTIONAL CLAIM. These clips are rendered on the production
field, which is what the manuscript reports. That field resolves 270 deg as the
worst of the three headings (p = 0.023 and 0.0052), but that verdict does NOT
survive a mesh refinement -- see GRID_INDEPENDENCE_2026-08-30.md. Read the
clips as "what the reported field looks like", not as evidence that the
headings differ.

Seeds are fixed, so the same fan renders every time.

Usage
-----
    python render_ensemble_fan.py --headings 90 --name wind_090_easterly
    python render_ensemble_fan.py --camera sideline --throw huck_50m \
        --field ../P3_CFD/vertical_fields.npz
    python render_ensemble_fan.py --preview          # one still, no video
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "P1_AnalyticalModel"))
sys.path.insert(0, str(ROOT / "P4_Coupling"))

import site_view as sv                                    # noqa: E402
import height_panel as hp                                 # noqa: E402

HEADINGS = (210, 270)
HEADING_LABEL = {90: "90 deg  easterly, off the lake",
                 210: "210 deg  SSW, prevailing",
                 270: "270 deg  westerly, over the campus"}

LAUNCH = sv.LAUNCH
BEARING_DEG = sv.EAST
SEED = 20260829
TRACK_SAMPLES = 72      # points kept per trajectory for rendering

DEFAULT_COLOUR = {90: (0.55, 0.82, 0.60),
                  210: (0.30, 0.69, 0.94),
                  270: (0.97, 0.45, 0.33)}
COLOUR = dict(DEFAULT_COLOUR)
BLOCKED_COLOUR = (0.55, 0.16, 0.20)

PITCH_X, PITCH_Y = sv.PITCH_X, sv.PITCH_Y


def _draw(rng, mean, std, rng_range):
    """Draw one release angle: uniform over `rng_range` if given, else normal.

    A Gaussian release scatter is the right model for ONE thrower repeating ONE
    throw, which is what the production sweep needed. It is the wrong model for
    showing what the site does to a RANGE of throws: 95% of its draws sit
    within +/-2 sigma of a single nominal angle, so a fan built from it shows
    one throw blurred rather than the envelope of throws a player actually has.
    A uniform sweep over a stated range is an envelope, and it is labelled as
    one in the caption so nobody reads it as a probability.
    """
    if rng_range is None:
        return rng.normal(mean, std)
    lo, hi = rng_range
    return rng.uniform(lo, hi)


def fly_ensemble(heading, n_throws, field_path, throw, launch, bearing,
                 seed=SEED, u_ref=None, elev_range=None, bearing_range=None,
                 roll_range=None):
    """Fly `n_throws` from the same launch point, returning full trajectories.

    Uses monte_carlo's own gust model and integrator settings and the release
    sampler's own distributions, so this ensemble is drawn from the same
    machinery as the production sweep rather than from a lookalike.
    """
    import disc_6dof as d6
    import monte_carlo as mc
    import synthetic_turbulence as st
    from throw_catalogue import build_catalogue

    cat = build_catalogue()
    if throw not in cat:
        raise SystemExit(f"unknown throw {throw!r}; have {sorted(cat)}")
    dist = cat[throw]

    field = mc._get_field(str(field_path), heading)
    if u_ref is not None:
        # Re-express the whole solution at a different reference wind speed.
        # Legal because P6_Theory/wind_scaling.py measured it: the identical
        # mesh re-solved at 2.5/5/10 m/s reproduces itself to a median 0.06-0.14%
        # of RMS over the pitch once the expected powers are divided out.
        sys.path.insert(0, str(ROOT / "P6_Theory"))
        from wind_scaling import ScaledField
        field = ScaledField(field, float(u_ref))
    rng = np.random.default_rng(seed)

    tracks, landings, blocked = [], [], []
    for _ in range(n_throws):
        gusts = st.CorrelatedGusts(field, rng, t_max=12.0, ref_pos=launch,
                                   model="sweeping",
                                   sigma_ratios=st.SIGMA_RATIOS)
        wind = mc.wind_with_gusts(field, gusts)
        r0, v0, q0, w0 = d6.build_initial_conditions(
            launch_pos=launch,
            speed=rng.normal(dist["speed_mean"], dist["speed_std"]),
            heading_deg=bearing + (0.0 if bearing_range is None
                                   else rng.uniform(*bearing_range)),
            elevation_deg=_draw(rng, dist["elevation_mean"],
                                dist["elevation_std"], elev_range),
            release_alpha_deg=rng.normal(dist["release_alpha_mean"],
                                         dist["release_alpha_std"]),
            roll_deg=_draw(rng, dist["roll_mean"], dist["roll_std"],
                           roll_range),
            spin_rad_s=max(20.0, rng.normal(dist["spin_mean"],
                                            dist["spin_std"])),
            handedness=dist["handedness"])
        res = d6.simulate(r0, v0, q0, w0,
                          coeffs=mc.DEFAULT_COEFFS, geom=mc.ULTRASTAR,
                          wind=wind, t_max=12.0,
                          rtol=mc.MC_RTOL, atol=mc.MC_ATOL,
                          compute_diag=False)
        r = np.asarray(res.r)
        # A tube through ~270 integrator samples carries several times the
        # geometry a smooth curve needs, and rebuilding it dominates the
        # per-frame render cost. Subsample to a fixed count, always keeping the
        # first and last point so the launch and the landing stay exact.
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


def check_z_coverage(ensembles, field_path, allow_clamped):
    """Refuse to render a throw that spends its flight above the sampled box.

    The field is clamped, not extrapolated, outside its z range, so a lofted
    disc would fly through wind frozen at the top sampled level. That is a
    silent corruption of exactly the quantity these clips are meant to show, so
    it is an error by default rather than a warning nobody reads.
    """
    z = np.load(field_path, allow_pickle=False)["z"]
    z_top = float(z.max())
    worst = 0.0
    n_pts = 0
    n_above = 0
    for tracks, _, _ in ensembles.values():
        for r in tracks:
            worst = max(worst, float(np.nanmax(r[:, 2])))
            n_pts += len(r)
            n_above += int((r[:, 2] > z_top).sum())
    frac_above = n_above / max(n_pts, 1)
    msg = (f"ensemble apex {worst:.1f} m vs field top z = {z_top:.1f} m; "
           f"{100 * frac_above:.1f}% of trajectory samples are above the "
           f"sampled box and fly on CLAMPED wind")
    if worst > z_top and not allow_clamped:
        raise SystemExit(
            "REFUSING TO RENDER: " + msg + ".\n"
            "  Use --field ../P3_CFD/vertical_fields.npz (same CFD solution, "
            "sampled to z = 100 m), or pass --allow-clamped if you accept it.")
    if worst > z_top:
        print(f"  !! CLAMPED (allowed): {msg}")
    else:
        print(f"  z coverage OK: apex {worst:.1f} m, field top {z_top:.1f} m")
    return worst, z_top, frac_above


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
    return circle @ np.diag(n_sigma * np.sqrt(np.maximum(vals, 0))) @ vecs.T + mu


def ground_speed_texture(field_path, heading, u_ref=None):
    """Horizontal wind speed at z = 2 m, as a structured grid for the ground.

    `u_ref` rescales it, for the same reason the sigma_u panel is rescaled: the
    ground colour and the trajectories have to describe the same wind.
    """
    import pyvista as pv
    from conventions import disc_altitude_index, horizontal_speed

    F = np.load(field_path, allow_pickle=False)
    x, y = F["x"], F["y"]
    zk, _, _ = disc_altitude_index(F["z"], 2.0)
    spd = horizontal_speed(F[f"dir_{int(heading):03d}__U"][:, :, zk, :])
    if u_ref is not None:
        spd = spd * (float(u_ref) / 5.0)      # U ~ Uref, exactly
    solid = F[f"dir_{int(heading):03d}__solid"]

    grid = pv.StructuredGrid()
    XX, YY = np.meshgrid(x, y, indexing="ij")
    ZZ = np.full_like(XX, 0.02)
    grid.points = np.c_[XX.ravel(order="F"), YY.ravel(order="F"), ZZ.ravel(order="F")]
    grid.dimensions = (len(x), len(y), 1)
    s = spd.copy()
    s[solid.astype(bool)] = np.nan
    grid["speed"] = s.ravel(order="F")
    return grid


def build_scene(args):
    import pyvista as pv

    pv.global_theme.allow_empty_mesh = True
    off = not args.interactive
    pl = pv.Plotter(off_screen=off, window_size=[args.width, args.height],
                    lighting="three lights")
    pl.set_background("#0b0f14", top="#1c2732")

    stl = ROOT / "P2_Geometry" / "seanEarl.stl"
    if stl.exists():
        pl.add_mesh(pv.read(str(stl)), color="#4a5361", specular=0.25,
                    specular_power=12, smooth_shading=False, ambient=0.28,
                    diffuse=0.75, name="buildings")
    else:
        print(f"  !! {stl} missing -- rendering without buildings")

    field_path = Path(args.field) if args.field else ROOT / "P3_CFD" / "flow_fields.npz"
    sv.assert_launch_clear(field_path, args.launch, args.headings)

    ground = ground_speed_texture(field_path, args.ground_heading, args.u_ref)
    # Clip the colour range at `--speed-max`. The full box reaches ~4.5 m/s in
    # accelerated corners well off the pitch; letting those set the scale
    # flattens the 0-2 m/s range the disc actually flies in, which is the range
    # the picture is about. The clip is stated on the bar rather than hidden:
    # everything above it is drawn in one flat colour.
    pl.add_mesh(ground, scalars="speed", cmap="cividis", opacity=0.92,
                clim=(0.0, args.speed_max), above_color="#fff3b0",
                nan_opacity=0.0, show_scalar_bar=True,
                scalar_bar_args=dict(title="2 m wind speed  (m/s)",
                                     color="#cfd8e3", n_labels=3,
                                     position_x=0.020, position_y=0.078,
                                     width=0.20, height=0.028,
                                     title_font_size=11, label_font_size=10,
                                     vertical=False))

    px, py = PITCH_X, PITCH_Y
    pl.add_mesh(pv.MultipleLines(points=np.array([
        [px[0], py[0], 0.05], [px[1], py[0], 0.05],
        [px[1], py[1], 0.05], [px[0], py[1], 0.05], [px[0], py[0], 0.05]])),
        color="#e8eef6", line_width=2.5, name="pitch")

    # A marker where the thrower stands, so the viewpoint is unambiguous.
    pl.add_mesh(pv.Cylinder(center=(args.launch[0], args.launch[1], 0.85),
                            direction=(0, 0, 1), radius=0.16, height=1.7),
                color="#e8b03a", ambient=0.45, name="thrower")

    ensembles = {}
    for h in args.headings:
        print(f"  flying {args.n_throws} {args.throw} at heading {h} deg ...")
        tracks, landings, blocked = fly_ensemble(
            h, args.n_throws, field_path, args.throw, args.launch,
            args.bearing, u_ref=args.u_ref, elev_range=args.elev_range,
            bearing_range=args.bearing_range, roll_range=args.roll_range)
        ensembles[h] = (tracks, landings, blocked)
        print(f"    {int((~blocked).sum())} delivered, {int(blocked.sum())} blocked")

    check_z_coverage(ensembles, field_path, args.allow_clamped)
    return pl, ensembles


_TUBE_CACHE = {}


def add_tracks(pl, ensembles, frac, args):
    """Draw each ensemble's trajectories up to fraction `frac` of the flight.

    Once the flight is fully drawn the geometry stops changing, so the finished
    tubes are built once and reused. Rebuilding them each frame dominated the
    render cost.
    """
    import pyvista as pv

    full = frac >= 1.0
    for h, (tracks, landings, blocked) in ensembles.items():
        col = COLOUR[h]
        segs = []
        for r, blk in zip(tracks, blocked):
            n = max(2, int(len(r) * frac))
            seg = r[:n]
            if len(seg) >= 2:
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
                        ambient=0.35, diffuse=0.8, name=f"tracks_{h}")
        tube = _tube("blk", True)
        if tube is not None:
            pl.add_mesh(tube, color=BLOCKED_COLOUR, opacity=0.5,
                        name=f"tracks_blocked_{h}")

        # Ground shadow: the same trajectories projected onto z = 0. Reading
        # the gap between a track and its shadow is the only cue to height the
        # behind-the-thrower view has, so it is drawn rather than left out.
        if args.shadows:
            shad = []
            for s, b in segs:
                q = np.array(s, dtype=float, copy=True)
                q[:, 2] = 0.06
                shad.append(pv.MultipleLines(points=q))
            merged = shad[0].merge(shad[1:]) if len(shad) > 1 else shad[0]
            pl.add_mesh(merged, color="#0a0d11", opacity=0.45, line_width=1.6,
                        name=f"shadow_{h}")

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


def caption(pl, ensembles, args, launch, axis):
    """Caption built from the ensemble actually rendered, not typed in."""
    lines = []
    for h, (tracks, landings, blocked) in ensembles.items():
        good = landings[~blocked]
        if len(good) >= 3:
            cov = np.cov(good[:, :2].T)
            sig = float(np.sqrt(0.5 * np.trace(cov)))
        else:
            sig = float("nan")
        down = ((landings[:, 0] - launch[0]) * axis[0]
                + (landings[:, 1] - launch[1]) * axis[1])
        apex = float(np.mean([r[:, 2].max() for r in tracks]))
        lines.append(f"{HEADING_LABEL[h]}    range {np.mean(down):.0f} m    "
                     f"apex {apex:.1f} m    scatter {sig:.2f} m    "
                     f"blocked {100 * blocked.mean():.0f}%")
    # Explicit pixel positions. `position="upper_left"` clips the first
    # character against the frame edge at 1920 px, and `lower_right` runs
    # underneath the height inset.
    W, H = args.width, args.height
    pl.add_text("\n".join(lines), position=(int(0.020 * W), int(0.950 * H)),
                font_size=11, color="#dfe7f0", name="caption")
    pl.add_text(sv.view_note(args.camera, args.launch),
                position=(int(0.020 * W), int(0.912 * H)),
                font_size=10, color="#8fa0b3", name="viewnote")
    src = Path(args.field).name if args.field else "flow_fields.npz"
    swept = []
    if args.elev_range:
        swept.append(f"elev {args.elev_range[0]:g} to {args.elev_range[1]:g} deg")
    if args.bearing_range:
        swept.append(f"aim {args.bearing_range[0]:+g} to "
                     f"{args.bearing_range[1]:+g} deg")
    if args.roll_range:
        swept.append(f"bank {args.roll_range[0]:+g} to "
                     f"{args.roll_range[1]:+g} deg")
    sweep_txt = ("  |  swept uniformly: " + ", ".join(swept)) if swept else ""
    u_txt = (f"  |  Uref {args.u_ref:g} m/s" if args.u_ref is not None
             else "  |  Uref 5 m/s")
    line1 = (f"{args.n_throws} x {args.throw}  |  6-DOF x RANS x sweeping "
             f"gusts  |  {src}{u_txt}")
    line2 = (f"launch ({args.launch[0]:+.0f}, {args.launch[1]:+.0f})  |  "
             f"camera: {args.camera}{sweep_txt}")
    # Two lines, and both above the colour bar: on one line this ran under the
    # height inset and the sweep description was invisible in exactly the clips
    # whose whole point was the sweep.
    pl.add_text(line1, position=(int(0.020 * W), int(0.150 * H)),
                font_size=9, color="#8fa0b3", name="sub")
    pl.add_text(line2, position=(int(0.020 * W), int(0.122 * H)),
                font_size=9, color="#8fa0b3", name="sub2")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-throws", type=int, default=60)
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--seconds", type=float, default=12.0)
    ap.add_argument("--width", type=int, default=1920)
    ap.add_argument("--height", type=int, default=1080)
    ap.add_argument("--tube-radius", type=float, default=0.32)
    ap.add_argument("--launch", nargs=2, type=float, default=None,
                    metavar=("X", "Y"),
                    help="release point, site-frame x y in metres "
                         "(default: site_view.LAUNCH, 3 m inside the western "
                         "goal line on the centreline)")
    ap.add_argument("--u-ref", type=float, default=None,
                    help="reference wind speed at 10 m, m/s. The CFD is solved "
                         "at 5.0; other values rescale it, which "
                         "P6_Theory/wind_scaling.py measured to be legal to "
                         "0.15%% over 2.5-10 m/s")
    ap.add_argument("--elev-range", nargs=2, type=float, default=None,
                    metavar=("LO", "HI"),
                    help="sweep the VERTICAL release angle uniformly over this "
                         "range, deg, instead of drawing it from the throw's "
                         "own Gaussian")
    ap.add_argument("--bearing-range", nargs=2, type=float, default=None,
                    metavar=("LO", "HI"),
                    help="sweep the HORIZONTAL release angle uniformly over "
                         "this offset range about --bearing, deg")
    ap.add_argument("--roll-range", nargs=2, type=float, default=None,
                    metavar=("LO", "HI"),
                    help="sweep the release bank angle (hyzer/anhyzer) "
                         "uniformly over this range, deg")
    ap.add_argument("--speed-max", type=float, default=None,
                    help="upper end of the ground colour scale, m/s. Default "
                         "is 3.5 at the CFD's own 5 m/s and scales with "
                         "--u-ref, so the colour means the same fraction of "
                         "the reference wind in every clip of a speed set")
    ap.add_argument("--ground-heading", type=int, default=None)
    ap.add_argument("--headings", nargs="+", type=int, default=list(HEADINGS))
    ap.add_argument("--field", default=None)
    ap.add_argument("--name", default="ensemble_fan")
    ap.add_argument("--throw", default="right_backhand",
                    help="entry in P1_AnalyticalModel/throw_catalogue.py")
    ap.add_argument("--bearing", type=float, default=BEARING_DEG,
                    help="throw bearing, deg CCW from +x (0 = due east)")
    ap.add_argument("--camera", default="thrower", choices=sorted(sv.CAMERAS))
    ap.add_argument("--cam-back", type=float, default=None)
    ap.add_argument("--cam-height", type=float, default=None)
    ap.add_argument("--focus-frac", type=float, default=None)
    ap.add_argument("--focus-z", type=float, default=None)
    ap.add_argument("--no-panel", action="store_true",
                    help="omit the height inset")
    ap.add_argument("--shadows", action="store_true", default=True)
    ap.add_argument("--no-shadows", dest="shadows", action="store_false")
    ap.add_argument("--allow-clamped", action="store_true",
                    help="render even if the flight leaves the sampled z range")
    ap.add_argument("--outdir", default="video")
    ap.add_argument("--preview", action="store_true", help="one still, no video")
    ap.add_argument("--interactive", action="store_true")
    args = ap.parse_args()

    args.launch = (np.array([args.launch[0], args.launch[1], LAUNCH[2]])
                   if args.launch else LAUNCH)
    if args.speed_max is None:
        # Scaled, not fixed. A speed set rendered on one fixed colour range
        # would saturate every clip above the lowest wind speed and the ground
        # would carry no information in four of them.
        args.speed_max = 3.5 * ((args.u_ref / 5.0) if args.u_ref else 1.0)
    for h in args.headings:
        COLOUR.setdefault(h, (0.65, 0.78, 0.55))
    if args.ground_heading is None:
        args.ground_heading = args.headings[0]

    outdir = HERE / args.outdir
    frames = outdir / "frames"
    frames.mkdir(parents=True, exist_ok=True)

    pl, ensembles = build_scene(args)

    lands = np.vstack([np.asarray(e[1])[~np.asarray(e[2])][:, :2]
                       for e in ensembles.values()
                       if (~np.asarray(e[2])).any()])
    lands3 = np.column_stack([lands, np.zeros(len(lands))])
    pos, focus, up, va, span, mean_land = sv.camera_for(
        args.camera, args.launch, lands3,
        overrides=dict(back=args.cam_back, height=args.cam_height,
                       focus_frac=args.focus_frac, focus_z=args.focus_z))
    pl.camera_position = [tuple(pos), tuple(focus), tuple(up)]
    pl.camera.view_angle = va
    axis, _, _ = sv.flight_axis(args.launch, lands3)
    print(f"  camera '{args.camera}': ({pos[0]:+.1f}, {pos[1]:+.1f}, {pos[2]:.1f}) m, "
          f"looking down a {span:.1f} m flight axis toward "
          f"({mean_land[0]:+.1f}, {mean_land[1]:+.1f})")

    panel_h = args.headings[0]

    def draw(frac, alpha):
        add_tracks(pl, ensembles, frac, args)
        if alpha > 0:
            add_landings(pl, ensembles, alpha, args)
        caption(pl, ensembles, args, args.launch, axis)

    def shoot(path, frac):
        pl.screenshot(str(path))
        if not args.no_panel:
            hp.composite(str(path),
                         hp.render_panel(ensembles[panel_h][0], frac,
                                         args.launch, axis, panel_h,
                                         COLOUR[panel_h], u_ref=args.u_ref))

    if args.preview:
        draw(1.0, 1.0)
        still = outdir / f"{args.name}_still.png"
        shoot(still, 1.0)
        print(f"wrote {still}")
        return

    n_frames = int(args.fps * args.seconds)
    print(f"  rendering {n_frames} frames at {args.width}x{args.height} ...")
    for i in range(n_frames):
        u = i / max(1, n_frames - 1)
        frac = min(1.0, u / 0.55)
        alpha = 0.0 if u < 0.55 else min(1.0, (u - 0.55) / 0.25)
        draw(frac, alpha)
        shoot(frames / f"f{i:04d}.png", frac)
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
