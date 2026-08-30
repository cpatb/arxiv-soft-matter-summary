"""Site frame, pitch geometry and camera placement for the Sean Earl renders.

Kept separate from the renderer so the geometry can be asserted in tests
without importing pyvista.

THE SITE FRAME (verified, not assumed)
--------------------------------------
`P2_Geometry/buildings_to_stl.py` fixes it: metres, **x East, y North, z up**,
origin at the Sean Earl Field centroid. The lake is east, beyond about
x = +150 m. Campus Towers -- the 45 m block the paper's headline is about --
has a footprint starting 9 m east of the field's eastern edge and spanning its
full width.

The sign convention was checked against the CFD rather than trusted. Taking the
box-mean horizontal velocity at z = 6 m over all non-solid cells:

    theta_wind = 90 deg  (easterly)   blows toward 302 deg true   (free stream 270)
    theta_wind = 210 deg (SSW)        blows toward   4 deg true   (free stream  30)
    theta_wind = 270 deg (westerly)   blows toward 114 deg true   (free stream  90)

Each is within ~30 deg of the free-stream direction it should carry, the
residual being wake deflection inside the canopy, and each has the right SIGN.
`theta_wind` is meteorological -- the direction the wind blows FROM.

THROW BEARING is a different convention and they are easy to confuse. It is
`disc_6dof.build_initial_conditions`'s `heading_deg`: the azimuth of the release
velocity measured from +x COUNTERCLOCKWISE, i.e. standard mathematical
convention. So:

    bearing 0 deg  = due EAST, toward Campus Towers and the lake
    bearing 90 deg = due NORTH

The production Monte Carlo sweep throws at bearing 0. `EAST` below is therefore
both "toward the lake" and "what the paper actually simulated".
"""

from __future__ import annotations

import numpy as np

# Playing surface, from the OSM footprint (78 x 25 m, long axis E-W). This is
# the same box `make_results_tex.py` and `throwability.py` average over; the
# regulation-pitch alternative is open question 2 on the status card and does
# not change any direction-level conclusion.
PITCH_X = (-31.0, 47.0)
PITCH_Y = (-10.0, 15.0)

EAST = 0.0        # throw bearing, degrees CCW from +x -- toward the lake
NORTH = 90.0

# Release point: on the pitch centreline, 3 m inside the western goal line, so
# the thrower is standing at the west end throwing the length of the field east.
# 1.5 m release height is `trajectories`' own default (roughly shoulder height).
# Asserted clear of every building in `assert_launch_clear`.
LAUNCH = np.array([PITCH_X[0] + 3.0,
                   0.5 * (PITCH_Y[0] + PITCH_Y[1]),
                   1.5])

# Camera presets. `back`/`height` are metres behind and above the release point
# along the reversed flight axis; `focus_frac` is how far down the flight the
# camera looks; `up` is the view-up vector.
CAMERAS = {
    # The view asked for: standing behind the thrower, looking east down the
    # field toward the lake. Eye height is raised above a real 1.7 m because at
    # true eye height a 3.8 m apex projects onto the horizon and the flight is
    # unreadable; 9 m is a lighting-rig platform, not a drone, and it keeps the
    # buildings in their true angular relationship to the flight.
    "thrower":  dict(back=22.0, height=9.0,  focus_frac=0.62, focus_z=3.0,
                     view_angle=36.0),
    # Same axis, further back and higher: the whole ensemble at once.
    "wide":     dict(back=38.0, height=25.0, focus_frac=0.55, focus_z=1.5,
                     view_angle=30.0),
    # From the south sideline looking north, perpendicular to the throw. This
    # is the view in which HEIGHT is legible -- from behind the thrower it is
    # not, because altitude projects onto the vertical axis of the frame
    # together with downrange distance.
    # `back` is negative here: the camera moves FORWARD of the release point so
    # it sits abeam the middle of the flight rather than abeam the thrower.
    # 16 m up because the south touchline is built up: at 9 m the camera sat
    # behind a block that occluded a third of the frame.
    "sideline": dict(back=-13.0, height=16.0, focus_frac=0.5, focus_z=5.0,
                     view_angle=30.0, offset_cross=-52.0),
}


def flight_axis(launch, landings):
    """Unit vector and span from the launch point to the MEAN landing.

    Measured from the ensemble actually flown rather than taken from the
    nominal bearing, because the wind deflects the flight by several degrees
    and by a different amount per heading. A camera placed on the nominal
    bearing is not behind the thrower for the throw the wind actually produced.
    """
    lands = np.asarray(landings, dtype=float)[:, :2]
    lands = lands[np.isfinite(lands).all(axis=1)]
    if len(lands) == 0:
        raise ValueError("no finite landings to measure the flight axis from")
    mean_land = lands.mean(axis=0)
    axis = mean_land - np.asarray(launch, dtype=float)[:2]
    span = float(np.linalg.norm(axis))
    return axis / max(span, 1e-9), span, mean_land


def camera_for(preset, launch, landings, overrides=None):
    """Return (position, focal_point, view_up, view_angle) for a preset."""
    cfg = dict(CAMERAS[preset])
    cfg.update({k: v for k, v in (overrides or {}).items() if v is not None})
    launch = np.asarray(launch, dtype=float)
    axis, span, mean_land = flight_axis(launch, landings)
    cross = np.array([-axis[1], axis[0]])       # 90 deg CCW: thrower's left

    focus = np.array([launch[0] + axis[0] * span * cfg["focus_frac"],
                      launch[1] + axis[1] * span * cfg["focus_frac"],
                      cfg["focus_z"]])
    pos = np.array([launch[0] - axis[0] * cfg["back"],
                    launch[1] - axis[1] * cfg["back"],
                    cfg["height"]])
    if cfg.get("offset_cross"):
        pos[:2] += cross * cfg["offset_cross"]
    return pos, focus, np.array([0.0, 0.0, 1.0]), cfg["view_angle"], span, mean_land


def assert_launch_clear(field_npz, launch=LAUNCH, headings=(90, 210, 270)):
    """Raise if the release point sits inside a building in any wind case.

    The solid mask is per-heading only because each case is meshed from its own
    rotated STL; the buildings are the same buildings, so a launch point clear
    in one heading and solid in another would mean a geometry bug, and this
    would catch it.
    """
    F = np.load(field_npz, allow_pickle=False)
    x, y = F["x"], F["y"]
    i = int(np.argmin(np.abs(x - launch[0])))
    j = int(np.argmin(np.abs(y - launch[1])))
    bad = []
    for h in headings:
        key = f"dir_{int(h):03d}__solid"
        if key in F.files and bool(F[key][i, j]):
            bad.append(h)
    if bad:
        raise ValueError(
            f"launch point {launch[:2]} is inside a building for heading(s) "
            f"{bad} of {field_npz}")
    return True


VIEW_NOTE = {
    "thrower":  "west end of Sean Earl Field  |  looking EAST toward the lake",
    "wide":     "behind the west end  |  looking EAST toward the lake",
    "sideline": "south touchline  |  looking NORTH across the flight, "
                "so height reads directly",
}
