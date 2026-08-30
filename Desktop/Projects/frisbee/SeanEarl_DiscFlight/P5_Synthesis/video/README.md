# The flight clips

Five 12-second clips, 1920×1080, 30 fps, 70 throws each. Regenerate all of them
with `bash ../render_all_conditions.sh`.

## The view

The thrower stands on the pitch centreline **3 m inside the western goal line**
and throws **due east**, the length of the field, toward Campus Towers and the
lake beyond it. The camera sits **behind them on the flight axis, 22 m back and
9 m up, and does not move for the whole clip** — nothing that moves in the
picture is the camera.

The site frame is x East, y North, origin at the field centroid
(`P2_Geometry/buildings_to_stl.py`). The throw bearing is 0°, measured
counterclockwise from +x, which is `disc_6dof.build_initial_conditions`'s
convention and is also the bearing the production 97,200-throw sweep used — so
these are the throws the paper's numbers describe, seen from behind.

The wind heading is a *different* convention and the two are easy to confuse:
`theta_wind` is meteorological, the direction the wind blows **from**. The sign
was checked against the CFD rather than assumed; see `../site_view.py`.

9 m is higher than a person. At true eye height a 4 m apex projects onto the
horizon and the flight is unreadable. It is a step-ladder, not a drone, and it
keeps the buildings in their real angular relationship to the flight.

| clip | wind | camera |
|---|---|---|
| `wind_090_easterly.mp4` | 90°, off the lake — a headwind aloft | behind the thrower |
| `wind_210_ssw.mp4` | 210°, the prevailing SSW | behind the thrower |
| `wind_270_westerly.mp4` | 270°, westerly over the campus | behind the thrower |
| `ensemble_fan.mp4` | 210 and 270 overlaid | behind the thrower |
| `sideline_090.mp4` | 90°, easterly | south touchline, looking north |

## Reading height

From behind the thrower, altitude and downrange distance project onto the same
axis of the frame, so **z is the one thing the main view cannot show**. Two
things carry it instead:

* the **ground shadow** under each trajectory — the same flight projected onto
  z = 0, so the gap between a track and its shadow is the height;
* the **inset panel**, bottom right. Its left half is the flight profile
  `z` against downrange, growing with the animation. Its right half is the
  site's own `σ_u(z)`, measured over the playing surface from the same CFD
  solution, on a **shared height axis** — so the height a throw reaches can be
  read straight across into the disturbance it meets up there. Both axes start
  at zero.

That pairing is the point. `σ_u` rises by 1.22× from 2 m to 6 m in all three
wind directions, and the integral length scale grows with height too, so a
lofted disc gets a stronger disturbance *and* less averaging of it. Throwing low
wins twice. `sideline_090.mp4` shows the same ensemble from the touchline, where
height is legible directly.

## What is data and what is drawing

Everything with a physical position is computed:

* buildings are `P2_Geometry/seanEarl.stl` — 116 OSM footprints, 6,324
  triangles, extruded to tagged or type-default heights;
* every trajectory is a fresh integration of the paper's 6-DOF model with the
  `potts_tabulated` coefficients, flown through the same time-averaged RANS
  field with the same `sweeping` gust model and the same release-scatter
  distribution as the production sweep;
* the landing scatter is those trajectories' own endpoints and the ellipse is
  fitted to them;
* the ground colour is the disc-altitude horizontal wind speed at z = 2 m from
  the same solution, and the inset's `σ_u(z)` is `sqrt(2k/3)` from it.

Camera, colour and opacity are the only choices made for looks. No trajectory is
smoothed, resampled for shape, or hand-placed. Throws that land inside a
building are drawn in a separate dark red rather than deleted — "blocked" is a
real outcome of the model and hiding it would flatter the picture.

The caption is generated from the ensemble that was actually flown, not typed.
The colour bar is clipped at 3.5 m/s, with everything above it in one flat
colour, because the full box reaches 4.5 m/s in accelerated corners well off the
pitch and letting those set the scale flattens the 0–2 m/s range the disc
actually flies in.

## Two things the clips are not evidence for

**The direction ranking.** These are rendered on the production field, which is
what the manuscript reports, and that field resolves 270° as the worst of the
three headings. **That verdict does not survive mesh refinement** — see
`../../../GRID_INDEPENDENCE_2026-08-30.md`. Read the clips as "what the reported
field looks like", not as evidence that the headings differ.

**The absolute scatter.** The release distribution is
`trajectories.RELEASE_DISTRIBUTIONS`, which its own docstring calls
"practitioner-typical" and which is not calibrated against a measured player. It
gives a 7.4–7.8 m landing scatter on a 27 m throw, and that floor — not the wind
— is most of the fan you are looking at. At 5 m/s the wind accounts for 4–19 %
of the landing error; see `THEORY_WIND_SPEED_2026-08-30.md` §2.

## The field, and why it is not the production one

`vertical_fields.npz`, not `flow_fields.npz`. Both are the same CFD solution;
the production sampling stops at **z = 6 m**, and part of even a 27 m backhand
ensemble goes above that — apex 6.6 m over 30 throws — which would fly those
throws on wind clamped at its 6 m value. The vertical sampling reaches z = 100 m
at 2 m horizontal spacing instead of 1 m: coarser in plan, but the flight never
leaves it. `render_ensemble_fan.py` refuses to render a clamped flight unless
explicitly told to, so this is enforced rather than remembered.
