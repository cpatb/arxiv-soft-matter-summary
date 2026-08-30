#!/bin/bash
# Three clip sets, all thrown EAST from behind the thrower.
#
#   SET Q  down the field in quarters  -- the same throw released from each
#          quarter point of the pitch, so what changes between the clips is
#          only where on the field you are standing. The last one is thrown
#          from inside the turbulent zone.
#
#   SET U  spanning wind speed         -- one launch point, one throw, at
#          Uref = 2, 5, 9 and 15 m/s. Legal because the CFD scales: see
#          P6_Theory/wind_scaling.py, which measured the rescaling to 0.15%
#          over a 4x range.
#
#   SET A  spanning release angle      -- vertical (elevation), horizontal
#          (aim) and bank swept one at a time so it is visible which release
#          angle is doing what, then all three together.
#
# EVERY CLIP USES theta_wind = 270 unless its name says otherwise. That is the
# westerly, and it is the heading whose disc-band turbulence is strongly
# non-uniform along the pitch -- sigma_u at z = 2 m runs 0.39 m/s at the west
# end to 0.77 m/s at x = +40, a factor 2.03, against 1.28x at 90 deg and a
# minimum-in-the-middle shape at 210. If you want to see the wind do something,
# this is the one.
#
# The RELEASE-ANGLE SWEEPS ARE UNIFORM, NOT GAUSSIAN, and the caption says so.
# The production release distribution is a Gaussian about one nominal angle,
# which is the right model for one thrower repeating one throw and the wrong
# one for showing the envelope of throws a player has. A uniform sweep is an
# envelope; reading it as a probability would overstate how often the extremes
# happen.
set -e
cd "$(dirname "$0")"

N=${N:-70}
SEC=${SEC:-12}
W=${W:-1920}
H=${H:-1080}
FIELD=${FIELD:-../P3_CFD/vertical_fields.npz}
HEADING=${HEADING:-270}

# Wide-but-real release envelope, used wherever a clip is meant to show the
# range of throws rather than the repeatability of one. Elevation 2-26 deg
# spans a flat drive to a high floating hammer-ish release; aim +/-12 deg is
# the width of a receiver's window at 25 m; bank -25 to +20 covers hyzer
# through anhyzer.
FAN="--elev-range 2 26 --bearing-range -12 12 --roll-range -25 20"

render () {   # render <name> <extra args...>
    name=$1; shift
    echo "############ $name ############"
    rm -rf video/frames
    python render_ensemble_fan.py --headings "$HEADING" --name "$name" \
        --field "$FIELD" --n-throws "$N" --seconds "$SEC" \
        --width "$W" --height "$H" "$@"
}

# ---- SET Q: quarters of the field, throwing east -------------------------
# The pitch runs x = -31 to +47, so the quarter points are -31, -11.5, +8,
# +27.5. Each launch is nudged 3 m inside so the thrower is on the field.
render q1_west_end        --launch -28.0  2.5 $FAN
render q2_quarter         --launch  -8.5  2.5 $FAN
render q3_midfield        --launch  11.0  2.5 $FAN
render q4_turbulent_zone  --launch  30.5  2.5 $FAN

# ---- SET U: spanning wind speed ------------------------------------------
# From the west end, so the flight crosses the whole field and the wind has
# the longest lever on it.
render u1_wind_2ms   --launch -28.0 2.5 --u-ref 2  $FAN
render u2_wind_5ms   --launch -28.0 2.5 --u-ref 5  $FAN
render u3_wind_9ms   --launch -28.0 2.5 --u-ref 9  $FAN
render u4_wind_15ms  --launch -28.0 2.5 --u-ref 15 $FAN

# ---- SET A: one release angle at a time -----------------------------------
render a1_vertical_angle   --launch -28.0 2.5 --elev-range 2 26
render a2_horizontal_angle --launch -28.0 2.5 --bearing-range -20 20
render a3_bank_angle       --launch -28.0 2.5 --roll-range -35 30
render a4_all_three        --launch -28.0 2.5 $FAN

rm -rf video/frames
echo "ALL SETS DONE"
ls -la video/*.mp4
