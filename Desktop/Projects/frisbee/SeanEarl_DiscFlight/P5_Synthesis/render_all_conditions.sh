#!/bin/bash
# Render the clip set: one behind-the-thrower view per wind condition, the
# two-heading comparison, and one sideline view where height reads directly.
#
# THE VIEW is fixed in `site_view.py`: the thrower stands on the pitch
# centreline 3 m inside the WESTERN goal line and throws EAST, the length of the
# field, toward Campus Towers and the lake. The camera sits behind them and does
# not move, so nothing that moves in the picture is the camera. It is derived
# per clip from the ensemble actually flown, so each condition is framed down
# its own mean flight axis rather than a shared nominal one.
#
# THE FIELD. `vertical_fields.npz` is used, not the production `flow_fields.npz`.
# Both are the same CFD solution; the production sampling stops at z = 6 m and
# part of even a 27 m backhand ensemble goes above that (apex 6.6 m over 30
# throws), which would fly those throws on clamped wind. The vertical sampling
# reaches z = 100 m at 2 m horizontal spacing instead of 1 m -- coarser in plan,
# but the flight never leaves it. `render_ensemble_fan.py` refuses to render a
# clamped flight unless told to, so this is enforced rather than remembered.
set -e
cd "$(dirname "$0")"

N=${1:-70}
SEC=${2:-12}
W=${3:-1920}
H=${4:-1080}
FIELD=${FIELD:-../P3_CFD/vertical_fields.npz}
THROW=${THROW:-right_backhand}

render () {   # render <name> <camera> <headings...>
    name=$1; cam=$2; shift 2
    echo "############ $name  camera=$cam  headings: $*  ############"
    rm -rf video/frames
    python render_ensemble_fan.py --headings "$@" --name "$name" \
        --camera "$cam" --throw "$THROW" --field "$FIELD" \
        --n-throws "$N" --seconds "$SEC" --width "$W" --height "$H"
}

render wind_090_easterly   thrower   90
render wind_210_ssw        thrower   210
render wind_270_westerly   thrower   270
render ensemble_fan        thrower   210 270
render sideline_090        sideline  90

rm -rf video/frames
echo "ALL CLIPS DONE"
ls -la video/*.mp4
