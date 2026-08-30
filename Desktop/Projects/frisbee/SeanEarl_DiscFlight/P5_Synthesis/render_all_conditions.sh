#!/bin/bash
# Render one clip per wind condition, plus the two-heading comparison, all with the
# same fixed behind-the-thrower camera so they can be watched back to back.
#
# The camera is derived per clip from the ensemble actually flown, so each condition is
# framed down its own mean flight axis rather than a shared nominal one.
set -e
cd "$(dirname "$0")"

N=${1:-70}
SEC=${2:-12}
W=${3:-1920}
H=${4:-1080}

render () {   # render <name> <headings...>
    name=$1; shift
    echo "############ $name  (headings: $*) ############"
    rm -rf video/frames
    python render_ensemble_fan.py --headings "$@" --name "$name" \
        --n-throws "$N" --seconds "$SEC" --width "$W" --height "$H"
}

render wind_090_easterly   90
render wind_210_ssw        210
render wind_270_westerly   270
render ensemble_fan        210 270

rm -rf video/frames
echo "ALL CLIPS DONE"
ls -la video/*.mp4
