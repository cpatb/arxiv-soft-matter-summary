#!/bin/bash
# Wait for the gridL4 solve to finish, then reconstruct, cycle-average over the SAME
# 800-iteration window production used, and sample onto the standard grid.
#
# Kept as a script rather than an inline `wsl.exe bash -lc '...'` because the nested
# quoting in an inline invocation mangles $(...) and tr arguments, and because a
# detached job launched from an inline invocation dies when the WSL session closes.
# Source BEFORE `set -e`: the OpenFOAM bashrc returns non-zero on this install, so
# sourcing it under `set -e` aborts the script instantly and silently -- the same class
# of hazard runSweep.sh documents for `set -u`.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

CASE=cases/dir_090_gridL4
NP=11
NAVG=800

echo "waiting for simpleFoam to finish in $CASE ..."
while pgrep -f "simpleFoam -parallel" > /dev/null 2>&1; do
    sleep 60
done
echo "solve finished: $(grep -E '^Time = ' $CASE/log.simpleFoam | tail -1)"

# The solve must actually have reached its endTime, not died early. An artifact that
# stops short would silently be averaged and compared as if it were converged.
LAST=$(grep -E '^Time = ' "$CASE/log.simpleFoam" | tail -1 | awk '{print $3}')
if [ "$LAST" != "1400" ]; then
    echo "ERROR: solve stopped at $LAST, expected 1400. Not averaging a short run."
    exit 1
fi
if ! grep -q "^End" "$CASE/log.simpleFoam"; then
    echo "ERROR: no clean End in log.simpleFoam; the solve did not terminate normally."
    exit 1
fi

echo "=== reconstructPar -latestTime ==="
( cd "$CASE" && reconstructPar -latestTime > log.reconstructPar 2>&1 )
echo "reconstructed: $(cd "$CASE" && ls -d [0-9]* | tr '\n' ' ')"

echo "=== averagePass ($NAVG iterations, matching production) ==="
bash averagePass.sh "$NP" "$NAVG" "$CASE"

echo "=== sample ==="
python3 sample_disc_altitude.py "$CASE" --out flow_fields_dir_090_gridL4.npz

echo "=== grid convergence, four members ==="
python3 grid_convergence.py --heading 90 --out grid_convergence.json

echo "ALL DONE $(date -Is)"
