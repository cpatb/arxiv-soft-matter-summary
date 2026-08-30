#!/bin/bash
# Mesh, add the porous cell zone, solve, cycle-average and sample ONE windbreak
# case. Identical to runCase.sh + averagePass.sh except for the topoSet stage,
# which must run AFTER meshing (the zone is a set of cells, and snappyHexMesh
# changes which cells exist) and BEFORE decomposePar (so every processor gets
# its share of the zone).
#
# THE CHECK THAT MATTERS is right after topoSet: a cell zone that selects ZERO
# cells is not an error in OpenFOAM. The solver runs, the fvOptions source
# applies to nothing, and the case completes looking exactly like a successful
# windbreak run while modelling no windbreak at all. This script fails loudly on
# an empty zone instead.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
CASE=${1:?usage: runWindbreakCase.sh <caseDir> [nproc] [navg]}
NP=${2:-11}
NAVG=${3:-800}
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"

T0=$SECONDS
echo "############ WINDBREAK $CASE  nproc=$NP  $(date -Is) ############"
(
cd "$CASE"
echo "=== blockMesh ==="            ; blockMesh            > log.blockMesh 2>&1
echo "=== surfaceFeatureExtract ===" ; surfaceFeatureExtract > log.sfe 2>&1
echo "=== snappyHexMesh ==="        ; snappyHexMesh -overwrite > log.snappy 2>&1
grep -E "^Snapped mesh : cells" log.snappy | tail -1
rm -rf 0 && cp -r 0.orig 0
echo "=== checkMesh ==="            ; checkMesh > log.checkMesh 2>&1 || true
grep -E "cells:|non-orthogonality Max|Max skewness" log.checkMesh | head -3

echo "=== topoSet (porous zone) ==="
topoSet > log.topoSet 2>&1
# topoSet writes "cellZoneSet windbreak now size 460" -- NOT "460 cells".
# The first version of this guard grepped for the latter, found nothing, and
# failed a case whose zone was perfectly correct. A guard that cries wolf gets
# switched off, so it parses the string OpenFOAM actually emits.
NCELL=$(grep -oE "cellZoneSet windbreak now size [0-9]+" log.topoSet | tail -1 | grep -oE "[0-9]+$")
echo "    windbreak zone: ${NCELL:-unknown} cells"
if [ -z "$NCELL" ] || [ "$NCELL" -lt 50 ]; then
    echo "!! FAIL: the porous zone selected ${NCELL:-0} cells. A zone this small"
    echo "!! models no barrier and the run would look successful anyway."
    tail -40 log.topoSet
    exit 3
fi

# decomposeParDict is inherited from the production case and hardcodes 11
# subdomains. decomposePar honours the DICT and mpirun honours -np, so running
# with any other count produces 11 processor directories and then a fatal
# "specifies 11 processors but job was started with N ranks" -- after the mesh
# and the porous zone have already been built. Reconcile them here rather than
# discover it 6 minutes in.
echo "=== decomposePar ($NP) ==="
sed -i "s/^numberOfSubdomains .*/numberOfSubdomains  $NP;/" system/decomposeParDict
grep -E "^numberOfSubdomains" system/decomposeParDict
rm -rf processor*
decomposePar -force  > log.decomposePar 2>&1
echo "=== simpleFoam ==="
mpirun -np "$NP" simpleFoam -parallel > log.simpleFoam 2>&1 || true
grep -E "^Time = " log.simpleFoam | tail -1
echo "=== reconstructPar ==="       ; reconstructPar -latestTime > log.reconstructPar 2>&1
)

echo "=== cycle-average over $NAVG ==="
bash averagePass.sh "$NP" "$NAVG" "$CASE"

TAG=$(basename "$CASE")
echo "=== sample (disc band) ==="
python3 sample_disc_altitude.py "$CASE" --out "flow_fields_${TAG}.npz"
echo "=== sample (tall) ==="
python3 sample_disc_altitude.py "$CASE" --out "vertical_fields_${TAG}.npz" \
    --dxy 2.0 --z 0.5 1 1.5 2 3 4 5 6 8 10 13 16 20 25 30 35 40 45 50 55 60 70 80 100

echo "############ DONE $CASE in $((SECONDS-T0)) s  $(date -Is) ############"
