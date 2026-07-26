#!/bin/bash
# Mesh and solve ONE case. Usage: bash runCase.sh cases/dir_090 [nproc]
#
# WHY MESHING IS SERIAL
# ---------------------
# The tempting pipeline is decomposePar -> parallel snappyHexMesh -> solve, but
# snappyHexMesh -overwrite changes the mesh after decomposition, so every
# processor's 0/ directory is stale. Re-seeding it by copying 0.orig fails in
# turn, because the copied dictionaries carry no entries for the processor
# boundary patches ("Cannot find patchField entry for procBoundary4to0").
#
# Meshing in serial and decomposing ONCE afterwards makes decomposePar generate
# the processor patch fields itself, which is always consistent. The solve --
# the part that actually dominates wall-clock -- still runs in parallel. At the
# cell counts used here serial meshing costs a few minutes per case.
set -e
CASE=${1:?usage: runCase.sh <caseDir> [nproc]}
NP=${2:-11}
cd "$CASE"
T0=$SECONDS

echo "=== [$CASE] blockMesh ==="
blockMesh > log.blockMesh 2>&1

if [ -f EMPTY_FETCH ]; then
    echo "=== [$CASE] empty fetch: no buildings, skipping snappyHexMesh ==="
else
    echo "=== [$CASE] surfaceFeatureExtract ==="
    surfaceFeatureExtract > log.sfe 2>&1
    echo "=== [$CASE] snappyHexMesh (serial) ==="
    snappyHexMesh -overwrite > log.snappy 2>&1
    grep -E "^Snapped mesh : cells" log.snappy | tail -1
fi

rm -rf 0 && cp -r 0.orig 0

echo "=== [$CASE] checkMesh ==="
checkMesh > log.checkMesh 2>&1 || true
grep -E "cells:|non-orthogonality Max|Max skewness" log.checkMesh | head -4

echo "=== [$CASE] decomposePar ($NP) ==="
decomposePar -force > log.decomposePar 2>&1

echo "=== [$CASE] simpleFoam ==="
mpirun -np "$NP" simpleFoam -parallel > log.simpleFoam 2>&1 || true
grep -E "^Time = " log.simpleFoam | tail -1
grep -E "SIMPLE solution converged" log.simpleFoam | tail -1 || \
    echo "  (did not hit residualControl; check convergence)"

echo "=== [$CASE] reconstructPar ==="
reconstructPar -latestTime > log.reconstructPar 2>&1
echo "=== [$CASE] done in $((SECONDS-T0)) s ==="
