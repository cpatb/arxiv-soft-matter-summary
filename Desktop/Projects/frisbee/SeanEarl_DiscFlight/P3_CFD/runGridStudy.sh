#!/bin/bash
# Grid-convergence study: mesh, solve, cycle-average and sample ONE case of the
# fieldBox refinement family, with every other setting held fixed.
#
# WHY THIS EXISTS
# ---------------
# `coarse` (fieldBox level 2) and `production` (level 3) are a clean factor-2
# pair, but two grids give a DIFFERENCE, not an ERROR. An observed order of
# convergence and a Richardson extrapolation need at least three geometrically
# similar grids. This script runs the extra members -- gridL1 (level 1) and
# gridL4 (level 4) -- through EXACTLY the pipeline the production cases went
# through, including the cycle-averaging pass, so the four solutions are
# comparable.
#
# The averaging window must match the production run (800 iterations) or the
# comparison measures the averaging, not the grid.
#
# Usage: bash runGridStudy.sh <caseDir> [nproc] [navg]
set -e

CASE=${1:?usage: runGridStudy.sh <caseDir> [nproc] [navg]}
NP=${2:-11}
NAVG=${3:-800}

source /usr/lib/openfoam/openfoam2606/etc/bashrc

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

echo "############################################################"
echo "# GRID STUDY: $CASE   nproc=$NP  navg=$NAVG"
echo "# started $(date -Is)"
echo "############################################################"

T0=$SECONDS

echo
echo ">>> stage 1/3: mesh + solve"
bash runCase.sh "$CASE" "$NP"

echo
echo ">>> stage 2/3: cycle-average over $NAVG iterations"
bash averagePass.sh "$NP" "$NAVG" "$CASE"

echo
echo ">>> stage 3/3: sample onto the standard disc-altitude grid"
TAG=$(basename "$CASE")
python3 sample_disc_altitude.py "$CASE" --out "flow_fields_${TAG}.npz"

echo
echo "############################################################"
echo "# DONE $CASE in $((SECONDS-T0)) s   $(date -Is)"
echo "# cells: $(grep -E '^    cells:' "$CASE/log.checkMesh" | head -1)"
echo "############################################################"
