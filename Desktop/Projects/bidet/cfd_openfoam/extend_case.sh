#!/usr/bin/env bash
# extend_case.sh <TAG> [nproc] [endTime] -- restart a capped case and let it reach
# residualControl.
#
# Generalises extend_r2.sh, which was hardcoded to the refine = 2 case. The need
# recurred immediately: the refine = 1.414 mesh also stopped at gen_case.py's default
# 3000-iteration endTime without converging, and a grid that has not met its own
# iterative criterion cannot be used to estimate DISCRETISATION error -- that is the whole
# finding of the `coarse` vs `coarse_conv` comparison in gci.py.
#
# Writes <TAG>_ext beside <TAG>. Restarting from the written latest time rather than from
# scratch is the cheap part: the fields are already close, so the residuals only have to
# fall the last half-decade.
#
# The OpenFOAM bashrc returns nonzero on this install, so it MUST be sourced before
# `set -e` or the script dies silently with no log.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
TAG="$1"
NPROC="${2:-4}"
ENDTIME="${3:-20000}"
SRC="$HERE/$TAG"
OUTTAG="${TAG}_ext"
WORK="$HOME/bidet_gci/$OUTTAG"

[ -d "$SRC" ] || { echo "no such case: $SRC"; exit 1; }

rm -rf "$WORK"
mkdir -p "$HOME/bidet_gci"
cp -r "$SRC" "$WORK"
cd "$WORK"
rm -rf processor* log.*

sed -i "s/^startFrom .*/startFrom       latestTime;/" system/controlDict
sed -i "s/^endTime .*/endTime         $ENDTIME;/" system/controlDict
grep -E "^(startFrom|endTime)" system/controlDict

cat > system/decomposeParDict <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
method          simple;
coeffs { n ($NPROC 1 1); }
EOF

decomposePar -latestTime > log.decomposePar 2>&1
mpirun -np "$NPROC" simpleFoam -parallel > log.simpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1

CONV=$(grep -c "SIMPLE solution converged" log.simpleFoam || true)
LASTT=$(grep -E "^Time = " log.simpleFoam | tail -1)
echo "[$OUTTAG] converged_flag=$CONV  $LASTT  nproc=$NPROC"

rm -rf processor*
rm -rf "$HERE/$OUTTAG"
cp -r "$WORK" "$HERE/$OUTTAG"

cd "$HERE"
echo "--- refit at the extended time ---"
python3 fit_closure.py "$OUTTAG"
echo "--- the capped value it replaces ---"
python3 fit_closure.py "$TAG"
