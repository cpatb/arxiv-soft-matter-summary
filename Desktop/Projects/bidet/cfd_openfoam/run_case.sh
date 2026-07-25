#!/usr/bin/env bash
# run_case.sh <Re_d> <H_over_d> [refine]
# Generates and runs one axisymmetric impinging-jet case, then fits the closure.
source /usr/lib/openfoam/openfoam2606/etc/bashrc   # must precede 'set -e' (bashrc returns nonzero)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

RE="$1"; HD="$2"; REF="${3:-1}"
python3 gen_case.py "$RE" "$HD" "$REF"

# case dir name mirrors gen_case.py
if [ "$REF" = "1" ] || [ "$REF" = "1.0" ]; then
  TAG="case_Re$(printf '%.0f' "$RE")_H$(printf '%.0f' "$HD")"
else
  TAG="case_Re$(printf '%.0f' "$RE")_H$(printf '%.0f' "$HD")_r$(printf '%.0f' "$REF")"
fi
cd "$HERE/$TAG"

blockMesh   > log.blockMesh   2>&1
checkMesh   > log.checkMesh   2>&1 || true
simpleFoam  > log.simpleFoam  2>&1
CONV=$(grep -c "SIMPLE solution converged" log.simpleFoam || true)
LASTT=$(grep -E "^Time =" log.simpleFoam | tail -1)
echo "[$TAG] converged_flag=$CONV  $LASTT"

cd "$HERE"
python3 fit_closure.py "$TAG"
