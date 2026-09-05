#!/usr/bin/env bash
# run_case.sh <Re_d> <H_over_d> [refine]
# Generates and runs one axisymmetric impinging-jet case, then fits the closure.
source /usr/lib/openfoam/openfoam2606/etc/bashrc   # must precede 'set -e' (bashrc returns nonzero)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

RE="$1"; HD="$2"; REF="${3:-1}"

# TAKE THE TAG FROM gen_case.py, DO NOT MIRROR IT. `printf '%.0f'` is round-half-to-even
# and gen_case.py uses Python's truncating `int()`: at refine = 1.5 they give 2 and 1, so
# the mirrored name pointed at a DIFFERENT case, which this script would then overwrite.
# See run_case_fast.sh for the full account (2026-09-03).
GEN_OUT="$(python3 gen_case.py "$RE" "$HD" "$REF")"
echo "$GEN_OUT"
TAG="$(printf '%s\n' "$GEN_OUT" | awk '/^wrote /{print $2}' | tail -1)"
if [ -z "$TAG" ] || [ ! -d "$HERE/$TAG" ]; then
  echo "could not read the case directory from gen_case.py output; refusing to guess" >&2
  exit 1
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
