#!/usr/bin/env bash
# run_case_fast.sh <Re_d> <H_over_d> <refine> [nproc]
#
# Same case as run_case.sh, but staged into the WSL filesystem and decomposed.
# Written because the refine = 4 mesh (491 520 cells) run serially under /mnt/c
# managed 18 iterations in six minutes against a 3000-iteration cap -- sixteen
# hours -- while the machine was carrying another campaign. Two changes:
#   * the case runs in ~/bidet_gci, not on the 9p mount;
#   * simpleFoam runs on $NPROC ranks via decomposePar/reconstructPar.
# The finished case is copied back beside the others so fit_closure.py and gci.py
# see it exactly as if run_case.sh had produced it.
#
# The OpenFOAM bashrc returns nonzero on this install, so it MUST be sourced
# before `set -e` or the script dies silently with no log.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
RE="$1"; HD="$2"; REF="$3"; NPROC="${4:-4}"

cd "$HERE"
# TAKE THE TAG FROM gen_case.py, DO NOT RECONSTRUCT IT.
#
# This script used to rebuild the directory name with `printf '%.0f'`, which is
# round-half-to-EVEN, while gen_case.py builds it with Python's `int()`, which TRUNCATES.
# They agree on integers and disagree otherwise -- at refine = 1.5, printf gives "2" and
# int() gives 1, so gen_case wrote case_..._r1 while this script then cd'd into
# case_..._r2 and, at the end, `rm -rf`'d and overwrote it. On 2026-09-03 that came within
# seconds of destroying the converged refine = 2 case the whole grid study rests on; what
# it did destroy was another case's metadata, recovered from git.
#
# gen_case.py already prints "wrote <dir>  (Re=..., refine=...)" as its last line. Parse
# that. One source of truth for the name, and it is the code that creates it.
GEN_OUT="$(python3 gen_case.py "$RE" "$HD" "$REF")"
echo "$GEN_OUT"
TAG="$(printf '%s\n' "$GEN_OUT" | awk '/^wrote /{print $2}' | tail -1)"
if [ -z "$TAG" ] || [ ! -d "$HERE/$TAG" ]; then
  echo "could not read the case directory from gen_case.py output; refusing to guess" >&2
  exit 1
fi

WORK="$HOME/bidet_gci/$TAG"
rm -rf "$WORK"
mkdir -p "$HOME/bidet_gci"
cp -r "$HERE/$TAG" "$WORK"
cd "$WORK"

cat > system/decomposeParDict <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
method          simple;
coeffs { n ($NPROC 1 1); }
EOF

blockMesh > log.blockMesh 2>&1
checkMesh > log.checkMesh 2>&1 || true
decomposePar > log.decomposePar 2>&1
mpirun -np "$NPROC" simpleFoam -parallel > log.simpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1

CONV=$(grep -c "SIMPLE solution converged" log.simpleFoam || true)
LASTT=$(grep -E "^Time = " log.simpleFoam | tail -1)
echo "[$TAG] converged_flag=$CONV  $LASTT  nproc=$NPROC"

# bring the finished case back beside its siblings, dropping the processor dirs
rm -rf processor*
rm -rf "$HERE/$TAG"
cp -r "$WORK" "$HERE/$TAG"

cd "$HERE"
python3 fit_closure.py "$TAG"
