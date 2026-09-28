#!/usr/bin/env bash
# run_so_lane.sh <nproc> <scheme:fo|so> <refine> [<scheme> <refine> ...]
#
# One compute lane for PREREG_so_triplet_2026-09-27.md: runs the listed cases one after
# another on <nproc> ranks each, staged in ~/bidet_so (NOT ~/lfdem or ~/pj, which belong
# to another session), niced, then copies each finished case back beside its siblings and
# runs fit_closure.py on it. Re_d = 1e4, H/d = 6, endTime 30000 for every case.
#
# The case name is PARSED from gen_case.py's output, never rebuilt here (see run_case.sh).
# The OpenFOAM bashrc returns nonzero, so it is sourced BEFORE set -e.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
NPROC="$1"; shift
mkdir -p "$HOME/bidet_so"

while [ $# -ge 2 ]; do
  SCH="$1"; REF="$2"; shift 2
  cd "$HERE"
  GEN_OUT="$(python3 gen_case.py 10000 6 "$REF" --scheme="$SCH" --end=30000)"
  echo "$GEN_OUT"
  TAG="$(printf '%s\n' "$GEN_OUT" | awk '/^wrote /{print $2}' | tail -1)"
  TAG="$(basename "$TAG")"
  if [ -z "$TAG" ] || [ ! -d "$HERE/$TAG" ]; then
    echo "could not read the case directory from gen_case.py output; refusing to guess" >&2
    exit 1
  fi
  WORK="$HOME/bidet_so/$TAG"
  rm -rf "$WORK"
  cp -r "$HERE/$TAG" "$WORK"
  cd "$WORK"
  # intermediate writes so a killed run can be restarted (purgeWrite 1 keeps the last)
  sed -i 's/^writeInterval .*/writeInterval   1000;/' system/controlDict
  echo "[$TAG] start $(date '+%F %T') nproc=$NPROC" | tee -a "$HOME/bidet_so/lane.log"
  blockMesh > log.blockMesh 2>&1
  checkMesh > log.checkMesh 2>&1 || true
  if [ "$NPROC" -gt 1 ]; then
    cat > system/decomposeParDict <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
method          simple;
coeffs { n ($NPROC 1 1); }
EOF
    decomposePar > log.decomposePar 2>&1
    nice -n 5 mpirun -np "$NPROC" simpleFoam -parallel > log.simpleFoam 2>&1 || true
    reconstructPar -latestTime > log.reconstructPar 2>&1 || true
    rm -rf processor*
  else
    nice -n 5 simpleFoam > log.simpleFoam 2>&1 || true
  fi
  CONV=$(grep -c "SIMPLE solution converged" log.simpleFoam || true)
  LASTT=$(grep -E "^Time = " log.simpleFoam | tail -1)
  echo "[$TAG] end $(date '+%F %T') converged_flag=$CONV $LASTT" | tee -a "$HOME/bidet_so/lane.log"
  rm -rf "$HERE/$TAG"
  cp -r "$WORK" "$HERE/$TAG"
  cd "$HERE"
  python3 fit_closure.py "$TAG" || echo "[$TAG] fit_closure failed" | tee -a "$HOME/bidet_so/lane.log"
done
echo "lane done $(date '+%F %T')" | tee -a "$HOME/bidet_so/lane.log"
