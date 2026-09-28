#!/usr/bin/env bash
# so_continue.sh <nproc> <case_tag> [<case_tag> ...]
#
# PREREG_so_triplet S6: continue a FINISHED case 400 iterations from its final state and
# write the result beside it as <case_tag>_cont, so the grader can compare every gate
# quantity before and after. The finished case itself is never modified.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
NPROC="$1"; shift
EXTRA=400
for TAG in "$@"; do
  SRC="$HERE/$TAG"
  [ -d "$SRC" ] || { echo "no $SRC" >&2; exit 1; }
  LAST=$(ls "$SRC" | grep -E '^[0-9]+$' | sort -n | tail -1)
  [ "$LAST" != "0" ] || { echo "$TAG has no final time" >&2; exit 1; }
  END=$((LAST + EXTRA))
  # NOT named case*: so_watch.py globs case* and would stop a continuation at once
  WORK="$HOME/bidet_so/cont_${TAG}"
  rm -rf "$WORK"; mkdir -p "$WORK"
  cp -r "$SRC/constant" "$SRC/system" "$SRC/case_info.txt" "$WORK/"
  cp -r "$SRC/$LAST" "$WORK/"
  rm -f "$WORK/$LAST/wallShearStress" "$WORK/$LAST/yPlus"
  cd "$WORK"
  sed -i "s/^startFrom .*/startFrom       latestTime;/; s/^stopAt .*/stopAt          endTime;/; s/^endTime .*/endTime         $END;/; s/^writeInterval .*/writeInterval   $EXTRA;/" system/controlDict
  echo "[${TAG}_cont] start $(date '+%F %T') from $LAST to $END nproc=$NPROC" | tee -a "$HOME/bidet_so/lane.log"
  if [ "$NPROC" -gt 1 ]; then
    cat > system/decomposeParDict <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
method          simple;
coeffs { n ($NPROC 1 1); }
EOF
    decomposePar -latestTime > log.decomposePar 2>&1
    nice -n 5 mpirun -np "$NPROC" simpleFoam -parallel > log.simpleFoam 2>&1 || true
    reconstructPar -latestTime > log.reconstructPar 2>&1 || true
    rm -rf processor*
  else
    nice -n 5 simpleFoam > log.simpleFoam 2>&1 || true
  fi
  LASTT=$(grep -E "^Time = " log.simpleFoam | tail -1)
  echo "[${TAG}_cont] end $(date '+%F %T') $LASTT" | tee -a "$HOME/bidet_so/lane.log"
  rm -rf "$HERE/${TAG}_cont"
  cp -r "$WORK" "$HERE/${TAG}_cont"
  cd "$HERE"
  python3 fit_closure.py "${TAG}_cont" || echo "[${TAG}_cont] fit_closure failed" | tee -a "$HOME/bidet_so/lane.log"
done
