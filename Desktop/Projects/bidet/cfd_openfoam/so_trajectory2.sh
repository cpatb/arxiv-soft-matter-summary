#!/usr/bin/env bash
# so_trajectory.sh <nproc> <case_tag> [<case_tag> ...]
#
# PREREG_so_triplet Amendment 4 (rule C1): continue a case 2400 iterations from its final
# state, writing every 400 (purgeWrite 0), regenerate the wall-shear sample at every written
# time with `simpleFoam -postProcess`, and copy the result back as <case_tag>_traj.
# The source case is never modified. Work dirs are named traj_* so so_watch.py (which globs
# case*) cannot stop them.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
NPROC="$1"; shift
TOTAL="${TOTAL:-2400}"  # Amendment 5: longer continuations for cases failing C1
EVERY=400
for TAG in "$@"; do
  SRC="$HERE/$TAG"
  [ -d "$SRC" ] || { echo "no $SRC" >&2; exit 1; }
  # the latest time with a FULL field write (an end time that is not a write time holds only
  # the end-of-run function-object output) -- copy of so_trajectory.sh, fixed 2026-09-28
  LAST=$(for x in $(ls "$SRC" | grep -E '^[0-9]+$' | sort -n); do [ -f "$SRC/$x/U" ] && echo $x; done | tail -1)
  [ "$LAST" != "0" ] || { echo "$TAG has no final time" >&2; exit 1; }
  END=$((LAST + TOTAL))
  WORK="$HOME/bidet_so/traj_${TAG}"
  rm -rf "$WORK"; mkdir -p "$WORK"
  cp -r "$SRC/constant" "$SRC/system" "$SRC/case_info.txt" "$WORK/"
  cp -r "$SRC/$LAST" "$WORK/"
  rm -f "$WORK/$LAST/wallShearStress" "$WORK/$LAST/yPlus"
  cd "$WORK"
  sed -i "s/^startFrom .*/startFrom       latestTime;/; s/^stopAt .*/stopAt          endTime;/; s/^endTime .*/endTime         $END;/; s/^writeInterval .*/writeInterval   $EVERY;/; s/^purgeWrite .*/purgeWrite      0;/" system/controlDict
  # residualControl must not end the run early: the trajectory needs every write
  sed -i 's/^\(\s*\)p\s\+1e-5;/\1p           1e-12;/; s/^\(\s*\)U\s\+1e-5;/\1U           1e-12;/; s/"(k|omega)" 1e-5;/"(k|omega)" 1e-12;/' system/fvSolution
  [ "$(grep -c "1e-12" system/fvSolution)" -eq 3 ] || { echo "fvSolution edit failed" >&2; exit 1; }
  echo "[${TAG}_traj] start $(date '+%F %T') from $LAST to $END every $EVERY nproc=$NPROC" | tee -a "$HOME/bidet_so/lane.log"
  if [ "$NPROC" -gt 1 ]; then
    cat > system/decomposeParDict <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
method          simple;
coeffs { n ($NPROC 1 1); }
EOF
    decomposePar -latestTime > log.decomposePar 2>&1
    nice -n 5 mpirun -np "$NPROC" simpleFoam -parallel > log.simpleFoam 2>&1 || true
    reconstructPar > log.reconstructPar 2>&1 || true
    rm -rf processor*
  else
    nice -n 5 simpleFoam > log.simpleFoam 2>&1 || true
  fi
  TIMES=$(ls | grep -E '^[0-9]+$' | sort -n | awk -v s="$LAST" '$1 > s' | paste -sd, -)
  simpleFoam -postProcess -time "$TIMES" > log.postProcess 2>&1 || true
  LASTT=$(grep -E "^Time = " log.simpleFoam | tail -1)
  echo "[${TAG}_traj] end $(date '+%F %T') $LASTT; post-processed times $TIMES" | tee -a "$HOME/bidet_so/lane.log"
  rm -rf "$HERE/${TAG}_traj"
  cp -r "$WORK" "$HERE/${TAG}_traj"
  cd "$HERE"
done
