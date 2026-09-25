#!/bin/bash
# rollchi shakeout (NOT graded data): off-grid seed 99. Times one cell per branch and checks
# that etaRoll is inert at mu_p = 0 (no tangential force -> spin cannot couple to translation).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out/rollchi_smoke
for f in packing.py phase2_shape.py phase2_phiJ.py shape.py psd.py measure.py; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
run () {  # tag mus etas n phats
  stdbuf -oL python3 -u phase2_shape.py --tag "rcsmoke_$1" --alphas 1.0 --mus "$2" --cvs 0.0 \
      --seeds 99 --phats "$5" --n "$4" --eta-rolls "$3" > ~/pj/out/rollchi_smoke/$1.log 2>&1
}
run m0 0.0 0.0,1.0 600 1e-3 &
run m03 0.3 0.1 1700 1e-3 &
run m3 3.0 1.0 1700 1e-3 &
run m03lo 0.3 0.1 1700 3e-5 &
wait
grep -h "phi=" ~/pj/out/rollchi_smoke/*.log
cp ~/pj/out/rollchi_smoke/*.log "$SRC/out/" 2>/dev/null
mkdir -p "$SRC/out/rollchi_smoke"; cp ~/pj/out/phase2_rcsmoke_* "$SRC/out/rollchi_smoke/"
echo smoke done
