#!/bin/bash
# shapechi -- PREREG_SHAPECHI_2026-09-22.md
#
# Is SHAPE a window-knob (moves chi = phim/phi0) or a position-knob (moves x = phi/phi0)?
# The entire 126-record shape campaign ran at mu_p = 0.3 ONLY, so it measures neither
# phi0 nor phim and cannot form chi at all. This runs the SAME shape grid at the two
# friction endpoints.
#
#   lognormal, CV = 0, N = 1700, alpha in {1.0, 2.0}, mu_p in {0.0, 0.3, 3.0},
#   seeds 1-3, 4 rungs = 72 records, PJ_MAX_ITER = 1000000.
#
# N = 1700 not 5400: alpha = 2 costs a median 5497 s at N = 5400 and 233 s at N = 1700,
# and chi is a RATIO AT FIXED N. The mu_p = 0.3 arm doubles as the wiring gate -- three
# existing N = 1700 records must come back bit-identical.
#
# THE WAIT IS LOAD-BEARING: this script must stay in the FOREGROUND of a held-open
# background call or WSL kills the worker group. Invoke BY FILE PATH from PowerShell.
set -u
NSHARDS="${1:-12}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000

echo "=== shapechi ==="
date -u +"start %Y-%m-%dT%H:%M:%SZ"
pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/shapechi_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase2_shape.py \
      --family lognormal --cvs 0.0 --n 1700 \
      --alphas 1.0,2.0 --mus 0.0,0.3,3.0 --seeds 1,2,3 \
      --tag shapechi --shard "$i" --nshards "$NSHARDS" > "$LOG" 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} workers"
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
echo "=== shapechi finished, $fail nonzero exits ==="

mkdir -p "$SRC/out/shape"
cp ~/pj/out/phase2*shapechi*.jsonl "$SRC/out/shape/" 2>/dev/null
cp ~/pj/out/shape/phase2*shapechi*.jsonl "$SRC/out/shape/" 2>/dev/null
cp ~/pj/out/shapechi_w*.log "$SRC/out/" 2>/dev/null
echo "=== synced ==="
cat "$SRC"/out/shape/phase2*shapechi*.jsonl 2>/dev/null | wc -l
