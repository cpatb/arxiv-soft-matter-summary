#!/bin/bash
# PREREG_LSPILOT_2026-09-26: 4 level-set cells + 2 FrictMat sphere references at N 100 (P0). Foreground.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
md5sum "$SRC/ls_pack.py" ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
for e in 1.0 0.5; do
  for m in 0.0 0.3; do
    python3 -u ls_pack.py --eps $e --mu $m --phat 1e-3 --n 100 --seed 1 --tag pilot --max-iter 300000 \
        > ~/pj/out/lsdem/pilot_e${e}_m${m}.log 2>&1 &
  done
done
stdbuf -oL python3 -u phase2_shape.py --tag lsref --family lognormal --n 100 --alphas 1.0 --mus 0.0,0.3 \
    --cvs 0.0 --seeds 1 --phats 1e-3 > ~/pj/out/lsdem/lsref.log 2>&1 &
wait
mkdir -p "$SRC/out/lsdem"
cp ~/pj/out/lsdem/ls_pilot.jsonl ~/pj/out/lsdem/*.log "$SRC/out/lsdem/" 2>/dev/null
cp ~/pj/out/phase2_lsref_lognormal.jsonl "$SRC/out/lsdem/" 2>/dev/null
echo "pilot done $(date -Is)"
