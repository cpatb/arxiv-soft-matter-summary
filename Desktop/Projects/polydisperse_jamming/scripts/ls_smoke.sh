#!/bin/bash
# ls_pack shakeout (not graded): does the level-set path run in a periodic PeriTriax cell? off-grid seed 99.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
timeout 1500 python3 -u ls_pack.py --eps 1.0 --mu 0.3 --phat 1e-3 --n 60 --seed 99 --tag smoke --max-iter 60000 2>&1 | tail -n 5 &
timeout 1500 python3 -u ls_pack.py --eps 0.5 --mu 0.3 --phat 1e-3 --n 60 --seed 99 --tag smoke --max-iter 60000 2>&1 | tail -n 5 &
wait
mkdir -p "$SRC/out/lsdem"; cp ~/pj/out/lsdem/ls_smoke.jsonl "$SRC/out/lsdem/" 2>/dev/null
echo smoke done
