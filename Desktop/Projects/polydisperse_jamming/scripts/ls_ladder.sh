#!/bin/bash
# PREREG_LSPILOT Amendment 2: convex rounded-cube ladder at 1000 surface nodes (7 cells).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
for cell in "1.0 0.0" "0.7 0.0" "0.7 0.3" "0.5 0.0" "0.5 0.3" "0.35 0.0" "0.35 0.3"; do
  set -- $cell
  python3 -u ls_pack.py --eps $1 --mu $2 --phat 1e-3 --n 100 --seed 1 --tag ladder --nodes 1000 --max-iter 300000 \
      > ~/pj/out/lsdem/ladder_e$1_m$2.log 2>&1 &
done
wait
cp ~/pj/out/lsdem/ls_ladder.jsonl ~/pj/out/lsdem/ls_nodes.jsonl ~/pj/out/lsdem/ls_pilot.jsonl "$SRC/out/lsdem/" 2>/dev/null
echo "ladder done $(date -Is)"
