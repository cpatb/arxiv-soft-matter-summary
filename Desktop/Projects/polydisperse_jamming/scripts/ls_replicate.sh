#!/bin/bash
# PREREG_LSPILOT Amendment 3: replicate seeds + N 200 finite-size check, 1000 nodes (16 cells).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
for s in 2 3; do for e in 1.0 0.7 0.5; do for m in 0.0 0.3; do
  python3 -u ls_pack.py --eps $e --mu $m --phat 1e-3 --n 100 --seed $s --tag repl --nodes 1000 --max-iter 300000 \
      > ~/pj/out/lsdem/repl_e${e}_m${m}_s$s.log 2>&1 &
done; done; done
for e in 1.0 0.5; do for m in 0.0 0.3; do
  python3 -u ls_pack.py --eps $e --mu $m --phat 1e-3 --n 200 --seed 1 --tag n200 --nodes 1000 --max-iter 300000 \
      > ~/pj/out/lsdem/n200_e${e}_m${m}.log 2>&1 &
done; done
wait
cp ~/pj/out/lsdem/ls_repl.jsonl ~/pj/out/lsdem/ls_n200.jsonl ~/pj/out/lsdem/ls_ladder.jsonl "$SRC/out/lsdem/" 2>/dev/null
echo "replicate done $(date -Is)"
