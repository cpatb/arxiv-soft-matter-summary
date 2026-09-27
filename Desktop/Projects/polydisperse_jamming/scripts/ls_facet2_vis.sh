#!/bin/bash
# Re-run the four frictionless seed-1 N-100 facet2 cells IDENTICALLY (repo ls_pack.py = queue copy + the LS_PROGRESS diagnostic only) with LS_PROGRESS
# logging, after killing their silent originals, so convergence toward the stop rule is visible (2026-09-27 ~12:30).
cd ~/pj || exit 1
tr -d "" < /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/ls_pack.py > ~/pj/ls_pack_live.py
source scripts/env.sh
export OMP_NUM_THREADS=1
for e in 0.9 0.8 0.7 0.5; do
  ps -eo pid,comm,args | awk -v e="--eps $e " '$2=="python3" && /facet2/ && /--n 100 / && /--seed 1 / && /--mu 0.0 / && index($0,e) {print $1}' | xargs -r kill
done
sleep 2
for e in 0.9 0.8 0.7 0.5; do
  LS_PROGRESS=2000 python3 -u ls_pack_live.py --n 100 --eps $e --mu 0.0 --seed 1 --post-iters 20000 --phat 1e-3 --tag facet2 \
      --nodes 1000 --max-iter 300000 --census > ~/pj/out/lsdem/facet2vis_n100_e${e}_m0.0_s1.log 2>&1 &
done
wait
echo "vis done $(date -Is)"
