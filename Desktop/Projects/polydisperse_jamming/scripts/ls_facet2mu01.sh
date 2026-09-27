#!/bin/bash
# PREREG_FACET Am. 4: mu_p 0.1 low-friction edge, N 100, seeds 1-3; eps 1 with the torque-free law. Also stops the
# unlogged frictionless superball facet2 cells (the watcher kills queued ones as they start).
cd ~/pj || exit 1
tr -d '\r' < /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/ls_pack.py > ~/pj/ls_pack_mu01.py
source scripts/env.sh
export OMP_NUM_THREADS=1
ps -eo pid,comm,args | awk '$2=="python3" && /ls_pack.py / && /facet2 / && /--mu 0.0 / {print $1}' | xargs -r kill
for s in 1 2 3; do
  python3 -u ls_pack_mu01.py --n 100 --eps 1.0 --mu 0.1 --seed $s --spherical-bodies 1 --phat 1e-3 --tag facet2mu01 \
      --nodes 1000 --max-iter 300000 --census > ~/pj/out/lsdem/facet2mu01_e1.0_s$s.log 2>&1 &
  for e in 0.9 0.8 0.7 0.5; do
    python3 -u ls_pack_mu01.py --n 100 --eps $e --mu 0.1 --seed $s --phat 1e-3 --tag facet2mu01 \
        --nodes 1000 --max-iter 300000 --census > ~/pj/out/lsdem/facet2mu01_e${e}_s$s.log 2>&1 &
  done
done
wait
echo "mu01 done $(date -Is)"
