#!/bin/bash
# PREREG_FACET Am. 5: node-refinement control of the noise confound (mu 0.3, seed 1, N 100, rotations live, logged).
cd ~/pj || exit 1
tr -d '\r' < /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/ls_pack.py > ~/pj/ls_pack_nodes.py
source scripts/env.sh
export OMP_NUM_THREADS=1
run() { LS_PROGRESS=1000 python3 -u ls_pack_nodes.py --n 100 --eps $1 --mu 0.3 --seed 1 --nodes $2 --spherical-bodies 0 \
          --phat 1e-3 --tag facet2nodes --max-iter 300000 --census > ~/pj/out/lsdem/facet2nodes_e$1_n$2.log 2>&1; }
run 0.5 2000 & run 0.5 4000 & run 0.9 4000 & run 1.0 4000 &
wait
echo "nodes done $(date -Is)"
