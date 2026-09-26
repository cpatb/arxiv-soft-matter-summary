#!/bin/bash
# PREREG_LSPILOT Amendment 1: node-convergence test on the eps 1.0, mu 0.3 cell (seed 1, N 100, p 1e-3).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
for nn in 400 1000; do
  python3 -u ls_pack.py --eps 1.0 --mu 0.3 --phat 1e-3 --n 100 --seed 1 --tag nodes --nodes $nn --max-iter 300000 \
      > ~/pj/out/lsdem/nodes_$nn.log 2>&1 &
done
wait
cp ~/pj/out/lsdem/ls_nodes.jsonl "$SRC/out/lsdem/" 2>/dev/null
echo "nodes done $(date -Is)"
