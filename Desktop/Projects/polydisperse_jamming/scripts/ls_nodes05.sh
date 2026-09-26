#!/bin/bash
# PREREG_LSPILOT Amendment 4: node convergence at eps 0.5, mu 0 (400 and 2000 surface nodes).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
for nn in 400 2000; do
  python3 -u ls_pack.py --eps 0.5 --mu 0.0 --phat 1e-3 --n 100 --seed 1 --tag nodes05 --nodes $nn --max-iter 300000 \
      > ~/pj/out/lsdem/nodes05_$nn.log 2>&1 &
done
wait
cp ~/pj/out/lsdem/ls_nodes05.jsonl "$SRC/out/lsdem/" 2>/dev/null
echo "nodes05 done $(date -Is)"
