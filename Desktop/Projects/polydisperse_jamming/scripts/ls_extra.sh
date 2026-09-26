#!/bin/bash
# PREREG_LSPILOT Amendment 5: extra seeds 4-5, eps 1.0/0.5 x mu 0/0.3, 1000 nodes.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
for s in 4 5; do for e in 1.0 0.5; do for m in 0.0 0.3; do
  python3 -u ls_pack.py --eps $e --mu $m --phat 1e-3 --n 100 --seed $s --tag extra --nodes 1000 --max-iter 300000 \
      > ~/pj/out/lsdem/extra_e${e}_m${m}_s$s.log 2>&1 &
done; done; done
wait
cp ~/pj/out/lsdem/ls_extra.jsonl "$SRC/out/lsdem/" 2>/dev/null
echo "extra done $(date -Is)"
