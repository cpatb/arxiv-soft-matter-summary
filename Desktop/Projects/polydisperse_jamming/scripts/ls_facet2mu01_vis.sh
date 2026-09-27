#!/bin/bash
# Execution-only restart (2026-09-27 ~13:20): the 12 mu 0.1 superball cells, identical, with LS_PROGRESS logging.
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
ps -eo pid,comm,args | awk '$2=="python3" && /ls_pack_mu01.py/ && !/--eps 1.0 / {print $1}' | xargs -r kill
sleep 2
for s in 1 2 3; do for e in 0.9 0.8 0.7 0.5; do
  LS_PROGRESS=1000 python3 -u ls_pack_mu01.py --n 100 --eps $e --mu 0.1 --seed $s --phat 1e-3 --tag facet2mu01 \
      --nodes 1000 --max-iter 300000 --census > ~/pj/out/lsdem/facet2mu01_e${e}_s$s.log 2>&1 &
done; done
wait
echo "mu01 vis done $(date -Is)"
