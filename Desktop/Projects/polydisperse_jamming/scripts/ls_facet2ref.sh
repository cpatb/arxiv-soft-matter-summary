#!/bin/bash
# PREREG_FACET Am. 3: N 200 sphere references, seeds 2-5, torque-free law (exact for spheres), census on.
cd ~/pj || exit 1
tr -d '\r' < /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/ls_pack.py > ~/pj/ls_pack_ref.py
source scripts/env.sh
export OMP_NUM_THREADS=1
for s in 2 3 4 5; do for m in 0.0 0.3; do
  python3 -u ls_pack_ref.py --n 200 --eps 1.0 --mu $m --seed $s --spherical-bodies 1 --phat 1e-3 --tag facet2ref \
      --nodes 1000 --max-iter 300000 --census > ~/pj/out/lsdem/facet2ref_n200_m${m}_s${s}.log 2>&1 &
done; done
wait
echo "ref done $(date -Is)"
