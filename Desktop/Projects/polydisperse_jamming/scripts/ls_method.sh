#!/bin/bash
# 2026-09-27 method check for the level-set frictionless edge (after Jiao et al. 2010 contradicted the pilot's phi_0 drop):
# frictionless, N 100, 1000 nodes, p 1e-3, then 40000 post-convergence iterations under the same stress control.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
md5sum "$SRC/ls_pack.py" ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
Q=~/pj/out/lsdem/method_queue.txt
: > "$Q"
for s in 1 2 3; do for e in 0.9 0.8; do echo "$e $s" >> "$Q"; done; done
for s in 1 2; do for e in 1.0 0.7 0.5; do echo "$e $s" >> "$Q"; done; done
echo "0.35 1" >> "$Q"
one() {
  set -- $1
  echo "start $1 $2 $(date -Is)"
  python3 -u ls_pack.py --eps "$1" --mu 0.0 --phat 1e-3 --n 100 --seed "$2" --tag method --nodes 1000 \
      --max-iter 300000 --post-iters 40000 --census > ~/pj/out/lsdem/method_e$1_s$2.log 2>&1
  echo "end $1 $2 rc=$? $(date -Is)"
}
export -f one
tr '\n' '\0' < "$Q" | xargs -0 -P 13 -I{} bash -c 'one "{}"'
echo "method done $(date -Is)"
