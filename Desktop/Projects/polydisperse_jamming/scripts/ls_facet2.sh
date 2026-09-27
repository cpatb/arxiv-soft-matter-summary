#!/bin/bash
# PREREG_FACET_2026-09-27 Amendment 1: rotations live (sphericalBodies False). 50 cells, 20 workers, priority queue.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
md5sum "$SRC/ls_pack.py" ~/pj/ls_pack.py
grep -n "^SPHERICAL_BODIES" ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
Q=~/pj/out/lsdem/facet2_queue.txt
: > "$Q"
for s in 1 2 3; do for e in 1.0 0.9 0.8 0.7 0.5; do for m in 0.0 0.3; do echo "100 $e $m $s 20000" >> "$Q"; done; done; done
for s in 1 2 3 4 5; do for e in 1.0 0.5; do for m in 0.0 0.3; do echo "200 $e $m $s 0" >> "$Q"; done; done; done
wc -l "$Q"
one() {
  set -- $1
  echo "start N$1 e$2 m$3 s$4 $(date -Is)"
  python3 -u ls_pack.py --n "$1" --eps "$2" --mu "$3" --seed "$4" --post-iters "$5" --phat 1e-3 --tag facet2 \
      --nodes 1000 --max-iter 300000 --census > ~/pj/out/lsdem/facet2_n$1_e$2_m$3_s$4.log 2>&1
  echo "end N$1 e$2 m$3 s$4 rc=$? $(date -Is)"
}
export -f one
tr '\n' '\0' < "$Q" | xargs -0 -P 20 -I{} bash -c 'one "{}"'
echo "facet2 done $(date -Is)"
