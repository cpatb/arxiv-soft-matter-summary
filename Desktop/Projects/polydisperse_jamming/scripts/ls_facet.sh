#!/bin/bash
# PREREG_FACET_2026-09-27: 38 level-set cells, N 200, 1000 nodes, p 1e-3, census on. Queue = priority order, 20 workers.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
md5sum "$SRC/ls_pack.py" ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
Q=~/pj/out/lsdem/facet_queue.txt
: > "$Q"
for s in 2 3 4 5; do for e in 1.0 0.5; do for m in 0.0 0.3; do echo "$e $m $s" >> "$Q"; done; done; done   # A
for e in 0.9 0.8; do for s in 1 2 3; do echo "$e 0.0 $s" >> "$Q"; done; done                            # B
for e in 0.9 0.8; do for s in 1 2 3; do echo "$e 0.3 $s" >> "$Q"; done; done                            # C
for m in 0.0 0.3; do for s in 1 2 3; do echo "0.7 $m $s" >> "$Q"; done; done                            # D
for e in 1.0 0.5; do for m in 0.0 0.3; do echo "$e $m 1" >> "$Q"; done; done                            # E
wc -l "$Q"
one() {
  set -- $1
  echo "start $1 $2 $3 $(date -Is)"
  python3 -u ls_pack.py --eps "$1" --mu "$2" --phat 1e-3 --n 200 --seed "$3" --tag facet --nodes 1000 \
      --max-iter 300000 --census > ~/pj/out/lsdem/facet_e$1_m$2_s$3.log 2>&1
  echo "end $1 $2 $3 rc=$? $(date -Is)"
}
export -f one
tr '\n' '\0' < "$Q" | xargs -0 -P 20 -I{} bash -c 'one "{}"'
echo "facet done $(date -Is)"
