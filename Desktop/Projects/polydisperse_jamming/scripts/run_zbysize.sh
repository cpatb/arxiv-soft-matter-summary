#!/bin/bash
# One per-grain dump cell (zbysize_dump.py), foreground; pulls the npz back to the repo. 2026-09-24.
set -u
cd ~/pj || exit 1
source scripts/env.sh >/dev/null 2>&1
export OMP_NUM_THREADS=1 PJ_MAX_ITER=1000000
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
sed 's/\r$//' $SRC/zbysize_dump.py > ~/pj/zbysize_dump.py
for spec in "$@"; do
  set -- $spec
  python3 -u zbysize_dump.py --cv $1 --mu $2 --seed $3 --phat $4 > out/zbysize_$1_$2_$3_$4.log 2>&1
  echo "rc=$? $spec"; tail -1 out/zbysize_$1_$2_$3_$4.log
done
mkdir -p $SRC/out/zbysize && cp -p out/zbysize/*.npz $SRC/out/zbysize/
