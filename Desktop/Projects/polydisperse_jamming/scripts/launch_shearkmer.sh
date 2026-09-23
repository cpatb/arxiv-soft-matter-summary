#!/bin/bash
# shearkmer -- PREREG_SHEARKMER_2026-09-22.md. 18 cells. Foreground; hold wsl.exe open.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out ~/pj/scripts ~/pj/out/shearkmer_dump
for f in shear.py phase9_constp.py packing.py measure.py psd.py shape.py; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
pids=()
for w in 0.0 0.5 1.0; do
  for i in 0 1 2 3 4 5; do
    stdbuf -oL python3 -u phase9_constp.py --cvs 0.0 --mus 0.0,0.3 --p-stars 3e-4 \
        --phi-starts 0.55 --gains 0.15 --seeds 1,2,3 --n 1000 --tag shearkmer_w$w \
        --shear-steps 102000 --doublet-w $w --shard $i --nshards 6 \
        --dump-dir ~/pj/out/shearkmer_dump > ~/pj/out/shearkmer_w${w}_$i.log 2>&1 &
    pids+=($!)
  done
done
date -u +"start %Y-%m-%dT%H:%M:%SZ"
for p in "${pids[@]}"; do wait "$p"; done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
mkdir -p "$SRC/out/shearkmer_dump"
cp ~/pj/out/phase9_shearkmer_*.jsonl "$SRC/out/"
cp ~/pj/out/shearkmer_dump/*.npz "$SRC/out/shearkmer_dump/"
echo "synced $(cat ~/pj/out/phase9_shearkmer_*.jsonl | wc -l) records"
