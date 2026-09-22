#!/bin/bash
# actshape finite-size control -- PREREG_ACTSHAPE Amendment 1. 6 cells at N = 4000.
# Foreground; the caller holds wsl.exe open. Restartable by key.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
DUMP=~/pj/out/actshape_dump
mkdir -p "$DUMP"
pids=()
for i in 0 1 2 3 4 5; do
  stdbuf -oL python3 -u phase9_constp.py \
      --cvs 0.5 --mus 0.0,0.3 --p-stars 3e-4 --phi-starts 0.55 \
      --gains 0.15 --seeds 1,2,3 --n 4000 --tag actshape_fs \
      --shear-steps 480000 --shard "$i" --nshards 6 \
      --dump-dir "$DUMP" > ~/pj/out/actshape_fs_w${i}.log 2>&1 &
  pids+=($!)
done
date -u +"start %Y-%m-%dT%H:%M:%SZ"
for p in "${pids[@]}"; do wait "$p"; done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
cp ~/pj/out/phase9_actshape_fs*.jsonl "$SRC/out/"
cp "$DUMP"/*N4000*.npz "$SRC/out/actshape_dump/"
cp ~/pj/out/actshape_fs_w*.log "$SRC/out/"
python3 scripts/w2_same_platform.py > out/actshape_w2wsl.txt && cp out/actshape_w2wsl.txt "$SRC/out/"
echo "synced: $(cat "$SRC"/out/phase9_actshape_fs*.jsonl | wc -l) records"
