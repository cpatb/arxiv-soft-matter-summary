#!/bin/bash
# kchain -- PREREG_KMER Amendment 1. 36 cells (straight-chain trimers). Foreground; hold wsl.exe open from PowerShell.
# Restartable: phase2_shape skips keys already written.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out ~/pj/out/kmer_dump ~/pj/scripts
for f in packing.py phase2_shape.py phase2_phiJ.py shape.py psd.py measure.py; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
NW=${1:-14}
KM="0:0.5:c;0:1.0:c"
pids=()
for i in $(seq 0 $((NW - 1))); do
  stdbuf -oL python3 -u phase2_shape.py --tag kchain --alphas 1.0 --mus 0.0,0.3,3.0 \
      --cvs 0.0 --seeds 1,2,3 --phats 1e-3,3e-5 --n 1700 --kmers "$KM" \
      --dump-dir ~/pj/out/kmer_dump --shard "$i" --nshards "$NW" \
      > ~/pj/out/kchain_w$(printf %02d "$i").log 2>&1 &
  pids+=($!)
done
date -u +"start %Y-%m-%dT%H:%M:%SZ"
for p in "${pids[@]}"; do wait "$p"; done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
mkdir -p "$SRC/out/shape" "$SRC/out/kmer_dump"
cp ~/pj/out/phase2_kchain_*.jsonl "$SRC/out/shape/"
cp ~/pj/out/kmer_dump/*.npz "$SRC/out/kmer_dump/"
echo "synced $(cat ~/pj/out/phase2_kchain_*.jsonl | wc -l) records"
