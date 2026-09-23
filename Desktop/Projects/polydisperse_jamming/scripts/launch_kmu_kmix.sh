#!/bin/bash
# A1 kmu + A2 kmix -- PREREG_KMIX_2026-09-23.md. 72 cells, 12 workers. Foreground.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out ~/pj/out/kmer_dump ~/pj/scripts
for f in packing.py phase2_shape.py phase2_phiJ.py shape.py psd.py measure.py; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1 PJ_MAX_ITER=1000000
pids=()
for i in 0 1 2 3 4 5; do
  stdbuf -oL python3 -u phase2_shape.py --tag kmu --alphas 1.0 --mus 0.6,1.0 --cvs 0.0 \
      --seeds 1,2,3 --phats 1e-3,3e-5 --n 1700 --kmers "0:0;0.5:0;1.0:0" \
      --dump-dir ~/pj/out/kmer_dump --shard $i --nshards 6 > ~/pj/out/kmu_w$i.log 2>&1 &
  pids+=($!)
  stdbuf -oL python3 -u phase2_shape.py --tag kmix --alphas 1.0 --mus 0.0,0.3 --cvs 0.0 \
      --seeds 1,2,3 --phats 1e-3,3e-5 --n 1700 --kmers "0.3:0.2:c;0.2:0.4:c;0.5:0.5:c" \
      --dump-dir ~/pj/out/kmer_dump --shard $i --nshards 6 > ~/pj/out/kmix_w$i.log 2>&1 &
  pids+=($!)
done
date -u +"start %Y-%m-%dT%H:%M:%SZ"
for p in "${pids[@]}"; do wait "$p"; done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
cp ~/pj/out/phase2_kmu_*.jsonl ~/pj/out/phase2_kmix_*.jsonl "$SRC/out/shape/"
cp ~/pj/out/kmer_dump/*.npz "$SRC/out/kmer_dump/"
echo "synced kmu $(cat ~/pj/out/phase2_kmu_*.jsonl | wc -l), kmix $(cat ~/pj/out/phase2_kmix_*.jsonl | wc -l)"
