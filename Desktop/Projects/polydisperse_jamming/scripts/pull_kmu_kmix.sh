#!/bin/bash
# Pull A1 kmu / A2 kmix raw shards and their per-contact dumps to the Windows side (raw only).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p "$SRC/out/shape" "$SRC/out/kmer_dump"
n=0
for f in ~/pj/out/phase2_kmu_*.jsonl ~/pj/out/phase2_kmix_*.jsonl; do
  [ -e "$f" ] || continue
  cp "$f" "$SRC/out/shape/" || exit 1
  n=$((n + 1))
done
cp ~/pj/out/kmer_dump/*.npz "$SRC/out/kmer_dump/" || exit 1
echo "pulled $n shard files; kmu $(cat ~/pj/out/phase2_kmu_*.jsonl | wc -l), kmix $(cat ~/pj/out/phase2_kmix_*.jsonl | wc -l); dumps $(ls ~/pj/out/kmer_dump/*.npz | wc -l)"
