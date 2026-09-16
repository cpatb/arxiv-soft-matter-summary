#!/bin/bash
# Pull the shape campaign's RAW shards back from WSL. From PowerShell, by file path:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/pull_shape.sh
#
# RAW RESULTS ONLY. `scripts/sync_back.sh` copies *.csv/*.txt as well and has already
# clobbered fifteen committed derived artifacts once (lessons.md 2026-09-08); it must
# never be used. Derived tables are RE-DERIVED on the Windows side, never copied.
#
# THE SHARDS LAND IN out/shape/, NOT out/. `out/phase2*.jsonl` is a SHARED GLOB read by
# ten scripts (analyse_phase2, analyse_finite_size, grade_capfix, measure_z_iso, ...).
# These cells carry an `alpha` field none of those readers know about, and their alpha=1
# rung is a REPRODUCTION of records already in that glob under a different key -- so
# dropping them into it would let a sphere analysis silently average a clump pack into
# its own ladder. A separate directory makes that impossible rather than merely
# unlikely, exactly as out/fabric/ does.
set -u
DST=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/out/shape
mkdir -p "$DST"
n=0
for f in ~/pj/out/phase2_shp*.jsonl; do
  [ -e "$f" ] || continue
  cp "$f" "$DST/$(basename "$f")" || exit 1
  n=$((n+1))
done
echo "pulled $n shard file(s) to $DST"
echo "--- line counts ---"
wc -l "$DST"/phase2_shp*.jsonl 2>/dev/null | tail -8
