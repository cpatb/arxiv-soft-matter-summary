#!/bin/bash
# Copy jamwin's per-cell spread files + RUN_STATUS + banner lines to the Windows side.
set -u
CAMP=${1:-jamwin}
SRC=/home/cpb/work/campaigns/${CAMP}_2026-09-22
DST=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/out/${CAMP}
mkdir -p "$DST"
n=0
for d in "$SRC"/*/; do
  lab=$(basename "$d")
  [ "$lab" = build ] && continue
  [ -f "$d/RUN_STATUS.json" ] || continue
  mkdir -p "$DST/$lab"
  cp "$d/RUN_STATUS.json" "$DST/$lab/"
  for f in "$d"/spread_*.dat; do [ -f "$f" ] && cp "$f" "$DST/$lab/spread.dat"; done
  grep -a '^#' "$d/err.log" > "$DST/$lab/banner.txt" 2>/dev/null
  n=$((n + 1))
done
REF=/home/cpb/work/campaigns/deltastiff_2026-09-03/DS_d1.3169_L9
mkdir -p "$DST/_ref_DS_d1.3169_L9"
for f in "$REF"/spread_*.dat; do cp "$f" "$DST/_ref_DS_d1.3169_L9/spread.dat"; done
echo "pulled $n finished cells (+ the deltastiff reference for G0)"
