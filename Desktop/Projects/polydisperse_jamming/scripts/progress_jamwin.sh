#!/bin/bash
# Print ONLY the latest t/tau of each running jamwin cell (health, not values).
R=/home/cpb/work/campaigns/${1:-jamwin}_2026-09-22
for d in "$R"/[ABN]_*/; do
  lab=$(basename "$d")
  [ -f "$d/RUN_STATUS.json" ] && continue
  f=$(ls "$d"/spread_*.dat 2>/dev/null | head -1)
  if [ -n "$f" ]; then
    t=$(grep -v '^#' "$f" | tail -1 | awk '{print $2}')
    echo "$lab running t/tau=$t"
  else
    echo "$lab running (no rows yet)"
  fi
done
