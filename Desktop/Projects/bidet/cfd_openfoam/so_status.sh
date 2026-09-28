#!/usr/bin/env bash
# so_status.sh -- progress of the PREREG_so_triplet lanes (read-only).
cat ~/bidet_so/lane.log 2>/dev/null
for d in ~/bidet_so/case*; do
  [ -d "$d" ] || continue
  t=$(grep -E '^Time = ' "$d/log.simpleFoam" 2>/dev/null | tail -1)
  c=$(grep -m1 'cells:' "$d/log.checkMesh" 2>/dev/null | awk '{print $2}')
  s=$(grep 'div(phi,k)' "$d/system/fvSchemes" | awk '{print $4,$5}')
  uz=$(grep 'Solving for Uz' "$d/log.simpleFoam" 2>/dev/null | tail -1 | sed 's/.*Initial residual = \([^,]*\).*/\1/')
  pr=$(grep 'Solving for p' "$d/log.simpleFoam" 2>/dev/null | tail -1 | sed 's/.*Initial residual = \([^,]*\).*/\1/')
  ex=$(grep ExecutionTime "$d/log.simpleFoam" 2>/dev/null | tail -1)
  echo "$(basename "$d") cells=$c k-scheme=[$s] $t Uz=$uz p=$pr | $ex"
done
uptime
