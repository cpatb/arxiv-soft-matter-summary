#!/bin/bash
# shape_status.sh -- what the shape campaign is actually doing right now.
#
# From PowerShell, by file path:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/shape_status.sh
#
# Written as a FILE rather than an inline `wsl bash -lc "..."` because quoting an awk
# program through PowerShell and Git Bash mangles `$1` (lessons.md 2026-09-09: a
# double-quoted wsl invocation lets the outer shell expand `$` first, and every field
# reference arrives EMPTY).
set -u

echo "=== worker processes, by age ==="
ps -eo pid,etimes,args | grep '[p]hase2_shape.py' \
  | awk '{printf "%-8s %6ss  %s\n", $1, $2, $0}' \
  | sed 's/--tag/\n        --tag/' | grep -c . > /dev/null
ps -eo etimes,args | grep -c '[p]hase2_shape.py' | xargs echo "  total workers:"
ps -eo etimes,args | grep '[p]hase2_shape.py' | awk '{print $1}' | sort -n \
  | awk 'NR==1{min=$1} {max=$1; n++} END {print "  age min/max:", min "s /", max "s  over", n, "workers"}'
echo "  distinct families in flight:"
ps -eo args | grep '[p]hase2_shape.py' | grep -o '\-\-family [a-z]*' | sort | uniq -c \
  | sed 's/^/    /'
echo "  distinct shard counts in flight (two different values = two campaigns running):"
ps -eo args | grep '[p]hase2_shape.py' | grep -o '\-\-nshards [0-9]*' | sort | uniq -c \
  | sed 's/^/    /'
echo "  distinct pressure ladders in flight (two values = the old 4-rung run survived):"
ps -eo args | grep '[p]hase2_shape.py' | grep -o '\-\-phats [0-9e.,-]*' | sort | uniq -c \
  | sed 's/^/    /'

echo
echo "=== records on disk ==="
cat ~/pj/out/phase2_shp_*.jsonl 2>/dev/null | wc -l | xargs echo "  phase2_shp_ records:"
for f in lognormal uniform bidisperse; do
  n=$(cat ~/pj/out/phase2_shp_"$f"_shard*.jsonl 2>/dev/null | wc -l)
  echo "    $f: $n"
done

echo
echo "=== wall time of the cells that finished ==="
grep -ho '[0-9.]*s$' ~/pj/out/shp_*_w*.log 2>/dev/null | tr -d 's' | sort -n \
  | awk '{v[n++]=$1} END {if(n){printf "  n=%d  min=%.0fs  median=%.0fs  max=%.0fs\n", n, v[0], v[int(n/2)], v[n-1]} else print "  none yet"}'
