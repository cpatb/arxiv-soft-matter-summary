#!/bin/bash
# Show RUN_STATUS, the err.log tail and the last spread rows for one campaign cell.
#   wsl bash .../scripts/diag_cell.sh <campaign_dir_name> <label>
d=/home/cpb/work/campaigns/$1/$2
echo "== $d"
cat "$d/RUN_STATUS.json"; echo
echo "-- err.log tail"
tail -8 "$d/err.log"
echo "-- last spread rows (cols t, ... )"
for f in "$d"/spread_*.dat; do tail -3 "$f"; done
