#!/bin/bash
# Block until at least $1 cells matching $2 (default B_*) have a RUN_STATUS.json. Echoes what it read.
N=${1:-3}; PAT=${2:-B_}
R=/home/cpb/work/campaigns/jamwin_2026-09-22
while true; do
  c=$(ls -d "$R"/${PAT}*/RUN_STATUS.json 2>/dev/null | wc -l)
  if [ "$c" -ge "$N" ]; then echo "read count=$c (>= $N) for $PAT"; break; fi
  sleep 30
done
