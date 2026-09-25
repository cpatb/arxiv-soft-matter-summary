#!/bin/bash
# Relaunch the unfinished lfsig cells with a worker cap read from ~/lfdem/MAXJOBS at every launch decision,
# so the cap can be raised or lowered mid-campaign (core sharing with the basilisk session). Same cell selection
# as relaunch_lfsig.sh: skip cells LF_DEM finished ("Job done ok") or whose data reaches strain >= 4; move any
# other existing run directory aside to <label>.killed_<time> and rerun it. Foreground; hold the wsl.exe call.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
R=~/lfdem/runs/lfsig
cp $SRC/lfdem/make_conf.py ~/pj/lfdem_make_conf.py
sed 's/\r$//' $SRC/lfdem/run_cell.sh > ~/lfdem/run_cell.sh.new && mv ~/lfdem/run_cell.sh.new ~/lfdem/run_cell.sh
[ -f ~/lfdem/MAXJOBS ] || echo 6 > ~/lfdem/MAXJOBS
mapfile -t LINES < <(sed 's/\r$//' $SRC/lfdem/queue_lfsig.txt | grep run_cell.sh)
TODO=()
for line in "${LINES[@]}"; do
  lab=$(echo "$line" | awk '{print $3}')
  d=$R/$lab
  if [ -d $d ]; then
    if tail -c 400 $d/run.log 2>/dev/null | grep -q "Job done ok"; then continue; fi
    s=$(grep -v '^#' $d/data_*.dat 2>/dev/null | tail -1 | awk '{print $2}')
    if [ -n "$s" ] && awk "BEGIN{exit !($s >= 3.95)}"; then continue; fi
    mv $d $d.killed_$(date +%H%M%S)
  fi
  TODO+=("$line")
done
echo "relaunch_dyn: ${#TODO[@]} cells, cap from ~/lfdem/MAXJOBS (now $(cat ~/lfdem/MAXJOBS)), $(date -Is)"
for line in "${TODO[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$(cat ~/lfdem/MAXJOBS)" ]; do sleep 20; done
  echo "start $(echo "$line" | awk '{print $3}') $(date -Is)"
  bash -c "$line" &
done
wait
echo "relaunch_dyn complete $(date -Is)"
