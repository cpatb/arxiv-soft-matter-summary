#!/bin/bash
# Relaunch what the lfsig2 queue still needs after a driver kill. A cell is SKIPPED if LF_DEM finished it
# ("Job done ok") or its data already reaches strain >= 3.95 (graded window complete). Anything else is
# (re)run from scratch; its run dir is moved aside to <label>.killed_<time> first. Deterministic cells
# cannot change their answer by being rerun (PREREG_LFSIG2 "What would change the plan"). Foreground only.
set -u
P=${P:-10}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
R=~/lfdem/runs/lfsig2
cp $SRC/lfdem/make_conf.py ~/pj/lfdem_make_conf.py
sed 's/\r$//' $SRC/lfdem/run_cell2.sh > ~/lfdem/run_cell2.sh.new && mv ~/lfdem/run_cell2.sh.new ~/lfdem/run_cell2.sh
Q=~/lfdem/queue_lfsig2_relaunch.txt
: > $Q
sed 's/\r$//' $SRC/lfdem/queue_lfsig2.txt | grep run_cell2.sh | while read -r line; do
  lab=$(echo "$line" | awk '{print $3}')
  d=$R/$lab
  if [ -d $d ]; then
    if tail -c 400 $d/run.log 2>/dev/null | grep -q "Job done ok"; then continue; fi
    s=$(grep -v '^#' $d/data_*.dat 2>/dev/null | tail -1 | awk '{print $2}')
    if [ -n "$s" ] && awk "BEGIN{exit !($s >= 3.95)}"; then continue; fi
    mv $d $d.killed_$(date +%H%M%S)
  fi
  echo "$line" >> $Q
done
echo "relaunch: $(wc -l < $Q) cells to run, P=$P, $(date -Is)"
cat $Q | awk '{print "  ", $3}'
xargs -P $P -I CMD bash -c 'CMD' < $Q
echo "relaunch complete $(date -Is)"
