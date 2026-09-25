#!/bin/bash
# Relaunch what the lfsig queue still needs (PREREG_LFSIG Amendments 4 and the ~23:05 / ~23:31 kills).
# A cell is SKIPPED if LF_DEM finished it ("Job done ok") or its data already reaches strain >= 4 (graded
# window complete). Anything else -- never started, or killed before strain 4 -- is (re)run from scratch:
# its run directory is moved aside to <label>.killed_<time> first, so no partial output is overwritten.
# Needs the author's OK (the harness reaped the original driver for low memory). Foreground only.
set -u
P=${P:-9}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
R=~/lfdem/runs/lfsig
cp $SRC/lfdem/make_conf.py ~/pj/lfdem_make_conf.py
sed 's/\r$//' $SRC/lfdem/run_cell.sh > ~/lfdem/run_cell.sh.new && mv ~/lfdem/run_cell.sh.new ~/lfdem/run_cell.sh
Q=~/lfdem/queue_lfsig_relaunch.txt
: > $Q
sed 's/\r$//' $SRC/lfdem/queue_lfsig.txt | grep run_cell.sh | while read -r line; do
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
