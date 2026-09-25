#!/bin/bash
# Relaunch the UNSTARTED part of the lfsig queue (PREREG_LFSIG Amendment 4). A cell whose run directory
# already exists is skipped: it is either DONE or still running as an orphan of the killed driver, and
# re-running it would duplicate work or clobber a live run. Foreground: run inside a held wsl.exe call.
set -u
P=${P:-9}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp $SRC/lfdem/make_conf.py ~/pj/lfdem_make_conf.py
sed 's/\r$//' $SRC/lfdem/run_cell.sh > ~/lfdem/run_cell.sh.new && mv ~/lfdem/run_cell.sh.new ~/lfdem/run_cell.sh
Q=~/lfdem/queue_lfsig_relaunch.txt
: > $Q
sed 's/\r$//' $SRC/lfdem/queue_lfsig.txt | grep run_cell.sh | while read -r line; do
  lab=$(echo "$line" | awk '{print $3}')
  [ -d ~/lfdem/runs/lfsig/$lab ] && continue
  echo "$line" >> $Q
done
echo "relaunch: $(wc -l < $Q) unstarted cells, P=$P, $(date -Is)"
xargs -P $P -I CMD bash -c 'CMD' < $Q
echo "relaunch complete $(date -Is)"
