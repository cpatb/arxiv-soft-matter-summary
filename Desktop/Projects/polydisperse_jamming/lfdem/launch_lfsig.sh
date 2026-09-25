#!/bin/bash
# lfsig queue (PREREG_LFSIG_2026-09-24.md) + the deferred cv05margin seeds 8,9, at most P jobs at once.
# Run in the FOREGROUND of a backgrounded wsl.exe call (processes die with their caller).
set -u
P=${P:-10}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp $SRC/lfdem/make_conf.py ~/pj/lfdem_make_conf.py
sed 's/\r$//' $SRC/lfdem/run_cell.sh > ~/lfdem/run_cell.sh
Q=~/lfdem/queue_lfsig.txt
sed 's/\r$//' $SRC/lfdem/queue_lfsig.txt > $Q
echo "lfsig: $(wc -l < $Q) jobs, P=$P, start $(date -Is)"
xargs -P $P -I CMD bash -c 'CMD' < $Q
echo "lfsig queue complete $(date -Is)"
