#!/bin/bash
# lfsig2 queue (PREREG_LFSIG2_2026-09-27.md), at most P jobs at once.
# Run in the FOREGROUND of a backgrounded wsl.exe call (processes die with their caller).
# Skips cells whose run dir already has DONE (run_cell2.sh does that itself).
set -u
P=${P:-10}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/lfdem/runs/lfsig2
cp $SRC/lfdem/make_conf.py ~/pj/lfdem_make_conf.py
sed 's/\r$//' $SRC/lfdem/run_cell2.sh > ~/lfdem/run_cell2.sh
Q=~/lfdem/queue_lfsig2.txt
sed 's/\r$//' $SRC/lfdem/queue_lfsig2.txt > $Q
md5sum ~/pj/lfdem_make_conf.py $SRC/lfdem/make_conf.py ~/lfdem/LF_DEM/LF_DEM/LF_DEM
echo "lfsig2: $(wc -l < $Q) jobs, P=$P, start $(date -Is), load $(cut -d' ' -f1-3 /proc/loadavg)"
xargs -P $P -I CMD bash -c 'CMD' < $Q
echo "lfsig2 queue complete $(date -Is)"
