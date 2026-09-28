#!/bin/bash
# lfsig2 monodisperse arm (PREREG_LFSIG2 Amendment 4), at most P jobs at once, beside the running lfsig2 driver.
# Run in the FOREGROUND of a detached wsl.exe call. Does NOT touch ~/lfdem/run_cell2.sh or ~/pj/lfdem_make_conf.py
# (both in use by the running driver; make_conf.py is md5-identical to the repo copy, checked below).
set -u
P=${P:-5}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/lfdem/runs/lfsig2m
md5sum ~/pj/lfdem_make_conf.py $SRC/lfdem/make_conf.py
sed 's/\r$//' $SRC/lfdem/run_cell2m.sh > ~/lfdem/run_cell2m.sh
Q=~/lfdem/queue_lfsig2m.txt
sed 's/\r$//' $SRC/lfdem/queue_lfsig2m.txt > $Q
echo "lfsig2m: $(wc -l < $Q) jobs, P=$P, start $(date -Is), load $(cut -d' ' -f1-3 /proc/loadavg)"
xargs -P $P -I CMD bash -c 'CMD' < $Q
echo "lfsig2m queue complete $(date -Is)"
