#!/bin/bash
# lfsig2 Amendment 7 binary arm. Foreground of a detached wsl.exe call.
set -u
P=${P:-12}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
sed 's/\r$//' $SRC/lfdem/run_cell_bi.sh > ~/lfdem/run_cell_bi.sh
Q=~/lfdem/queue_lfsig2bi.txt
sed 's/\r$//' $SRC/lfdem/queue_lfsig2bi.txt > $Q
echo "lfsig2bi: $(wc -l < $Q) jobs, P=$P, start $(date -Is), load $(cut -d' ' -f1-3 /proc/loadavg)"
xargs -P $P -I CMD bash -c 'CMD' < $Q
echo "lfsig2bi queue complete $(date -Is)"
