#!/bin/bash
# Reverse-order helper for the lfsig2 queue (execution only). Uses the claim-protected run_cell2.sh, so it and the main
# driver can never run the same cell. Foreground of a detached wsl.exe call.
set -u
P=${P:-5}
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
Q=~/lfdem/queue_lfsig2_rev.txt
sed 's/\r$//' $SRC/lfdem/queue_lfsig2_rev.txt > $Q
echo "helper: $(wc -l < $Q) queue lines (reversed), P=$P, start $(date -Is)"
xargs -P $P -I CMD bash -c 'CMD' < $Q
echo "helper complete $(date -Is)"
