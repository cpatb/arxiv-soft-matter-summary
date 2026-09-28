#!/bin/bash
# Pull lfsig2 raw outputs (META, conf.json, params, data_*, int_*, conf.dat, run.log, DONE) from WSL into
# the repo's out/lfsig2/. Raw results flow back; nothing derived is copied. 2026-09-27.
set -u
DST=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/out/lfsig2
mkdir -p $DST
cd ~/lfdem/runs/lfsig2 || exit 1
for d in */; do
  d=${d%/}; mkdir -p $DST/$d
  cp -p $d/META $d/conf.json $d/params.txt $DST/$d/ 2>/dev/null
  cp -p $d/data_*.dat $DST/$d/ 2>/dev/null
  cp -p $d/conf.dat $d/int_*.dat $DST/$d/ 2>/dev/null
  cp -p $d/run.log $DST/$d/ 2>/dev/null
  [ -f $d/DONE ] && touch $DST/$d/DONE
done
echo "synced $(ls $DST | wc -l) cell dirs; DONE: $(ls $DST/*/DONE 2>/dev/null | wc -l); $(date -Is)"
