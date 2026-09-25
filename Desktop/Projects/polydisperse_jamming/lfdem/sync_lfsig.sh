#!/bin/bash
# Pull lfsig raw outputs (META, conf.json, data_*, st_*) from WSL into the repo's out/lfsig/.
# Raw results flow back; nothing derived is copied. 2026-09-24.
set -u
DST=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/out/lfsig
mkdir -p $DST
cd ~/lfdem/runs/lfsig || exit 1
for d in */; do
  d=${d%/}; mkdir -p $DST/$d
  cp -p $d/META $d/conf.json $d/params.txt $DST/$d/ 2>/dev/null
  cp -p $d/data_*.dat $DST/$d/ 2>/dev/null
  cp -p $d/conf.dat $d/int_*.dat $DST/$d/ 2>/dev/null
  [ -f $d/DONE ] && touch $DST/$d/DONE
done
ls $DST | wc -l
