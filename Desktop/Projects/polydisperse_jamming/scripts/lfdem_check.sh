#!/bin/bash
# Regression check of the LF_DEM build against the repo's own reference outputs. 2026-09-24.
set -u
cd ~/lfdem/LF_DEM/LF_DEM/check
mkdir -p run && cd run
B=~/lfdem/LF_DEM/LF_DEM/LF_DEM
( time $B -f -r 0.2r -n ../confs/D3N200VF0.45Bidi1.4_0.5Cubic_1.bin ../repulsive.txt > rate.log 2>&1 ) 2>&1 | grep real
( time $B -f -s 5r ../confs/D3N200VF0.45Bidi1.4_0.5Cubic_1.dat ../repulsive_stress.txt > stress.log 2>&1 ) 2>&1 | grep real
ls
for f in data_*.dat; do echo "== diff $f"; diff -I '^# LF_DEM.*' $f ../reference/$f | head -8; echo "lines $(wc -l < $f) ref $(wc -l < ../reference/$f)"; done
tail -3 stress.log
