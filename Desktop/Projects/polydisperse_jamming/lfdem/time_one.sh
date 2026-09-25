#!/bin/bash
# Timing/smoke run of LF_DEM on a polydisperse config. 2026-09-24.
set -u
source ~/miniforge3/etc/profile.d/conda.sh
cd ~/pj && source scripts/env.sh >/dev/null 2>&1
B=~/lfdem/LF_DEM/LF_DEM/LF_DEM
W=~/lfdem/runs/smoke; mkdir -p $W; cd $W
cp /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/lfdem/make_conf.py ~/pj/lfdem_make_conf.py
sed -e "s/MU/${MU:-0.3}/g" -e "s/KN/${KN:-10000}/" -e "s/STRAIN/${STRAIN:-1}/" /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/lfdem/params_base.txt > $W/params.txt; grep -E 'mu_|kn|time_end' $W/params.txt
( time python3 ~/pj/lfdem_make_conf.py --n ${N:-500} --cv ${CV:-0.2} --phi ${PHI:-0.54} --seed 1 --out $W/conf.dat ) 2>&1 | tail -4
( time timeout 540 $B -s ${SIG:-3r} conf.dat params.txt > run.log 2>&1 ) 2>&1 | grep real
tail -5 run.log
ls
f=$(ls data_*.dat | head -1); grep -v '^#' $f | wc -l; grep '^#' $f | head -40 | tail -35; grep -v '^#' $f | tail -3 | cut -c1-200
