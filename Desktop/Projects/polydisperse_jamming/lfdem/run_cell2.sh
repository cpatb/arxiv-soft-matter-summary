#!/bin/bash
# One lfsig2 cell (PREREG_LFSIG2_2026-09-27.md): polydisperse config, LF_DEM stress-controlled, outputs kept.
# Usage: run_cell2.sh LABEL CV PHI SIGMA MU SEED STRAIN
# Same construction as run_cell.sh (lfsig) except: no STRAIN_OVERRIDE file, its own runs dir, interaction
# dumps on for every cell. kn rule (registered): kn = max(1e4, 1000*sigma) F_R/a. Wall cap TMAX (default 4 h):
# rc 124 = censored after the window at best, never imputed.
set -u
LABEL=$1; CV=$2; PHI=$3; SIG=$4; MU=$5; SEED=$6; STRAIN=$7
N=${N:-500}
cd ~/pj && source scripts/env.sh >/dev/null 2>&1
B=~/lfdem/LF_DEM/LF_DEM/LF_DEM
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/lfdem
W=~/lfdem/runs/lfsig2/$LABEL
mkdir -p $W && cd $W || exit 1
if [ -f DONE ]; then echo "$LABEL already done"; exit 0; fi
KN=$(python3 -c "print(int(max(1e4, 1000*$SIG)))")
sed -e "s/MU/$MU/g" -e "s/KN/$KN/" -e "s/STRAIN/$STRAIN/" $SRC/params_base.txt | sed 's/\r$//' > params.txt
sed -i 's/^out_data_interaction = false;/out_data_interaction = true;/' params.txt
cd ~/pj && python3 ~/pj/lfdem_make_conf.py --n $N --cv $CV --phi $PHI --seed $SEED --out $W/conf.dat > $W/conf.log 2>&1 || { echo "conf failed"; exit 2; }
cd $W
echo "LABEL=$LABEL CV=$CV PHI=$PHI SIG=$SIG MU=$MU SEED=$SEED STRAIN=$STRAIN N=$N KN=$KN start=$(date -Is) md5bin=$(md5sum $B | cut -c1-12)" > META
timeout ${TMAX:-14400} $B -s ${SIG}r conf.dat params.txt > run.log 2>&1
rc=$?
echo "rc=$rc end=$(date -Is)" >> META
[ $rc -eq 0 ] && touch DONE
exit $rc
