#!/bin/bash
# One lfsig2 Amendment 7 binary cell. Usage: run_cell_bi.sh LABEL RATIO XSMALL PHI SIGMA MU SEED STRAIN
# Same LF_DEM binary, params_base.txt, kn rule, interaction dumps and atomic claim as run_cell2.sh; config from make_conf_bi.py.
set -u
LABEL=$1; RATIO=$2; XS=$3; PHI=$4; SIG=$5; MU=$6; SEED=$7; STRAIN=$8
N=${N:-500}
cd ~/pj && source scripts/env.sh >/dev/null 2>&1
B=~/lfdem/LF_DEM/LF_DEM/LF_DEM
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/lfdem
mkdir -p ~/lfdem/runs/lfsig2bi
W=~/lfdem/runs/lfsig2bi/$LABEL
if ! mkdir $W 2>/dev/null; then echo "$LABEL exists; skipped"; exit 0; fi
cd $W || exit 1
KN=$(python3 -c "print(int(max(1e4, 1000*$SIG)))")
sed -e "s/MU/$MU/g" -e "s/KN/$KN/" -e "s/STRAIN/$STRAIN/" $SRC/params_base.txt | sed 's/\r$//' > params.txt
sed -i 's/^out_data_interaction = false;/out_data_interaction = true;/' params.txt
cp $SRC/make_conf.py $SRC/make_conf_bi.py $W/ 2>/dev/null   # provenance copy only
cd ~/pj && python3 $SRC/make_conf_bi.py --n $N --ratio $RATIO --xsmall $XS --phi $PHI --seed $SEED --out $W/conf.dat > $W/conf.log 2>&1 || { echo "conf failed"; exit 2; }
cd $W
echo "LABEL=$LABEL CV=bi RATIO=$RATIO XSMALL=$XS PHI=$PHI SIG=$SIG MU=$MU SEED=$SEED STRAIN=$STRAIN N=$N KN=$KN start=$(date -Is) md5bin=$(md5sum $B | cut -c1-12)" > META
timeout ${TMAX:-14400} $B -s ${SIG}r conf.dat params.txt > run.log 2>&1
rc=$?
echo "rc=$rc end=$(date -Is)" >> META
[ $rc -eq 0 ] && touch DONE
exit $rc
