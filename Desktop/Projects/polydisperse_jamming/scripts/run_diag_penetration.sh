#!/bin/bash
set -u
cd ~/pjdt || exit 1
tr -d '\r' < /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/diag_penetration.py > diag_penetration.py
source scripts/env.sh
export OMP_NUM_THREADS=1 PJ_MAX_ITER=1000000
python3 diag_penetration.py 2>&1 | tail -12
