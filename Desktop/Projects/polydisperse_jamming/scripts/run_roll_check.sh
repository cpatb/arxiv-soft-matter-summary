#!/bin/bash
# Three gates on the rolling-resistance path -- PREREG_ROLLING_2026-09-13.md R-A/R-B/R-C.
# Four cells at N = 600, a few minutes. From PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_roll_check.sh
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py check_roll_path.py packing_pre_roll_src.txt"
mkdir -p ~/pj/out
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
stdbuf -oL python3 -u check_roll_path.py --sidecar packing_pre_roll_src.txt 2>&1 \
    | tee ~/pj/out/roll_check.log
echo "exit: ${PIPESTATUS[0]}"
