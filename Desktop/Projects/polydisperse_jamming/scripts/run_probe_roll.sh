#!/bin/bash
# Probe this Yade build's rolling-resistance API. Seconds. From PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_probe_roll.sh
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
for f in probe_roll_api.py; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
python3 -u probe_roll_api.py
echo "exit: $?"
