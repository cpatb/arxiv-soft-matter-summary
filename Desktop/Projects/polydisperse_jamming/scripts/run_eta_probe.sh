#!/bin/bash
# Scope the rolling-resistance lever. ~1 min. From PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_eta_probe.sh
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
for f in packing.py psd.py measure.py probe_eta_ladder.py; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
stdbuf -oL python3 -u probe_eta_ladder.py 2>&1 | tee ~/pj/out/eta_probe.log
echo "exit: ${PIPESTATUS[0]}"
