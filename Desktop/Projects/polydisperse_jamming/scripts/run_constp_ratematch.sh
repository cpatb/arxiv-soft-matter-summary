#!/bin/bash
# Sync, then run the rate-matched ladder in the foreground. Invoke from PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_constp_ratematch.sh
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="shear.py phase9_constp.py grade_constp.py grade_constp_prod.py
scripts/launch_constp_ratematch.sh"
mkdir -p ~/pj/out
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
exec bash scripts/launch_constp_ratematch.sh
