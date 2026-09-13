#!/bin/bash
# Sync, then run the strain-clock campaign in the foreground. Invoke from PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_strainclock.sh
# (never from Git Bash: it rewrites /mnt/c/... before wsl.exe sees it).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="shear.py packing.py psd.py measure.py phase10_strainclock.py
scripts/launch_strainclock.sh"
mkdir -p ~/pj/out
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
exec bash scripts/launch_strainclock.sh
