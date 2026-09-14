#!/bin/bash
# Sync code into WSL, then run one ratio of the mixture campaign. From PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_mixture.sh 3.0 mixr3
#
# Invoke BY FILE PATH from PowerShell, never inline from Git Bash: Git Bash rewrites a
# /mnt/c/... argument into a Windows path before wsl.exe sees it (lessons.md 2026-09-12).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py phase12_mixture.py scripts/launch_mixture.sh"
mkdir -p ~/pj/out ~/pj/scripts
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
exec bash scripts/launch_mixture.sh "$@"
