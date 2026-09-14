#!/bin/bash
# Sync code into WSL, then run the ternary additivity campaign. From PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_ternary.sh tern
#
# Invoke BY FILE PATH from PowerShell, never inline from Git Bash: Git Bash rewrites a
# /mnt/c/... argument into a Windows path before wsl.exe sees it (lessons.md 2026-09-12).
#
# The file list is FIXED and explicit, exactly as run_mixture.sh's is. A glob would copy
# whatever happens to be in the directory, which is how a WSL tree silently drifts from
# the repo; naming the files means a new dependency fails LOUDLY at import rather than
# running against a stale copy.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py phase12_mixture.py phase12_ternary.py scripts/launch_ternary.sh"
mkdir -p ~/pj/out ~/pj/scripts
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
exec bash scripts/launch_ternary.sh "$@"
