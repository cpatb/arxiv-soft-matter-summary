#!/bin/bash
# Sync code into WSL, then run the K-ladder campaign. From PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_kladder.sh kl
#   wsl bash /mnt/c/.../scripts/run_kladder.sh kl3k --n 3000
#
# Invoke BY FILE PATH from PowerShell, never inline from Git Bash: Git Bash rewrites a
# /mnt/c/... argument into a Windows path before wsl.exe sees it (lessons.md 2026-09-12).
#
# The file list is FIXED and explicit, exactly as run_ternary.sh's is. A glob would copy
# whatever happens to be in the directory, which is how a WSL tree silently drifts from
# the repo; naming the files means a new dependency fails LOUDLY at import rather than
# running against a stale copy.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py phase12_mixture.py kladder.py phase12_kladder.py scripts/launch_kladder.sh"
mkdir -p ~/pj/out ~/pj/scripts
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
# The grid is the claim; print its md5 on the WSL side so it can be compared with the
# Windows one rather than assumed equal (lessons.md 2026-09-09).
echo "kladder.py md5 in WSL: $(md5sum kladder.py | cut -d' ' -f1)"
exec bash scripts/launch_kladder.sh "$@"
