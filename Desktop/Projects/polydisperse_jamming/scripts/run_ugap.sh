#!/bin/bash
# Sync code into WSL, then run the U-GAP campaign. From PowerShell, BY FILE PATH:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_ugap.sh
#   wsl bash /mnt/c/.../scripts/run_ugap.sh --n 3000        (arm A only; see below)
#
# Git Bash rewrites a /mnt/c/... argument into a Windows path before wsl.exe sees it
# (lessons.md 2026-09-12), so this must be invoked from PowerShell by file path, never
# inline from Git Bash.
#
# The file list is FIXED and explicit. A glob copies whatever happens to be in the
# directory, which is how a WSL tree silently drifts from the repo; naming the files
# makes a new dependency fail LOUDLY at import rather than run against a stale copy.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py phase12_mixture.py ugap.py scripts/launch_ugap.sh"
mkdir -p ~/pj/out ~/pj/scripts
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1

# The grid is the claim. Print its md5 on the WSL side so it can be compared with the
# Windows one rather than assumed equal (lessons.md 2026-09-09), and take the composition
# list from ugap.py itself rather than re-typing it here -- one source of truth.
echo "ugap.py md5 in WSL: $(md5sum ugap.py | cut -d' ' -f1)"
XG=$(python3 -c "import ugap; print(ugap.xgrid_arg())") || exit 1
echo "xgrid from ugap.py: $XG"

if [ "${1:-}" = "--n" ] && [ "${2:-}" = "3000" ]; then
  # Arm A only, on the IDENTICAL cell list, as a labelled finite-size control (U5).
  echo "=== N = 3000, arm A only (tag ug3k) ==="
  bash scripts/launch_ugap.sh ug3k 1.4142135623730951 "$XG" --n 3000
  exit $?
fi

rc=0
for R in 1.4142135623730951 1.25 1.5; do
  echo "=== N = 1500, ratio $R (tag ug) ==="
  bash scripts/launch_ugap.sh ug "$R" "$XG" || rc=1
done
echo "RUN_UGAP DONE rc=$rc"
exit $rc
