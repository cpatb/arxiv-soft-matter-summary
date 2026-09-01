#!/bin/bash
# Sync the Windows tree, then run a python script with output tee'd to a pollable log.
#   ~/rl.sh <script.py> [args...]
# Log lands at ~/pj/out/<script>.log and is flushed line-by-line (python -u).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
SCRIPT="$1"; shift
LOG=~/pj/out/"$(basename "$SCRIPT" .py)".log
echo "=== $(date -Is) start $SCRIPT $* ===" > "$LOG"
stdbuf -oL -eL python3 -u "$SCRIPT" "$@" 2>&1 | stdbuf -oL tee -a "$LOG"
echo "=== $(date -Is) exit=$? ===" >> "$LOG"
