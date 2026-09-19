#!/bin/bash
# PREREG_WCBIDI_2026-09-19.md -- the frictionless bidisperse block.
#
# 72 records = 6 CV x 3 seeds x 4 p_hat rungs, family=bidisperse, mu_p=0.0 ONLY, N=1500,
# tag `wcbidi` (matches no existing tag, so these can never pool with the mu_p=0.3 block).
#
# THE WAIT IS LOAD-BEARING, same as launch_phase2.sh: `wsl.exe -e bash -lc '... &'` returns
# as soon as the launcher shell exits and WSL kills the process group. This script must stay
# in the FOREGROUND of a held-open background Bash call.
#
# Invoke BY FILE PATH from PowerShell, never inline from Git Bash (Git Bash rewrites a
# /mnt/c argument into a Windows path before wsl.exe sees it).
set -u
NSHARDS="${1:-12}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

echo "=== wcbidi: bidisperse, mu_p=0.0, CV 0.0-0.5, seeds 1-3, 4 rungs, N=1500 ==="
date -u +"start %Y-%m-%dT%H:%M:%SZ"

pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/wcbidi_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase2_phiJ.py \
      --family bidisperse \
      --mus 0.0 \
      --cvs 0.0,0.1,0.2,0.3,0.4,0.5 \
      --seeds 1,2,3 \
      --tag wcbidi \
      --shard "$i" --nshards "$NSHARDS" \
      > "$LOG" 2>&1 &
  pids+=($!)
  echo "worker $i -> $LOG (pid ${pids[-1]})"
done
echo "launched ${#pids[@]} workers, OMP_NUM_THREADS=$OMP_NUM_THREADS; waiting..."
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
echo "=== all workers finished, $fail nonzero exits ==="

# Sync the shards back to the Windows tree. RAW results only -- derived tables are
# re-derived, never copied (2026-09-08 lesson: sync_back.sh clobbered fifteen committed
# outputs with stale versions).
cp ~/pj/out/phase2_wcbidi*.jsonl "$SRC/out/" 2>/dev/null
cp ~/pj/out/wcbidi_w*.log "$SRC/out/" 2>/dev/null
ls -la "$SRC"/out/phase2_wcbidi*.jsonl 2>/dev/null | head -20
echo "=== synced ==="
