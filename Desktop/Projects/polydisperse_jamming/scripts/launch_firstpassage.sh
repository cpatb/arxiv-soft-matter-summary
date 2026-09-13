#!/bin/bash
# THE FIRST-PASSAGE EXTENSION -- PREREG_FIRSTPASSAGE_2026-09-13.md §2.
#
# 54 cells, same driver / protocol / rate / N as the strain-clock campaign, so the 45
# committed records pool with these unchanged:
#
#   A  the NEW rung x = +0.05, seeds 1-6            3 CV x 6 = 18 cells   tag fpA
#   B  seeds 4-6 at x = +0.10, 0.35, 0.60, 0.85     3 CV x 12 = 36 cells  tag fpB
#
# giving 6 seeds at every rung from x = +0.05 to +0.85.
#
# 18 shards per arm over a grid ordered CV-outermost, so each worker draws one cell per
# CV and no shard is all-CV-0.5 (which costs ~1400 s/cell against ~290 s at CV 0).
# Arms run SEQUENTIALLY, 18 workers each: 22 cores, and running both at once would
# oversubscribe. Expected wall ~13 min (A) + ~38 min (B).
#
# ⚠ Foreground. Nothing spawned by a wsl.exe call survives it exiting, so launch this in
# the background of one call and poll from the next. Restartable: the driver skips keys
# already in its shard file.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

run_arm () {
  tag="$1"; xs="$2"; seeds="$3"; nsh="$4"
  pids=()
  for i in $(seq 0 $((nsh - 1))); do
    LOG=~/pj/out/firstpassage_${tag}_w${i}.log
    stdbuf -oL python3 -u phase10_strainclock.py \
        --tag "$tag" --shard "$i" --nshards "$nsh" \
        --x-rungs "$xs" --seeds "$seeds" > "$LOG" 2>&1 &
    pids+=($!)
  done
  echo "arm $tag: ${#pids[@]} workers -- xs=$xs seeds=$seeds"
  for pid in "${pids[@]}"; do
    wait "$pid" || echo "arm $tag worker $pid exited nonzero"
  done
  echo "arm $tag DONE"
}

run_arm fpA 0.05 1,2,3,4,5,6 18
run_arm fpB 0.10,0.35,0.60,0.85 4,5,6 18
echo "ALL DONE"
