#!/bin/bash
# THE STRAIN-CLOCK CAMPAIGN -- PREREG_STRAINCLOCK_2026-09-13.md.
#
# 45 cells: lognormal, CV {0, 0.2, 0.5} x 5 rungs x 3 seeds, N = 1000, frictionless
# preparation (mu_prep = 0) then friction on then constant-volume shear to gamma = 3 at
# a matched gdot = 1.17e-04. Step counts are set PER CV inside phase10_strainclock.py
# (102k / 200k / 480k) from the measured dt, exactly as the rate-matched constant-
# pressure ladder did.
#
#     CV 0.0   102 000 steps   ~420 s/cell
#     CV 0.2   200 000 steps   ~650 s/cell
#     CV 0.5   480 000 steps   ~850 s/cell
#
# 15 shards, round-robin over a grid ordered CV-outermost, so every worker draws exactly
# one cell from each CV and no shard can be all-CV-0.5. Expected wall ~32 min on 22 cores.
#
# ⚠ Foreground, and the caller must keep the call alive: nothing spawned by a wsl.exe
# call survives it exiting. Launch in the background of one call and poll from the next.
# Restartable: phase10_strainclock skips keys already in its shard file.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

NSHARDS=15
pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/strainclock_w${i}.log
  stdbuf -oL python3 -u phase10_strainclock.py \
      --tag sc --shard "$i" --nshards "$NSHARDS" > "$LOG" 2>&1 &
  pids+=($!)
done
echo "strainclock: ${#pids[@]} workers -- pids: ${pids[*]}"
for pid in "${pids[@]}"; do
  wait "$pid" || echo "worker $pid exited nonzero"
done
echo "ALL DONE"
