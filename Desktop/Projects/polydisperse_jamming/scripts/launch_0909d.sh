#!/bin/bash
# The E4 REPAIR LADDER. E2's pressure root follows the bottom of its own fitting window
# (grade_root_pinning.py); E4, the isostatic crossing of z_nr, does not, because it
# INTERPOLATES inside the ladder and never touches the stress. The P4 frictionless
# anchor is E4's credential: there the crossing of Z_ISO = 6 lands on rung 3 of 5, well
# clear of both ends, and returns 0.6347 against an independently measured 0.6332.
#
# E4 cannot currently measure the production arm because the ladders start AT the
# published isotropic phi_J, where z_nr is already 3.9-4.6 -- so the crossing of the
# FRICTIONAL Z_ISO = 4 sits at or below the bottom rung and E4 returns NaN on 5 of 9
# production groups and on every cv = 0.5 one.
#
# This runs the same cells with the ladder shifted DOWN (offsets -0.06..+0.02), which
# brackets z_nr = 4 at every CV. cv = 0 already exists as `downlad`; this adds cv = 0.2
# and cv = 0.5, 3 seeds, 5 rungs = 30 cells. Tagged `downlad2`; the analysis glob
# `phase8_downlad*.jsonl` matches both it and the cv = 0 `downlad` run, so the three
# CV levels pool, while group_key keeps every knob separate.
set -u
NSHARDS=15
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/dl2_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase8_criticalstate.py --shard "$i" --nshards "$NSHARDS" \
      --cvs 0.2,0.5 --seeds 1,2,3 --offsets -0.06,-0.04,-0.02,0.0,0.02 --tag downlad2 \
      > "$LOG" 2>&1 &
  pids+=($!)
  echo "dl2 worker $i -> $LOG (pid ${pids[-1]})"
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== downlad2 DONE, $fail nonzero exits ==="
