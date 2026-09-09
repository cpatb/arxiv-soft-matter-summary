#!/bin/bash
# The down-ladder that decides PREREG_ROOT_PINNING_2026-09-09.md: is E2's root a
# material critical state, or the bottom of whatever window it is handed?
#
# Same cell as prod2's cv = 0 arm in every respect (family, N, mu_p, seeds, gamma_max,
# shear_steps, protocol) except the ladder OFFSETS, which run -0.06 to +0.02 instead of
# 0.00 to +0.08. Tag `downlad`, which matches no other analysis glob.
#
# ⚠ This does NOT go through scripts/launch_phase8.sh, and deliberately: that launcher
# re-copies $SRC/*.py into ~/pj on every call, and launch_0909.sh is running its own
# stages concurrently. Two concurrent `cp` of the same tree is a race that can hand a
# starting worker a half-written module. This runs against the copy already in ~/pj.
set -u
NSHARDS=8
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/downlad_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase8_criticalstate.py --shard "$i" --nshards "$NSHARDS" \
      --cvs 0.0 --seeds 1,2,3 --offsets -0.06,-0.04,-0.02,0.0,0.02 --tag downlad \
      > "$LOG" 2>&1 &
  pids+=($!)
  echo "downlad worker $i -> $LOG (pid ${pids[-1]})"
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== downlad DONE, $fail nonzero exits ==="
