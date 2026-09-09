#!/bin/bash
# The cv = 0.5 half of the ladder extension. NEXT_STEPS_2026-09-09 task 0.2 asked for
# +0.10 and +0.12 at cv = 0 and 0.2 only, to rescue the GATED fits at the low-CV end --
# but adding rungs to two of the three CV levels makes the extension itself CV-
# asymmetric, so a contrast computed on the pooled set mixes "more rungs" with "more
# CV". This adds the same two rungs at cv = 0.5 so all three ladders are 7 rungs.
#
# Tag `prod2ext` so it joins the same groups; runs against the copy already in ~/pj,
# for the same no-concurrent-cp reason as launch_0909b.sh.
set -u
NSHARDS=6
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/ext05_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase8_criticalstate.py --shard "$i" --nshards "$NSHARDS" \
      --cvs 0.5 --seeds 1,2,3 --offsets 0.10,0.12 --tag prod2ext \
      > "$LOG" 2>&1 &
  pids+=($!)
  echo "ext05 worker $i -> $LOG (pid ${pids[-1]})"
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== ext05 DONE, $fail nonzero exits ==="
