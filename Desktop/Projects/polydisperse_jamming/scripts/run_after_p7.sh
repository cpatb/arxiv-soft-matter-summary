#!/bin/bash
# After the phase-7 rerun drains: re-run the scale-invariance test on the fixed force,
# then thicken the CV = 0.5 cohesion cell, which is the thinnest evidence under the
# headline. In a script file because inline waits sent through wsl.exe lose their
# command substitution.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

n=0
while pgrep -f "phase7_confined.py --shard" > /dev/null; do
  n=$((n+1)); echo "[$n] phase7 still running: $(pgrep -c -f phase7_confined)"; sleep 60
done
echo "=== phase7 drained after ~$n min ==="

echo "=== scale invariance on the fixed force ==="
OMP_NUM_THREADS=4 python3 -u test_scale_invariance.py > ~/pj/out/scale_invariance.log 2>&1
tail -30 ~/pj/out/scale_invariance.log

echo "=== phase 6: more seeds at the CV = 0.5 cohesion cell ==="
pids=()
for i in $(seq 0 9); do
  stdbuf -oL python3 -u phase6_cohesion.py --shard "$i" --nshards 10 --tag deep \
      --ratios 100,300,1000,3000 --cvs 0.5 --seeds 1,2,3,4,5,6 \
      > ~/pj/out/p6deep_w$(printf '%02d' "$i").log 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || echo "worker $p nonzero"; done
echo "=== ALL QUEUED WORK COMPLETE ==="
