#!/bin/bash
# After the phase-6 deep run drains: repeat the scale-invariance test at CV = 0.5 with
# SIX seeds instead of two, to separate a real stiffness dependence from the seed noise
# the CV = 0.5 impact is already known to carry.
set -u
cd ~/pj || exit 1
source scripts/env.sh
n=0
while pgrep -f "phase6_cohesion.py --shard" > /dev/null; do
  n=$((n+1)); echo "[$n] phase6 deep still running: $(pgrep -c -f phase6_cohesion)"; sleep 60
done
echo "=== phase6 deep drained after ~$n min ==="
rm -f out/scale_invariance.jsonl
OMP_NUM_THREADS=2 python3 -u test_scale_invariance.py --cvs 0.5 --seeds 1,2,3,4,5,6 \
    > ~/pj/out/scale_invariance_deep.log 2>&1
tail -24 ~/pj/out/scale_invariance_deep.log
echo "=== DEEP INVARIANCE COMPLETE ==="
