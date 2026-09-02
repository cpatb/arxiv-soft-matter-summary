#!/bin/bash
# Wait for the phase-5 workers to drain, then run phase 7 on the freed cores.
# In a script FILE: an inline `until [ $(...) ]` sent through `wsl.exe -e bash -lc`
# loses its command substitution to quote mangling and spins forever.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

n=0
while pgrep -f "phase5_phase_diagram.py --shard" > /dev/null; do
  n=$((n+1))
  echo "[$n] phase5 still running: $(pgrep -c -f phase5_phase_diagram) workers"
  sleep 60
done
echo "phase5 drained after $n minutes; launching phase 7"

pids=()
for i in $(seq 0 5); do
  stdbuf -oL python3 -u phase7_confined.py --shard "$i" --nshards 6 \
      > ~/pj/out/phase7_w$(printf '%02d' "$i").log 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} phase7 workers"
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "=== phase7 complete ==="
