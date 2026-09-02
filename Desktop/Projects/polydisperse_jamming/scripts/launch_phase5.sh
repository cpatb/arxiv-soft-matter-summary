#!/bin/bash
# The (J, CV) impact phase diagram. 12 single-threaded workers; see env.sh for why
# many single-threaded workers beat few multi-threaded ones here.
set -u
NSHARDS="${1:-12}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  stdbuf -oL python3 -u phase5_phase_diagram.py --shard "$i" --nshards "$NSHARDS" \
      > ~/pj/out/phase5_w$(printf '%02d' "$i").log 2>&1 &
  pids+=($!)
done
echo "launched ${#pids[@]} phase-5 workers"
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "=== phase 5 complete ==="
