#!/bin/bash
# Run two sweeps side by side, 6 single-threaded workers each on a 16-core box.
#   ~/lb.sh <scriptA> <nA> <tagA> <scriptB> <nB> <tagB> [-- extra args for B]
# Waits for both, so the whole group survives (WSL kills orphaned process groups).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()
launch () {
  local script="$1" n="$2" tag="$3"; shift 3
  for i in $(seq 0 $((n - 1))); do
    stdbuf -oL python3 -u "$script" --shard "$i" --nshards "$n" "$@" \
        > ~/pj/out/"${tag}_w$(printf '%02d' "$i")".log 2>&1 &
    pids+=($!)
  done
  echo "launched $n workers for $tag ($script $*)"
}

launch phase5_phase_diagram.py 6 phase5
launch phase6_cohesion.py 6 phase6

echo "waiting on ${#pids[@]} workers..."
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "=== batch complete ==="
