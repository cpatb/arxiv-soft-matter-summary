#!/bin/bash
# Chain the remaining compute so it proceeds without supervision:
#
#   1. wait for the running Phase-2 lognormal sweep to drain
#   2. Phase 2b: the UNIFORM (zero-skewness) family, to separate the CV axis from the
#      skewness axis that Desmond & Weeks (PRE 90, 022204) show can dominate it
#   3. Phase 4: the drop-like boundary
#   4. full analysis chain + figures + sync back to Windows
#
# Each stage is resumable (cells already in the jsonl are skipped), so killing this
# script and re-running it loses at most the cells in flight.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

echo "=== $(date -Is) stage 0: waiting for the phase-2 lognormal sweep ==="
while pgrep -f "phase2_phiJ.py --shard" > /dev/null; do sleep 30; done
echo "=== $(date -Is) lognormal sweep drained ==="

run_sharded () {
  local script="$1"; shift
  local n="$1"; shift
  local tag="$1"; shift
  local pids=()
  for i in $(seq 0 $((n - 1))); do
    stdbuf -oL python3 -u "$script" --shard "$i" --nshards "$n" "$@" \
        > ~/pj/out/"${tag}_w$(printf '%02d' "$i")".log 2>&1 &
    pids+=($!)
  done
  echo "  launched ${#pids[@]} workers for $tag"
  for p in "${pids[@]}"; do wait "$p" || echo "  worker $p exited nonzero"; done
  echo "=== $(date -Is) $tag finished ==="
}

echo "=== $(date -Is) stage 1: phase 2b, UNIFORM family (zero skewness) ==="
# mu_p = 0.3 only: the uniform family exists to answer ONE question (is the
# effect CV or skewness?) and it only needs the frictional column the lab
# actually cares about. 6 CV x 1 mu x 3 seeds x 4 p_hat = 72 cells.
run_sharded phase2_phiJ.py 12 phase2b --family uniform --mus 0.3

echo "=== $(date -Is) stage 2: phase 4, the drop-like boundary ==="
run_sharded phase4_impact.py 10 phase4 --n 1500

echo "=== $(date -Is) stage 3: analysis + figures ==="
bash ~/aa.sh

echo "=== $(date -Is) ALL STAGES COMPLETE ==="
