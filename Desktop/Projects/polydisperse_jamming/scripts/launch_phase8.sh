#!/bin/bash
# Launch the Phase-8 critical-state (constant-volume shear) sweep across NSHARDS
# workers and WAIT. Same shape and same reasoning as launch_phase2.sh: `wsl.exe -e
# bash -lc '... &'` returns and WSL kills the process group, so the wait is what keeps
# the workers alive. Extra args after NSHARDS are passed through to every worker.
#
#   scripts/launch_phase8.sh 12 --cvs 0.0,0.5 --seeds 1 --tag locate
set -u
NSHARDS="${1:-12}"; shift || true
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/phase8_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase8_criticalstate.py --shard "$i" --nshards "$NSHARDS" "$@" \
      > "$LOG" 2>&1 &
  pids+=($!)
  echo "worker $i -> $LOG (pid ${pids[-1]})"
done
echo "launched ${#pids[@]} workers, OMP_NUM_THREADS=$OMP_NUM_THREADS; waiting..."
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== all workers finished, $fail nonzero exits ==="
