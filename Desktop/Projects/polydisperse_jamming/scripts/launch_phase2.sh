#!/bin/bash
# Launch the Phase-2 phi_J sweep across NSHARDS parallel Yade workers, and WAIT.
#
# The wait is load-bearing. `wsl.exe -e bash -lc '... &'` returns as soon as the
# launcher shell exits, and WSL then kills the process group -- three of nine workers
# died that way on the first attempt (the same failure mode that killed the userspace
# Yade install). Keeping this script in the foreground of a held-open background Bash
# call keeps the whole group alive.
#
# 16 physical cores (Ultra 7 155H: 6 P + 8 E + 2 LP-E), 22 threads. MEASURED: 8 workers
# x 2 OMP threads gave load 42 and 0.43 cells/min, because idle OpenMP threads spin.
# 12 single-threaded workers with OMP_WAIT_POLICY=PASSIVE is the right shape for a
# sweep of many small independent cells.
#
# Each worker writes its OWN jsonl: appending to one file from 8 processes would
# interleave partial lines. analyse_phase2.py globs them back together.
set -u
NSHARDS="${1:-12}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/phase2_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase2_phiJ.py --shard "$i" --nshards "$NSHARDS" \
      > "$LOG" 2>&1 &
  pids+=($!)
  echo "worker $i -> $LOG (pid ${pids[-1]})"
done
echo "launched ${#pids[@]} workers, OMP_NUM_THREADS=$OMP_NUM_THREADS; waiting..."
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== all workers finished, $fail nonzero exits ==="
