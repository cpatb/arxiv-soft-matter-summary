#!/bin/bash
# Re-run Phase 7 on the FIXED _com_vel (mass-weighted centre-of-mass velocity).
#
# Every COR at CV > 0 in the previous Phase-7 records was a number-average of grain
# velocities rather than a centre-of-mass velocity -- a CV-dependent error that was
# identically zero at CV = 0 and therefore invisible in every monodisperse control
# (REPORT 7i). The packing channels (spread, Z, z_retention) and ke_ratio are unaffected,
# but the COR column is the one the rebound null is quoted from, so it is re-measured
# rather than left provisional.
#
# The superseded records are MOVED, not deleted, so the before/after comparison stays
# available -- the same treatment the thrust-bug records got.
set -u
NSHARDS="${1:-8}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp "$SRC"/*.py ~/pj/ 2>/dev/null
find ~/pj -maxdepth 1 -name '*.py' | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

mkdir -p out/superseded_numberavg_comvel
mv -f out/phase7_shard*.jsonl out/superseded_numberavg_comvel/ 2>/dev/null

pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/p7rerun_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase7_confined.py --shard "$i" --nshards "$NSHARDS" \
      > "$LOG" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== phase7 rerun finished, $fail nonzero exits ==="
python3 analyse_phase7.py 2>&1 | tail -30
