#!/bin/bash
# mustarcap -- PREREG_MUSTARCAP_2026-09-22.md
#
# The mu ladder that grade_mustar.py fits is a MIXED-CAP ladder: mu in {0.6, 3.0} ran at the
# default MAX_ITER = 200000 while mu in {0.0, 1.0} ran at 1000000, and the ladder's single
# monotonicity violation (CV 0.5, mu 0.6 -> 1.0, +0.01633 = 7.8 x PHI_SCALE) sits exactly at
# that cap boundary, with the sign under-compression predicts.
#
# Re-runs lognormal mu_p in {0.0, 0.6, 3.0} x CV {0, 0.2, 0.5} x 3 seeds x 4 rungs = 108
# records at PJ_MAX_ITER = 1000000. mu = 0.3 and mu = 1.0 are deliberately NOT re-run (see
# the prereg's section 1). M1, the identity gate, is primary and gates everything else.
#
# THE WAIT IS LOAD-BEARING (see launch_phase2.sh): this script must stay in the FOREGROUND
# of a held-open background call or WSL kills the worker group. Invoke BY FILE PATH from
# PowerShell, never inline from Git Bash.
set -u
NSHARDS="${1:-12}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out
cp -r "$SRC"/*.py "$SRC"/scripts ~/pj/ 2>/dev/null
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

run_block () {
  local label="$1"; shift
  echo "=== $label ==="
  date -u +"start %Y-%m-%dT%H:%M:%SZ"
  local pids=()
  for i in $(seq 0 $((NSHARDS - 1))); do
    local LOG=~/pj/out/${label}_w$(printf "%02d" "$i").log
    stdbuf -oL python3 -u phase2_phiJ.py "$@" \
        --shard "$i" --nshards "$NSHARDS" > "$LOG" 2>&1 &
    pids+=($!)
  done
  echo "launched ${#pids[@]} workers"
  local fail=0
  for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
  date -u +"end %Y-%m-%dT%H:%M:%SZ"
  echo "=== $label finished, $fail nonzero exits ==="
}

export PJ_MAX_ITER=1000000
run_block mustarcap \
    --family lognormal --mus 0.0,0.6,3.0 --cvs 0.0,0.2,0.5 --seeds 1,2,3 \
    --tag mustarcap
cp ~/pj/out/phase2_mustarcap*.jsonl "$SRC/out/" 2>/dev/null
cp ~/pj/out/mustarcap_w*.log "$SRC/out/" 2>/dev/null

echo "=== mustarcap synced ==="
cat "$SRC"/out/phase2_mustarcap*.jsonl 2>/dev/null | wc -l
