#!/bin/bash
# Two campaigns, SEQUENTIALLY so they never compete for the 16 physical cores.
#
#   1) wcbidicap -- PREREG_WCBIDI_2026-09-19.md Amendment 3. The full bidisperse mu_p = 0
#      grid re-run at PJ_MAX_ITER = 1000000, to rescue the 14 cells that capped at 200000.
#      Primary gate C1 is an IDENTITY gate: the 58 cells that already converged must come
#      back BIT-IDENTICAL (a raised limit cannot move a cell that never reached it).
#
#   2) mustar   -- PREREG_MUSTAR_2026-09-19.md. mu_p in {0.6, 3.0} x CV {0, 0.2, 0.5},
#      lognormal, to separate (M1) a CV-dependent friction saturation scale mu* from
#      (M2) a CV-dependent amplitude.
#
# THE WAIT IS LOAD-BEARING (see launch_phase2.sh): this script must stay in the FOREGROUND
# of a held-open background Bash call or WSL kills the worker group. Invoke BY FILE PATH
# from PowerShell, never inline from Git Bash.
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

# --- 1) the cap re-run. PJ_MAX_ITER is read by phase2_phiJ from the environment. ---------
export PJ_MAX_ITER=1000000
run_block wcbidicap \
    --family bidisperse --mus 0.0 --cvs 0.0,0.1,0.2,0.3,0.4,0.5 --seeds 1,2,3 \
    --tag wcbidicap
cp ~/pj/out/phase2_wcbidicap*.jsonl "$SRC/out/" 2>/dev/null

# --- 2) mustar. Back to the DEFAULT cap: this campaign is not about caps. -----------------
unset PJ_MAX_ITER
run_block mustar \
    --family lognormal --mus 0.6,3.0 --cvs 0.0,0.2,0.5 --seeds 1,2,3 \
    --tag mustar
cp ~/pj/out/phase2_mustar*.jsonl "$SRC/out/" 2>/dev/null
cp ~/pj/out/wcbidicap_w*.log ~/pj/out/mustar_w*.log "$SRC/out/" 2>/dev/null

echo "=== both campaigns synced ==="
ls -la "$SRC"/out/phase2_wcbidicap*.jsonl "$SRC"/out/phase2_mustar*.jsonl 2>/dev/null | wc -l
