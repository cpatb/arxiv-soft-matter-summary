#!/bin/bash
# cv05margin -- PREREG_CV05MARGIN_2026-09-24.md. 2 cells x seeds {1 (control), 4..9} at p_hat 3e-5, cap 1e6.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
pids=()
for arm in "0.3 0.5 m03c05" "0.0 0.0 m00c00"; do
  set -- $arm
  for s in 1 4 5 6 7 8 9; do
    LOG=~/pj/out/cv05margin_${3}_s${s}.log
    stdbuf -oL python3 -u phase2_phiJ.py --family lognormal --mus "$1" --cvs "$2" --seeds "$s" \
        --phats 3e-5 --tag "cv05margin_${3}_s${s}" > "$LOG" 2>&1 &
    pids+=($!)
  done
done
echo "cv05margin: ${#pids[@]} workers"
for p in "${pids[@]}"; do wait "$p"; done
echo "cv05margin complete $(date -Is)"
