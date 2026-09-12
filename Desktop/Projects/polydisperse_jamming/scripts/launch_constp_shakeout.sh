#!/bin/bash
# THE CONSTANT-PRESSURE SERVO SHAKEOUT -- PREREG_CONSTP_2026-09-12.md §4.
#
# 6 cells at N = 600: phi_start in {0.55, 0.64} x gain in {0.25, 0.5, 1.0}, all at
# cv = 0, mu_p = 0.3, p* = 3e-04 (a rung of the isotropic ladder). Nothing here is a
# measurement of phi_c; the run decides whether the servo is admissible at all:
#
#   S2  phi_ss must NOT move with the gain     (else it tracks the controller)
#   S3  phi_ss must NOT remember phi_start     (else gamma_max = 3 is too little strain)
#
# One shard per cell, so a cell that hangs cannot hold up the other five.
#
# The _prepare refactor is verified separately by scripts/run_prepare_check.sh, which
# runs the PRE-refactor shear.py (from a sidecar generated out of git) and the current
# one on the same cell and requires bit-identical output on the continuous channels.
#
# ⚠ NO nohup. A nohup'd process here does NOT survive: WSL2 tears the distro down
# shortly after the last session closes, and the 06:47 attempt on 2026-09-12 left a
# zero-byte launch log and no workers. The caller must hold the session open for the
# duration -- this script waits on every worker and is meant to be run in the
# foreground.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()

for i in 0 1 2 3 4 5; do
  LOG=~/pj/out/constp_shake_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase9_constp.py --shakeout --shard "$i" --nshards 6 \
      > "$LOG" 2>&1 &
  pids+=($!)
done
echo "shakeout: 6 workers launched"

echo "pids: ${pids[*]}"
for pid in "${pids[@]}"; do
  wait "$pid" || echo "worker $pid exited nonzero"
done
echo "ALL DONE"
