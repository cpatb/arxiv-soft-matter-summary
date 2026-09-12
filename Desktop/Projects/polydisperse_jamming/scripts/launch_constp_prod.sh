#!/bin/bash
# THE CONSTANT-PRESSURE PRODUCTION LADDER -- PREREG_CONSTP_PROD_2026-09-12.md §2.
#
# 18 cells: CV {0, 0.2, 0.5} x phi_start {0.55, 0.64} x seeds {1,2,3}, all at
# N = 600, mu_p = 0.3, p* = 3e-04, K = 0.15 (inside the usable gain window the
# shakeout measured), gamma_max = 3, 200000 steps.
#
# ⚠ N = 600, NOT the 1500 of PREREG_CONSTP_PROD_2026-09-12.md §2 -- see its §2a
# AMENDMENT, written before this run. A WSL process does not survive the wsl.exe
# invocation that spawned it (not under nohup, not under setsid nohup, and not with a
# second session held open -- tested with a bare `sleep 300`), and one agent call caps
# at 600 s while an N = 1500 cell costs ~900-1000 s. So no N = 1500 cell can complete
# in this environment. P1 is a CONTRAST at fixed N and is still answerable; P2 is
# DEFERRED.
#
# The grid contains NO phi_J: the same two absolute preparations are used at every CV,
# so nothing is centred on anything. That is the whole point -- delta phi_c from a route
# that takes no input from the isotropic ladder.
#
# 18 workers on 22 cores. Restartable: phase9_constp skips keys already in the shard
# file, so re-running the launcher resumes rather than duplicating.
#
# ⚠ Foreground, no nohup -- WSL2 tears the distro down after the last session closes.
set -u
NSHARDS=18
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()
for i in $(seq 0 $((NSHARDS - 1))); do
  LOG=~/pj/out/constp_prod_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase9_constp.py \
      --cvs 0.0,0.2,0.5 --mus 0.3 --p-stars 3e-4 --phi-starts 0.55,0.64 \
      --gains 0.15 --seeds 1,2,3 --n 1000 --tag prod \
      --shard "$i" --nshards "$NSHARDS" > "$LOG" 2>&1 &
  pids+=($!)
done
echo "constp prod: $NSHARDS workers launched -- pids: ${pids[*]}"
for pid in "${pids[@]}"; do
  wait "$pid" || echo "worker $pid exited nonzero"
done
echo "ALL DONE"
