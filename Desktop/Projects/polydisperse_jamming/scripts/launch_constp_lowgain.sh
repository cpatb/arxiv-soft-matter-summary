#!/bin/bash
# LOW-GAIN EXTENSION of the constant-pressure shakeout.
#
# PREREG_CONSTP_2026-09-12.md §5 S1 says, in advance: "Fewer than 5 means the gain range
# is wrong and the grid is re-chosen, not the bar." It came back 2 of 6, and the
# recorded series say why -- at K >= 0.5 the cell PERIODICALLY UNJAMS (p_hat down to
# 1e-08, six decades below target, and z_nr down to 0.02), so those cells are not
# droopy, they are unstable, and they are not measurements of anything.
#
# That also re-opens S2: gain-independence passed over a set that INCLUDES the unstable
# cells, and only one gain (0.25) was stable, so it cannot be read on stable cells
# alone. This adds K = 0.05 and 0.1 at both preparations so that S2 has three STABLE
# gains to be read on.
#
# 4 cells, one shard each. Tag `lowgain`; the grader reads phase9_shakeout* and
# phase9_lowgain* as one set and says which cells came from which.
#
# ⚠ Foreground, no nohup -- see launch_constp_shakeout.sh.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()
for i in 0 1 2 3; do
  LOG=~/pj/out/constp_lowgain_w$(printf "%02d" "$i").log
  stdbuf -oL python3 -u phase9_constp.py \
      --cvs 0.0 --mus 0.3 --p-stars 3e-4 --phi-starts 0.55,0.64 \
      --gains 0.05,0.1 --seeds 1 --n 600 --tag lowgain \
      --shard "$i" --nshards 4 > "$LOG" 2>&1 &
  pids+=($!)
done
echo "lowgain: 4 workers launched -- pids: ${pids[*]}"
for pid in "${pids[@]}"; do
  wait "$pid" || echo "worker $pid exited nonzero"
done
echo "ALL DONE"
