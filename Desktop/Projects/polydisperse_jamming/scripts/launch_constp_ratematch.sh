#!/bin/bash
# THE RATE-MATCHED CONSTANT-PRESSURE LADDER -- PREREG_CONSTP_RATEMATCH_2026-09-12.md.
#
# gdot is an OUTPUT of the packing (gdot = gamma_max/(steps*dt), dt set by the SMALLEST
# grain), so at a common step count it ran 4.70x across the CV ladder and every CV 0.5
# cell failed the quasi-static bar. Here `shear_steps` is set PER CV from the measured
# dt so all three levels shear at gdot = 1.1747e-04:
#
#     CV 0.0   dt 0.250000   102000 steps   I ~ 6.8e-03   ~410 s
#     CV 0.2   dt 0.127959   200000 steps   I ~ 6.7e-03   ~765 s   <-- the CONTROL
#     CV 0.5   dt 0.053180   480000 steps   I ~ 6.1e-03   ~1410 s
#
# CV 0.2 keeps the production run's exact step count, so its six cells are a BIT-IDENTITY
# control on a protocol whose determinism has never been checked -- and the isotropic
# protocol IS non-deterministic at CV 0 (15 duplicated keys, worst dphi 5.0e-03).
#
# One worker per cell, 18 total. Tag `rm`. Restartable: phase9_constp skips keys already
# in its shard file.
#
# ⚠ Foreground, and the caller must keep polling: nothing spawned by a wsl.exe call
# survives it exiting. The longest cell here is ~1410 s, i.e. longer than one 600 s call,
# so launch this in the background of one call and hold the distro alive with a polling
# loop from the next.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()
launch_cv () {
  cv="$1"; steps="$2"; nick="$3"
  for i in 0 1 2 3 4 5; do
    LOG=~/pj/out/constp_rm_${nick}_w${i}.log
    stdbuf -oL python3 -u phase9_constp.py \
        --cvs "$cv" --mus 0.3 --p-stars 3e-4 --phi-starts 0.55,0.64 \
        --gains 0.15 --seeds 1,2,3 --n 1000 --tag "rm_${nick}" \
        --shear-steps "$steps" --shard "$i" --nshards 6 > "$LOG" 2>&1 &
    pids+=($!)
  done
  echo "cv=$cv steps=$steps -> 6 workers"
}

launch_cv 0.0 102000 c00
launch_cv 0.2 200000 c20
launch_cv 0.5 480000 c50

echo "ratematch: ${#pids[@]} workers -- pids: ${pids[*]}"
for pid in "${pids[@]}"; do
  wait "$pid" || echo "worker $pid exited nonzero"
done
echo "ALL DONE"
