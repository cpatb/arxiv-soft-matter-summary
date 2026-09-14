#!/bin/bash
# THE COMPOSITION CAMPAIGN -- PREREG_MIXTURE_2026-09-13.md.
# Binary packings at FIXED size ratio across the fines volume fraction, one pressure
# rung (p_hat = 1e-3), N = 1500, mu_p = 0.3, 3 seeds. Isotropic compression at this
# rung runs ~60-140 s a cell, so one ratio is a few minutes on 12 workers.
#
# Called with the ratio and tag:  bash scripts/launch_mixture.sh 3.0 mixr3
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
RATIO="${1:?usage: launch_mixture.sh RATIO TAG [EXTRA_ARGS]}"
TAG="${2:?usage: launch_mixture.sh RATIO TAG [EXTRA_ARGS]}"
shift 2
NSH=12
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase12_mixture.py --tag "$TAG" --ratio "$RATIO" \
      --shard "$i" --nshards "$NSH" "$@" \
      > ~/pj/out/${TAG}_w${i}.log 2>&1 &
  pids+=($!)
done
echo "mixture ratio=$RATIO tag=$TAG: ${#pids[@]} workers"
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "ALL DONE"
