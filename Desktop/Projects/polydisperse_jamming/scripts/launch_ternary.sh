#!/bin/bash
# THE TERNARY ADDITIVITY CAMPAIGN -- PREREG_TERNARY_2026-09-14.md.
# Three-class packings at d = 1, 2, 4 across six registered compositions plus a
# monodisperse protocol control, one pressure rung (p_hat = 1e-3), N = 1500,
# mu_p = 0.3, 3 seeds = 21 cells. Isotropic compression at this rung runs ~60-140 s a
# cell in the binary campaign, so 21 cells on 12 workers is a few minutes.
#
# Called with the tag:  bash scripts/launch_ternary.sh tern
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
TAG="${1:?usage: launch_ternary.sh TAG [EXTRA_ARGS]}"
shift 1
NSH=12
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase12_ternary.py --tag "$TAG" \
      --shard "$i" --nshards "$NSH" "$@" \
      > ~/pj/out/${TAG}_w${i}.log 2>&1 &
  pids+=($!)
done
echo "ternary tag=$TAG: ${#pids[@]} workers"
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "ALL DONE"
