#!/bin/bash
# THE U-GAP CAMPAIGN -- PREREG_UGAP_2026-09-16.md.
# Binary packings at three size ratios (4^(1/4) = 1.41421, 1.25, 1.5) across a NINE-point
# composition grid that straddles the model's branch crossover, p_hat = 1e-3, mu_p = 0.3,
# 3 seeds. 27 cells, 81 runs at N = 1500.
#
#   N = 1500: 81 runs, ~60-140 s a run  -> ~7-10 min on 20 workers (SPLIT PER RATIO)
#   N = 3000: arm A only, 27 runs       -> ~5 min on 20 workers
#
# Nothing spawned by a wsl.exe call survives that call exiting (lessons.md 2026-09-12),
# so every invocation must finish inside its own call. The runner is idempotent -- it
# skips keys already present in its own shard file -- so re-invoking with the same tag
# RESUMES rather than repeating.
#
# The composition grid is passed on the command line rather than hardcoded here, because
# `ugap.py` is the single source of truth for it and two places that construct the same
# list will eventually disagree (lessons.md 2026-09-03).
#
#   bash scripts/launch_ugap.sh ug 1.4142135623730951 "$XG"
#   bash scripts/launch_ugap.sh ug3k 1.4142135623730951 "$XG" --n 3000
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
TAG="${1:?usage: launch_ugap.sh TAG RATIO XGRID [EXTRA_ARGS]}"
RATIO="${2:?usage: launch_ugap.sh TAG RATIO XGRID [EXTRA_ARGS]}"
XG="${3:?usage: launch_ugap.sh TAG RATIO XGRID [EXTRA_ARGS]}"
shift 3
NSH="${NSH:-20}"
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase12_mixture.py --tag "$TAG" --ratio "$RATIO" \
      --xgrid "$XG" --no-anchor \
      --shard "$i" --nshards "$NSH" "$@" \
      > ~/pj/out/${TAG}_r${RATIO}_w${i}.log 2>&1 &
  pids+=($!)
done
echo "ugap tag=$TAG ratio=$RATIO nshards=$NSH xgrid=$XG: ${#pids[@]} workers, args: $*"
rc=0
for p in "${pids[@]}"; do wait "$p" || { echo "worker $p exited nonzero"; rc=1; }; done
echo "ALL DONE rc=$rc"
exit $rc
