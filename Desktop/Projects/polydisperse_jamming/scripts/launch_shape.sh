#!/bin/bash
# THE SHAPE CAMPAIGN -- PREREG_SHAPE_2026-09-18.md.
#
#   bash scripts/launch_shape.sh TAG FAMILY ALPHAS PHATS SEEDS CVS N
#
# One family per invocation. Nothing spawned by a wsl.exe call survives that call
# exiting (lessons.md 2026-09-12), so every invocation must finish inside its own call;
# `phase2_shape.py` is idempotent by cell key, so re-invoking with the same tag RESUMES
# rather than repeating, and the campaign can be chopped into calls of any length.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
TAG="${1:?usage: launch_shape.sh TAG FAMILY ALPHAS PHATS SEEDS CVS N}"
FAM="${2:?}"
ALPHAS="${3:?}"
PHATS="${4:?}"
SEEDS="${5:?}"
CVS="${6:?}"
NG="${7:?}"
NSH="${NSH:-20}"
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase2_shape.py --tag "$TAG" --family "$FAM" \
      --alphas "$ALPHAS" --phats "$PHATS" --seeds "$SEEDS" --cvs "$CVS" \
      --mus 0.3 --n "$NG" --shard "$i" --nshards "$NSH" \
      > ~/pj/out/${TAG}_${FAM}_w${i}.log 2>&1 &
  pids+=($!)
done
echo "shape tag=$TAG family=$FAM alphas=$ALPHAS phats=$PHATS seeds=$SEEDS cvs=$CVS N=$NG nshards=$NSH"
rc=0
for p in "${pids[@]}"; do wait "$p" || { echo "worker $p exited nonzero"; rc=1; }; done
echo "LAUNCH_SHAPE DONE tag=$TAG family=$FAM rc=$rc"
exit $rc
