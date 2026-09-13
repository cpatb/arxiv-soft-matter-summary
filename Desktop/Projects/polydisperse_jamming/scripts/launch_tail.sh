#!/bin/bash
# THE TAIL CAMPAIGN -- PREREG_TAIL_2026-09-13.md.
# 72 isotropic-compression cells at N = 1500: four truncated lognormals matched in
# REALISED NUMBER CV = 0.500 but cut at d/d10 <= 2, 3, 5, 8, plus a CV = 0 control at
# the two extreme windows. Isotropic compression is cheap (seconds to a minute a cell),
# so this is ~5 min on 12 workers.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
NSH=12
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase11_tail.py --tag tail --shard "$i" --nshards "$NSH" \
      > ~/pj/out/tail_w${i}.log 2>&1 &
  pids+=($!)
done
echo "tail: ${#pids[@]} workers"
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "ALL DONE"
