#!/bin/bash
# TAIL CAMPAIGN, CAP REPAIR -- PREREG_TAIL_2026-09-13.md T1.
#
# The first run refused all four CV 0.5 arms, and the diagnosis is clean: EVERY arm
# converges 3/3 on the p_hat = 1e-3, 3e-4 and 1e-4 rungs and 0/3 on 3e-5, where every
# cell sits exactly at the 200 000 iteration cap with unbalanced 0.014-0.043 against a
# 0.02 tolerance. The cells are slow, not wrong. This is the `capfix` precedent
# (2026-09-09) applied to the same defect.
#
# PJ_MAX_ITER = 1e6 and the SAME tag, so phase11_tail's own skip-by-key resumes: only
# the missing keys run, and the already-good rungs are not recomputed. The cap is
# recorded in every record (`max_iter`), so the grader can see which rung ran under
# which cap.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
NSH=12
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase11_tail.py --tag tailcap --shard "$i" --nshards "$NSH" \
      > ~/pj/out/tailcap_w${i}.log 2>&1 &
  pids+=($!)
done
echo "tailcap: ${#pids[@]} workers at PJ_MAX_ITER=$PJ_MAX_ITER"
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "ALL DONE"
