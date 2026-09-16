#!/bin/bash
# THE FABRIC CAMPAIGN -- PREREG_FABRIC_2026-09-18.md.
#
# One family per invocation, at ONE p_hat rung, matched mu_p = 0.3, N = 1500. A single
# rung is NOT a phi_J measurement and the prereg says so; this campaign reads fabric
# against CV at matched pressure.
#
#   bash scripts/launch_fabric.sh TAG FAMILY PHAT SEEDS CVS
#   bash scripts/launch_fabric.sh fab lognormal 1e-3 1,2,3,4,5,6,7,8,9,10,11,12 0.0,0.1,0.2,0.3,0.4,0.5
#
# Nothing spawned by a wsl.exe call survives that call exiting (lessons.md 2026-09-12),
# so every invocation must finish inside its own call. The runner is idempotent -- it
# skips keys already present in its own shard file -- so re-invoking with the same tag
# RESUMES rather than repeating.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
TAG="${1:?usage: launch_fabric.sh TAG FAMILY PHAT SEEDS CVS}"
FAM="${2:?usage: launch_fabric.sh TAG FAMILY PHAT SEEDS CVS}"
PHAT="${3:?usage: launch_fabric.sh TAG FAMILY PHAT SEEDS CVS}"
SEEDS="${4:?usage: launch_fabric.sh TAG FAMILY PHAT SEEDS CVS}"
CVS="${5:?usage: launch_fabric.sh TAG FAMILY PHAT SEEDS CVS}"
NSH="${NSH:-20}"
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase2_phiJ.py --tag "$TAG" --family "$FAM" \
      --phats "$PHAT" --seeds "$SEEDS" --cvs "$CVS" --mus 0.3 --n 1500 \
      --shard "$i" --nshards "$NSH" \
      > ~/pj/out/${TAG}_${FAM}_p${PHAT}_w${i}.log 2>&1 &
  pids+=($!)
done
echo "fabric tag=$TAG family=$FAM p_hat=$PHAT seeds=$SEEDS cvs=$CVS nshards=$NSH"
rc=0
for p in "${pids[@]}"; do wait "$p" || { echo "worker $p exited nonzero"; rc=1; }; done
echo "LAUNCH_FABRIC DONE tag=$TAG family=$FAM p_hat=$PHAT rc=$rc"
exit $rc
