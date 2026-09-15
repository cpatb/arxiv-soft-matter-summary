#!/bin/bash
# THE K-LADDER / y_S-WIDENING CAMPAIGN -- PREREG_KLADDER_2026-09-15.md.
# Geometric diameter ladders of fixed span 4 at K = 2, 3, 4, 5, crossed with a fines
# volume fraction y_S in {0.03, 0.05, 0.08, 0.13, 0.18}, plus a monodisperse protocol
# control, a degenerate-split control and the NA2 exception at six seeds.
# One pressure rung (p_hat = 1e-3), mu_p = 0.3.
#
#   N = 1500: 20 cells, 63 runs, ~40 s a run  -> ~130 s on 20 workers
#   N = 3000: 23 cells, 72 runs, ~155 s a run -> ~560 s on 20 workers (SPLIT IN TWO)
#
# The 600 s agent cap is per CALL and nothing spawned by a wsl.exe call survives that
# call exiting (lessons.md 2026-09-12), so every invocation must FINISH inside the cap.
# The runner is idempotent -- it skips keys already in its own shard file -- so a second
# call with the same tag resumes rather than repeating.
#
# Called with the tag and any extra runner args:
#   bash scripts/launch_kladder.sh kl
#   bash scripts/launch_kladder.sh kl3k --n 3000 --only-arm K
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
TAG="${1:?usage: launch_kladder.sh TAG [EXTRA_ARGS]}"
shift 1
NSH="${NSH:-20}"
pids=()
for i in $(seq 0 $((NSH-1))); do
  stdbuf -oL python3 -u phase12_kladder.py --tag "$TAG" \
      --shard "$i" --nshards "$NSH" "$@" \
      > ~/pj/out/${TAG}_w${i}.log 2>&1 &
  pids+=($!)
done
echo "kladder tag=$TAG nshards=$NSH: ${#pids[@]} workers, args: $*"
for p in "${pids[@]}"; do wait "$p" || echo "worker $p exited nonzero"; done
echo "ALL DONE"
