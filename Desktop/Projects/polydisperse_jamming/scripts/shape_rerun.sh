#!/bin/bash
# shape_rerun.sh -- the CV = 0.5 arm under the RAISED cap, PREREG_SHAPE amendment 5.10.
#
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/shape_rerun.sh
#
# Invoked by file path from PowerShell (never inline: Git Bash rewrites /mnt/c and an
# inline `-lc` string loses its command substitutions). Everything this script spawns
# dies when the calling session ends, and a cell writes its record ONLY on completion,
# so the longest cells are launched FIRST and the script is re-runnable: the runner
# resumes by key.
#
# PJ_MAX_ITER is exported HERE, not in the PowerShell environment -- `$env:X` does not
# cross into WSL (it did not on 2026-09-16: NSH=14 gave 19 workers).
set -u
cd ~/pj || exit 1
export PJ_MAX_ITER=1000000     # amendment 5.10: 12/12 predicted cells with 30 % margin
TAG=shpcap                     # production (not in GATE_TAGS); files match phase2_shp*
NG=5400
SEEDS=1,2,3
LOG=~/pj/out/shape_rerun.log

echo "shape_rerun: PJ_MAX_ITER=$PJ_MAX_ITER tag=$TAG start $(date -Is)" | tee -a "$LOG"

# The registered identity check (5.9/5.10): two CV = 0 cells at the new cap, which must
# come back bit-identical to the default-cap records. Cheapest cells on the grid.
NSH=2 bash scripts/launch_shape.sh "$TAG" lognormal 1.0,2.0 1e-3 1 0.0 "$NG" \
    >> "$LOG" 2>&1 &
pid_ident=$!

# Longest first. lognormal alpha 2, all rungs (seed 3 never ran; seeds 1-2 sat at cap).
NSH=9 bash scripts/launch_shape.sh "$TAG" lognormal 2.0 1e-3,3e-4,1e-4 "$SEEDS" 0.5 "$NG" \
    >> "$LOG" 2>&1 &
pid_ln2=$!
# uniform alpha 2, all rungs (worst predicted cell, 772k bias-corrected).
NSH=9 bash scripts/launch_shape.sh "$TAG" uniform 2.0 1e-3,3e-4,1e-4 "$SEEDS" 0.5 "$NG" \
    >> "$LOG" 2>&1 &
pid_un2=$!
# uniform alpha 1 at the two rungs that hit the cap.
NSH=6 bash scripts/launch_shape.sh "$TAG" uniform 1.0 3e-4,1e-4 "$SEEDS" 0.5 "$NG" \
    >> "$LOG" 2>&1 &
pid_un1=$!
# lognormal alpha 1 at 1e-4 (needed 203k against a cap of 200k).
NSH=3 bash scripts/launch_shape.sh "$TAG" lognormal 1.0 1e-4 "$SEEDS" 0.5 "$NG" \
    >> "$LOG" 2>&1 &
pid_ln1=$!

echo "shape_rerun: 29 workers launched $(date -Is)" | tee -a "$LOG"

# bidisperse (18 cells, cheapest, never ran) starts once the identity pair and the
# short lognormal alpha-1 set are done, at 9 workers so the box stays near 1 per core.
rc=0
wait "$pid_ident" || { echo "identity launch exited nonzero" | tee -a "$LOG"; rc=1; }
wait "$pid_ln1"   || { echo "lognormal a1 launch exited nonzero" | tee -a "$LOG"; rc=1; }
echo "shape_rerun: identity + lognormal a1 done $(date -Is); starting bidisperse" | tee -a "$LOG"
NSH=9 bash scripts/launch_shape.sh "$TAG" bidisperse 1.0,2.0 1e-3,3e-4,1e-4 "$SEEDS" 0.5 "$NG" \
    >> "$LOG" 2>&1 &
pid_bi=$!

for p in "$pid_ln2" "$pid_un2" "$pid_un1" "$pid_bi"; do
  wait "$p" || { echo "launch $p exited nonzero" | tee -a "$LOG"; rc=1; }
done
echo "--- records under tag $TAG ---" | tee -a "$LOG"
for fam in lognormal uniform bidisperse; do
  n=$(cat ~/pj/out/phase2_${TAG}_"$fam"_shard*.jsonl 2>/dev/null | wc -l)
  echo "  $fam: $n" | tee -a "$LOG"
done
echo "SHAPE_RERUN DONE rc=$rc $(date -Is)" | tee -a "$LOG"
exit $rc
