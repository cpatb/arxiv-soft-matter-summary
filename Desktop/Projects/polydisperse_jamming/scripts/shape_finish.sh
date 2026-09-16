#!/bin/bash
# shape_finish.sh -- wait for the in-flight family, then run the remaining two.
#
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/shape_finish.sh
#
# Written as a FILE and invoked by file path from PowerShell. An inline
# `wsl bash -lc "until [ $(pgrep -c ...) -eq 0 ]; do ..."` loses its command
# substitution to quote mangling and becomes `until [ -eq 0 ]`, which never fires
# (lessons.md 2026-09-01: 18 shells polled for 14 hours and none could have fired).
# The loop below ECHOES the value it actually read, for the same reason.
#
# Families run ONE AT A TIME at 20 workers. Running two families concurrently puts 40
# DEM workers on 22 cores and halves the throughput of both, which is exactly what
# happened at 07:30 today and took 20 minutes to notice.
set -u
cd ~/pj || exit 1

# CV = 0.5 ONLY for the second and third families.
#
# At CV = 0 every PSD family is the SAME delta distribution: `diag_shape_cv0.py` shows
# `psd.make_psd` returns bit-identical radii for all three (sum 2700.0, min = max = 0.5)
# and the nine CV = 0 cells present in two families returned phi bit-identically
# (max |dphi| = 0.000e+00), while at CV = 0.5 the three families genuinely differ --
# which is the discriminating half, without which the identity would only mean the
# family flag does nothing. So running CV = 0 once per family is running ONE cell three
# times, and a third of the grid was redundant.
#
# The grader ALIASES the CV = 0 rung across families under gate S5 and prints that it
# did. The nine already-run duplicates are the gate's evidence and are kept.
CVS=0.5
ALPHAS=1.0,2.0
SEEDS=1,2,3
PHATS=1e-3,3e-4,1e-4
NG=5400

echo "shape_finish: waiting for the in-flight workers"
i=0
while true; do
  n=$(ps -eo args | grep -c '[p]hase2_shape.py')
  echo "  poll $i: $n worker(s) in flight"
  if [ "$n" -eq 0 ]; then break; fi
  i=$((i+1))
  sleep 60
done
echo "shape_finish: in-flight set is empty; starting the remaining families"

rc=0
for fam in uniform bidisperse; do
  n_before=$(cat ~/pj/out/phase2_shp_"$fam"_shard*.jsonl 2>/dev/null | wc -l)
  echo "=== $fam: $n_before record(s) already on disk (the runner resumes by key) ==="
  NSH=20 bash scripts/launch_shape.sh shp "$fam" "$ALPHAS" "$PHATS" "$SEEDS" "$CVS" "$NG" || rc=1
  n_after=$(cat ~/pj/out/phase2_shp_"$fam"_shard*.jsonl 2>/dev/null | wc -l)
  echo "=== $fam done: $n_before -> $n_after records ==="
done

echo "--- final counts ---"
for fam in lognormal uniform bidisperse; do
  n=$(cat ~/pj/out/phase2_shp_"$fam"_shard*.jsonl 2>/dev/null | wc -l)
  echo "  $fam: $n"
done
echo "SHAPE_FINISH DONE rc=$rc"
exit $rc
