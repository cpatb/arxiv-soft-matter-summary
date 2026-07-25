#!/usr/bin/env bash
# run_sweep.sh -- parameter sweep (Re_d x H/d) + grid-independence, run in
# parallel batches. Collects each fit_closure JSON line into sweep_results.jsonl.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
OUT="$HERE/sweep_results.jsonl"
: > "$OUT"
: > "$HERE/sweep_console.log"

MAXJOBS="${MAXJOBS:-6}"

one () {  # Re Hd refine
  bash run_case.sh "$1" "$2" "${3:-1}" >> "$HERE/sweep_console.log" 2>&1
}

# build job list: main sweep + grid-independence at nominal
JOBS=()
for RE in 3000 10000 30000; do
  for HD in 1 6 20 60; do JOBS+=("$RE $HD 1"); done
done
JOBS+=("10000 6 2")

echo "launching ${#JOBS[@]} cases, up to $MAXJOBS concurrent..."
for J in "${JOBS[@]}"; do
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do sleep 2; done
  echo "  start: Re=$J"
  one $J &
done
wait

# collect all fits (nominal already in file? no -- recompute all for consistency)
: > "$OUT"
for d in case_Re*_H*; do
  [ -d "$d/postProcessing" ] && python3 fit_closure.py "$d" >> "$OUT" 2>>"$HERE/sweep_console.log" || true
done
echo "sweep done. $(wc -l < "$OUT") fits -> $OUT"
