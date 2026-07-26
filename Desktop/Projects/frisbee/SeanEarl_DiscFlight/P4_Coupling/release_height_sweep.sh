#!/bin/bash
# Paired release-height experiment: how much does turbulence cost you as a
# function of how high you throw?
#
# Two arms per height, sharing seed and therefore sharing release-parameter
# draws:
#   gust arm  -- RANS mean + synthetic correlated gusts scaled by local k
#   mean arm  -- RANS mean only (--no-gusts)
#
# The mean arm is the control for a confound that would otherwise dominate:
# a disc released higher stays airborne longer and disperses more even in
# perfectly smooth air. The DIFFERENCE between the arms is the part of the
# scatter attributable to turbulence, which is the uncorrectable part.
#
# Uses vertical_fields.npz (z up to 100 m), not flow_fields.npz (z <= 6 m):
# CFDField clamps rather than extrapolates, so on the 6 m field every release
# above 6 m would silently reuse the 6 m wind and the trend would be an
# artefact of the sampling ceiling.
cd "$(dirname "$0")"
set -e

FIELD=../P3_CFD/vertical_fields.npz
HEIGHTS="1.0 2.0 3.0 4.0 6.0 10.0"
NTHROWS=${1:-64}
JOBS=${2:-10}
OUTDIR=height_sweep
mkdir -p $OUTDIR

for z in $HEIGHTS; do
  for arm in gust mean; do
    flag=""
    [ "$arm" = "mean" ] && flag="--no-gusts"
    out=$OUTDIR/h${z}_${arm}.npz
    if [ -f "$out" ]; then echo "skip $out (exists)"; continue; fi
    echo "#### release z = $z m, arm = $arm"
    python3 monte_carlo.py \
        --field $FIELD \
        --launch-grid -20 30 25 -5 10 15 \
        --launch-z "$z" \
        --headings 90 210 270 \
        --throws right_backhand right_forehand \
        --throw-bearings 0 90 180 270 \
        --n-throws "$NTHROWS" \
        --jobs "$JOBS" \
        $flag \
        --out "$out" 2>&1 | tail -3
  done
done
echo "#### HEIGHT SWEEP COMPLETE ####"
