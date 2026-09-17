#!/bin/bash
# The controls for probe_seed3_crash.sh: the same lines on the cells that RAN.
for c in "lognormal 0.5 1.0 3" "lognormal 0.5 2.0 1" "lognormal 0.5 2.0 2" "uniform 0.5 2.0 3" "bidisperse 0.5 2.0 3" "bidisperse 0.5 2.0 1" "bidisperse 0.5 2.0 2"; do
  bash ~/pj/scripts/probe_seed3_crash.sh $c 2>&1 | grep -v "^    iter" | grep -v "Using python"
  bash ~/pj/scripts/probe_seed3_crash.sh $c 2>&1 | grep "^    iter" | tail -1
done
