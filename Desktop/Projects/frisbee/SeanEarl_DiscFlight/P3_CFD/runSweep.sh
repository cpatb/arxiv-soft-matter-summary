#!/bin/bash
# Run a set of cases SEQUENTIALLY. Each case already uses all cores for the
# solve, so running cases concurrently would only contend for memory bandwidth.
# Usage: bash runSweep.sh [nproc] case1 case2 ...
# NOTE: no `set -u`. The OpenFOAM bashrc references unset variables
# (WM_PROJECT_DIR among them) and aborts under nounset.
NP=${1:?usage: runSweep.sh <nproc> <case> [case...]}
shift
source /usr/lib/openfoam/openfoam2606/etc/bashrc
for c in "$@"; do
    echo "############ $c ############"
    if bash runCase.sh "$c" "$NP"; then
        echo "OK  $c"
    else
        echo "FAILED $c (continuing)"
    fi
done
echo "############ SWEEP COMPLETE ############"
