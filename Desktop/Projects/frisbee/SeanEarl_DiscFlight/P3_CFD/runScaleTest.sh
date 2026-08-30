#!/bin/bash
# Reynolds-number-independence test: the SAME gridL1 case solved at three
# reference wind speeds. If the RANS solution is self-similar, U/Uref, k/Uref^2
# and omega/Uref collapse onto one curve and a wind-speed theory can reuse a
# single CFD solve. If it is not, the residual is the error every scaled field
# inherits, and it must be reported.
#
# gridL1 (260k cells) is used deliberately: the test is about the PHYSICS of
# scaling, and it costs a quarter of the production mesh. Whether the residual
# is mesh-dependent is a separate question, flagged rather than assumed.
# NOTE: source the OpenFOAM bashrc BEFORE `set -e`. It returns a nonzero status
# on this install, so sourcing it under `set -e` kills the script silently --
# no output, no log file, exit 1.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
for C in cases/dir_090_gridL1_u10 cases/dir_090_gridL1_u2p5; do
    echo "############ $C  $(date -Is) ############"
    bash runCase.sh "$C" 11
    bash averagePass.sh 11 800 "$C"
    TAG=$(basename "$C")
    python3 sample_disc_altitude.py "$C" --out "flow_fields_${TAG}.npz"
    echo "############ DONE $C  $(date -Is) ############"
done
