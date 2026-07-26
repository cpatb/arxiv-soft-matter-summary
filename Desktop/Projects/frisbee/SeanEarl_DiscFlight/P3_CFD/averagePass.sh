#!/bin/bash
# Second pass: restart each converged-as-far-as-it-gets case and TIME-AVERAGE
# the solution over a further N iterations.
#
# WHY THIS IS NEEDED
# ------------------
# Steady SIMPLE does not reach a fixed point on this geometry. The field sits in
# the separated near-wake of a 45 m building, the physical flow there is
# unsteady, and the solver settles into a limit cycle: disc-altitude velocity
# probes oscillate by 13-19% peak-to-peak with no downward trend, while the
# pressure residual plateaus near 5e-3. Taking the final iterate as "the" mean
# field would therefore report an arbitrary phase of that cycle.
#
# Averaging U, k and omega over the limit cycle gives a genuine mean, which is
# what both the throwability metric and the trajectory coupling actually want.
# The averaging window is reported with the results.
NP=${1:-11}
NAVG=${2:-800}
shift 2 2>/dev/null || shift $# 
source /usr/lib/openfoam/openfoam2606/etc/bashrc
for c in "$@"; do
    [ -d "$c" ] || { echo "skip missing $c"; continue; }
    echo "############ averaging $c ############"
    (
    cd "$c" || exit 1
    LATEST=$(foamListTimes -latestTime 2>/dev/null | tail -1)
    [ -z "$LATEST" ] && { echo "no time dirs in $c"; exit 1; }
    END=$((LATEST + NAVG))
    sed -i "s/^startFrom .*/startFrom       latestTime;/" system/controlDict
    sed -i "s/^endTime .*/endTime         $END;/" system/controlDict
    # Append the averaging function object once.
    if ! grep -q "fieldAverageWindow" system/controlDict; then
        python3 - <<PY
import re
p='system/controlDict'; s=open(p).read()
fo = '''
    fieldAverageWindow
    {
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        timeStart       $LATEST;
        fields
        (
            U     { mean on; prime2Mean off; base time; }
            k     { mean on; prime2Mean off; base time; }
            omega { mean on; prime2Mean off; base time; }
        );
    }
'''
i = s.rstrip().rfind('}')
s = s[:i] + fo + s[i:]
open(p,'w').write(s)
PY
    fi
    echo "  restart from $LATEST, averaging to $END"
    decomposePar -force -fields > log.decomposePar.avg 2>&1
    mpirun -np "$NP" simpleFoam -parallel > log.simpleFoam.avg 2>&1 || true
    reconstructPar -latestTime > log.reconstructPar.avg 2>&1
    ls -d [0-9]* 2>/dev/null | tail -2
    ) || echo "FAILED $c"
done
echo "############ AVERAGING COMPLETE ############"
