#!/bin/bash
# Run the remaining two headings at gridL4 (fieldBox level 4, 1.0 x 0.25 m), so the
# heading COMPARISON -- the paper's one resolved directional claim -- can be tested at a
# resolution where every major building clears the COST 732 guideline.
#
# Sequential: each solve already uses all 11 ranks, so running them concurrently would
# only contend for memory bandwidth.
#
# Source BEFORE `set -e`: OpenFOAM's bashrc returns non-zero on this install and would
# abort the script instantly and silently under `set -e`.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

NP=11
NAVG=800

for H in 210 270; do
    CASE="cases/dir_${H}_gridL4"
    echo "############################################################"
    echo "# heading $H at gridL4   $(date -Is)"
    echo "############################################################"
    T0=$SECONDS

    ( cd "$CASE" && blockMesh > log.blockMesh 2>&1 )
    ( cd "$CASE" && surfaceFeatureExtract > log.sfe 2>&1 )
    ( cd "$CASE" && snappyHexMesh -overwrite > log.snappy 2>&1 )
    ( cd "$CASE" && rm -rf 0 && cp -r 0.orig 0 )
    ( cd "$CASE" && checkMesh > log.checkMesh 2>&1 || true )
    echo "  mesh: $(grep -E '^    cells:' "$CASE/log.checkMesh" | head -1)"

    ( cd "$CASE" && decomposePar -force > log.decomposePar 2>&1 )
    ( cd "$CASE" && mpirun -np "$NP" simpleFoam -parallel > log.simpleFoam 2>&1 || true )

    LAST=$(grep -E '^Time = ' "$CASE/log.simpleFoam" | tail -1 | awk '{print $3}')
    echo "  solve ended at Time = $LAST"
    if [ "$LAST" != "1400" ] || ! grep -q "^End" "$CASE/log.simpleFoam"; then
        echo "  ERROR: heading $H did not reach 1400 cleanly; skipping (not averaging a short run)"
        continue
    fi

    ( cd "$CASE" && reconstructPar -latestTime > log.reconstructPar 2>&1 )
    bash averagePass.sh "$NP" "$NAVG" "$CASE"
    python3 sample_disc_altitude.py "$CASE" --out "flow_fields_dir_${H}_gridL4.npz"
    echo "  heading $H done in $((SECONDS-T0)) s"
done

echo "ALL HEADINGS DONE $(date -Is)"
