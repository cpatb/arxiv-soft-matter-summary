#!/usr/bin/env bash
# extend_r2.sh -- restart the refine = 2 case and let it reach its own residual criterion.
#
# WHY. gci.py --family coarse reports that the refine = 2 mesh, which is the paper's fine
# grid and the one every GCI is quoted on, stopped at the 3000-iteration CAP with a worst
# final residual of 3e-5 against a 1e-5 target -- while the refine = 1 and refine = 0.5
# meshes both satisfied residualControl and stopped on their own. So of the three grids in
# the study, the one carrying the answer is the only one that is not iteratively converged.
#
# That matters more than it sounds. A GCI attributes ALL of the difference between grids to
# discretisation. If part of the refine=2 vs refine=1 difference is unconverged iteration,
# then the paper's "C_P changes 12.3 % between grids" and this study's GCI_fine = 6.06 % are
# both measuring the iteration count as well as the mesh. The test is direct: restart from
# the written time-3000 fields, let it run to the criterion, and re-fit. If C_P and
# tau_w(15 mm) do not move, iterative error is negligible and the GCI stands. If they move,
# the grid study has to be redone with a converged fine grid.
#
# Restarting from t = 3000 rather than from scratch is the cheap part: the fields are
# already close, so the residuals have only to fall the last half-decade.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="$HERE/case_Re10000_H6_r2"
TAG="case_Re10000_H6_r2_ext"
WORK="$HOME/bidet_gci/$TAG"
NPROC="${1:-4}"
ENDTIME="${2:-15000}"

rm -rf "$WORK"
mkdir -p "$HOME/bidet_gci"
cp -r "$SRC" "$WORK"
cd "$WORK"
rm -rf processor* log.*

sed -i "s/^startFrom .*/startFrom       latestTime;/" system/controlDict
sed -i "s/^endTime .*/endTime         $ENDTIME;/" system/controlDict
grep -E "^(startFrom|endTime)" system/controlDict

cat > system/decomposeParDict <<EOF
FoamFile { version 2.0; format ascii; class dictionary; object decomposeParDict; }
numberOfSubdomains $NPROC;
method          simple;
coeffs { n ($NPROC 1 1); }
EOF

decomposePar -latestTime > log.decomposePar 2>&1
mpirun -np "$NPROC" simpleFoam -parallel > log.simpleFoam 2>&1
reconstructPar -latestTime > log.reconstructPar 2>&1

CONV=$(grep -c "SIMPLE solution converged" log.simpleFoam || true)
LASTT=$(grep -E "^Time = " log.simpleFoam | tail -1)
echo "[$TAG] converged_flag=$CONV  $LASTT"

rm -rf processor*
rm -rf "$HERE/$TAG"
cp -r "$WORK" "$HERE/$TAG"

cd "$HERE"
echo "--- refit at the extended time ---"
python3 fit_closure.py "$TAG"
echo "--- the value the study currently uses (t = 3000) ---"
python3 fit_closure.py case_Re10000_H6_r2
