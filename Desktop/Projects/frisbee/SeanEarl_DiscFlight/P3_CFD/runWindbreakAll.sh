#!/bin/bash
# The windbreak design sweep: three architectures against the same baseline.
#
# Run at PRODUCTION resolution, not on the cheap grid. The grid-independence
# study (GRID_INDEPENDENCE_2026-08-30.md) found that the coarse meshes do not
# resolve Campus Towers -- 3 to 6 cells across its 24 m narrow side against a
# COST 732 minimum of 10 -- and the whole point of this study is what happens in
# that building's windward stagnation region. A design sweep on a mesh that
# cannot see the building would be measuring the mesh.
#
#   W1  shelterbelt        x = +21 m, H = 5 m,  k_r = 3.00  (beta ~ 0.5)
#   W2  tower deflector    x = +52 m, H = 10 m, k_r = 3.00
#   W3  low open spoiler   x = +30 m, H = 3 m,  k_r = 1.78  (beta ~ 0.6)
#
# All three at theta_wind = 270, the heading whose east/west sigma_u ratio is
# 1.80. Whichever wins is then re-run at 90 and 210 to check it does not make
# those worse -- a barrier that fixes one wind direction and ruins the other two
# is not a fix.
source /usr/lib/openfoam/openfoam2606/etc/bashrc
set -e
HERE=$(cd "$(dirname "$0")" && pwd); cd "$HERE"
NP=${1:-10}
for C in cases/dir_270_W1 cases/dir_270_W2 cases/dir_270_W3; do
    bash runWindbreakCase.sh "$C" "$NP" 800
done
echo "ALL WINDBREAK CASES DONE"
