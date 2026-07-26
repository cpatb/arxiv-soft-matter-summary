#!/bin/bash
# Runs after the main sweep: the matched-roughness ABL validation case, then
# the time-averaging pass over every production case.
cd "$(dirname "$0")"
until grep -q "SWEEP COMPLETE" cases/sweep.log 2>/dev/null; do sleep 20; done
echo "### main sweep finished; starting follow-up ###"
source /usr/lib/openfoam/openfoam2606/etc/bashrc
echo "### 1. matched-roughness ABL homogeneity case ###"
bash runCase.sh cases/dir_090_abl 11 || echo "ABL case FAILED"
echo "### 2. time-averaging pass ###"
bash averagePass.sh 11 800 cases/dir_090 cases/dir_090_coarse cases/dir_210 cases/dir_270
echo "### FOLLOW-UP COMPLETE ###"
