#!/bin/bash
# Full analysis chain: phase2 -> phase3 -> figures -> sync back to Windows.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cp -r "$SRC"/*.py ~/pj/ 2>/dev/null
find ~/pj -maxdepth 1 -name "*.py" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
echo "########## PHASE 2 ##########"; python3 analyse_phase2.py
echo; echo "########## PHASE 3 ##########"; python3 analyse_phase3.py
echo; echo "########## SKEWNESS vs CV ##########"; python3 analyse_skewness.py
echo; echo "########## PHASE 4 ##########"; python3 analyse_phase4.py 2>/dev/null || echo "(phase4 analysis not present / no data)"
echo; echo "########## FIGURES ##########"; python3 make_figures.py --check
bash ~/sb.sh
