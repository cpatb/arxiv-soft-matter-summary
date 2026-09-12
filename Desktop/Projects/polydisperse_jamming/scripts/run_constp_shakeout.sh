#!/bin/bash
# Sync the Windows tree into ~/pj, then run the shakeout IN THE FOREGROUND.
# Foreground on purpose: see the nohup note in launch_constp_shakeout.sh.
#
# Invoke from PowerShell, never from Git Bash: Git Bash rewrites a /mnt/c/... argument
# into a Windows path before wsl.exe sees it (measured 2026-09-12: it became
# "C:/Program Files/Git/mnt/c/...".) And never through `wsl bash -lc "..."` with double
# quotes, which eats every $.
#
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_constp_shakeout.sh
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="shear.py phase9_constp.py check_prepare_refactor.py measure_z_iso.py
analyse_phase8.py packing.py measure.py psd.py shear_pre_refactor_src.txt
scripts/launch_constp_shakeout.sh scripts/run_prepare_check.sh"
mkdir -p ~/pj/out
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
exec bash scripts/launch_constp_shakeout.sh
