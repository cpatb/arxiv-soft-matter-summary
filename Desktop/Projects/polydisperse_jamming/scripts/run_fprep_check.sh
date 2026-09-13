#!/bin/bash
# Bit-identity of the DEFAULT path after the 2026-09-12 `mu_prep` / `n_samples` change,
# plus the knob's non-default arm (mu_prep = 0). ~4 min, foreground. From PowerShell:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_fprep_check.sh
# The sidecar `shear_head_src.txt` is `git show HEAD:...shear.py` made on the Windows
# side, because ~/pj is a file copy, not a checkout.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="shear.py check_prepare_refactor.py shear_head_src.txt"
mkdir -p ~/pj/out
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
stdbuf -oL python3 -u check_prepare_refactor.py --sidecar shear_head_src.txt \
    --old-has-prepare --fprep-smoke 2>&1 | tee ~/pj/out/fprep_check.log
echo "exit: ${PIPESTATUS[0]}"
