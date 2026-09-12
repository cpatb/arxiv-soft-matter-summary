#!/bin/bash
# The _prepare refactor check, in the foreground. ~2 min.
set -u
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
exec stdbuf -oL python3 -u check_prepare_refactor.py
