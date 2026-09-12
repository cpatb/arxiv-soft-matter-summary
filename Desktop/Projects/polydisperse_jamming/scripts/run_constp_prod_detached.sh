#!/bin/bash
# Sync, then launch the production ladder DETACHED and return immediately.
#
# ⚠ A nohup'd WSL process only survives while the DISTRO does, and WSL2 tears the distro
# down shortly after its last session closes -- an earlier attempt today left a
# zero-byte log and no workers. So the caller must keep one session open for the
# duration, e.g. by holding a `wsl bash -lc "sleep 570"` alongside and renewing it.
# Use this variant when the run is longer than one foreground call can hold; use
# run_constp_prod.sh (foreground) otherwise.
#
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_constp_prod_detached.sh
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="shear.py phase9_constp.py grade_constp.py grade_constp_prod.py
scripts/launch_constp_prod.sh"
mkdir -p ~/pj/out
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
setsid nohup bash scripts/launch_constp_prod.sh \
    > ~/pj/out/constp_prod_launch.log 2>&1 < /dev/null &
echo "launched pid $!  (setsid + nohup; hold a session open)"
