#!/bin/bash
# Sync the code into ~/pj (CR-stripped), then run the actshape campaign in the foreground.
# Invoke from PowerShell by file path, never inline from Git Bash:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_actshape.sh
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out ~/pj/scripts
for f in shear.py phase9_constp.py packing.py measure.py psd.py scripts/launch_actshape.sh scripts/launch_actshape_fs.sh; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
for f in shear.py phase9_constp.py scripts/launch_actshape.sh scripts/launch_actshape_fs.sh; do
  a=$(tr -d '\r' < "$SRC/$f" | md5sum | cut -c1-32); b=$(md5sum < ~/pj/"$f" | cut -c1-32)
  [ "$a" = "$b" ] || { echo "MD5 MISMATCH $f"; exit 1; }
done
echo "synced, md5 ok"
cd ~/pj || exit 1
exec bash scripts/${1:-launch_actshape.sh}
