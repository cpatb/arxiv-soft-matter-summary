#!/bin/bash
# Compare the CR-stripped repo shear.py/packing.py with the copies the WSL campaigns ran.
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
for f in shear.py packing.py phase9_constp.py phase2_shape.py; do
  a=$(tr -d '\r' < "$SRC/$f" | md5sum | cut -c1-12)
  b=$(md5sum < ~/pj/$f 2>/dev/null | cut -c1-12)
  echo "$f repo $a wsl ${b:-missing} $([ "$a" = "$b" ] && echo SAME || echo DIFF)"
done
