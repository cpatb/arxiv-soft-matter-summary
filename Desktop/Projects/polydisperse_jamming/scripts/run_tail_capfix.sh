#!/bin/bash
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py phase11_tail.py scripts/launch_tail_capfix.sh"
mkdir -p ~/pj/out
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
exec bash scripts/launch_tail_capfix.sh
