#!/bin/bash
# Copy the dsttype launcher into ~/work (CR-stripped), verify md5 both sides, then exec it,
# so the file bash is reading can never be edited under it.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/dsttype_2026-09-26.sh
DST=/home/cpb/work/dsttype_2026-09-26.sh
tr -d '\r' < "$SRC" > "$DST"
a=$(tr -d '\r' < "$SRC" | md5sum | cut -c1-32); b=$(md5sum < "$DST" | cut -c1-32)
[ "$a" = "$b" ] || { echo "MD5 MISMATCH"; exit 1; }
bash -n "$DST" || { echo "SYNTAX"; exit 1; }
echo "launcher md5 $b"
exec bash "$DST"
