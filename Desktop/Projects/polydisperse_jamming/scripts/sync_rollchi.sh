#!/bin/bash
# Pull rollchi raw records (jsonl) and clump contact dumps back to the repo. Raw results only.
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p "$SRC/out/rollchi" "$SRC/out/achi_dump"
cp ~/pj/out/phase2_rchi*_lognormal*.jsonl ~/pj/out/phase2_achi*_lognormal*.jsonl "$SRC/out/rollchi/" 2>/dev/null
cp ~/pj/out/achi_dump/*.npz "$SRC/out/achi_dump/" 2>/dev/null
echo "synced $(cat "$SRC"/out/rollchi/*.jsonl | wc -l) records, $(ls "$SRC"/out/achi_dump | wc -l) dumps"
