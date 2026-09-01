#!/bin/bash
# Copy results OUT of WSL back to the Windows project tree.
# Code flows Windows -> WSL (run_logged.sh / launch_phase2.sh); results flow back here.
set -u
DST=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p "$DST/out" "$DST/figures"
cp -f ~/pj/out/*.jsonl ~/pj/out/*.csv ~/pj/out/*.txt ~/pj/out/*.log "$DST/out/" 2>/dev/null
cp -f ~/pj/figures/* "$DST/figures/" 2>/dev/null
echo "synced: $(ls -1 "$DST/out" 2>/dev/null | wc -l) files in out/, $(ls -1 "$DST/figures" 2>/dev/null | wc -l) in figures/"
