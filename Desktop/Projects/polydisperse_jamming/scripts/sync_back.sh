#!/bin/bash
# Copy RAW RESULTS out of WSL back into the Windows project tree.
# Code flows Windows -> WSL (run_logged.sh / launch_phase*.sh); results flow back here.
#
# 🔴 2026-09-08: this script used to copy *.csv and *.txt as well, and doing so
# CLOBBERED FIFTEEN COMMITTED ANALYSIS OUTPUTS with stale WSL copies -- including
# out/phase2_summary.csv, whose WSL version predates the 2026-09-08 N-grouping fix and
# is missing the `n_requested` column entirely. Nothing failed loudly; it was caught
# only because test_shear.py reads that file and raised KeyError. A summary table is a
# DERIVED artifact and belongs to whichever side last ran the analysis, so this script
# now moves only the things WSL is the source of truth for: raw jsonl records and run
# logs. Re-derive summaries by running the analysers, never by copying them.
set -u
DST=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p "$DST/out" "$DST/figures"
cp -f ~/pj/out/*.jsonl "$DST/out/" 2>/dev/null
cp -f ~/pj/out/*.log   "$DST/out/" 2>/dev/null
cp -f ~/pj/figures/*   "$DST/figures/" 2>/dev/null
echo "synced RAW results only (jsonl + log): $(ls -1 "$DST"/out/*.jsonl 2>/dev/null | wc -l) jsonl in out/"
echo "NOTE: *.csv and *.txt are DERIVED and are deliberately not copied -- run the analysers."
