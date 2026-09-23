#!/bin/bash
# Overnight status: record counts per campaign (no values).
date
echo "kmu  $(cat ~/pj/out/phase2_kmu_*.jsonl 2>/dev/null | wc -l)/36"
echo "kmix $(cat ~/pj/out/phase2_kmix_*.jsonl 2>/dev/null | wc -l)/36"
echo "kmermap $(ls ~/work/campaigns/kmermap_2026-09-23/*/RUN_STATUS.json 2>/dev/null | wc -l)/18"
echo "kmerdrop2 $(ls ~/work/campaigns/kmerdrop2_2026-09-22/*/RUN_STATUS.json 2>/dev/null | wc -l)/2"
echo "failed DEM cells: $(grep -h 'CELL FAILED' ~/pj/out/kmu_w*.log ~/pj/out/kmix_w*.log 2>/dev/null | wc -l)"
uptime
