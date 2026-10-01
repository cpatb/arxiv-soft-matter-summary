#!/usr/bin/env bash
# Episode 1 restoration, reproducible from one command (Git Bash on Windows):
#   bash episodes/ep01_market_street/restore.sh
# Source: archive.org MarketStreet19064KScan20181016,
#   MarketStreet_4K_to_2K_cropped_higher_contrast.mov (2,649,104,587 bytes)
# Outputs (E:, per the outputs-on-the-data-drive rule):
#   bw_master.mp4     phase-locked 24->16 fps decimation, deflicker am/7, 1440x1080
#   color_master.mp4  DDColor-tiny on the iGPU + flow-warped temporal chroma median
cd "$(dirname "$0")/../.." || exit 9
E=/e/silver_to_color/ep01
SRC=/e/silver_to_color/sources/market_2k_cropped.mov
mkdir -p "$E"
DIFFS=()
[ -f "$E/diffs.npz" ] && DIFFS=(--diffs "$E/diffs.npz")
{
  python -m studio.colorize prepare "$SRC" "$E/bw_master.mp4" "${DIFFS[@]}" &&
  python -m studio.colorize color "$E/bw_master.mp4" "$E/color_master.mp4"
  echo "rc=$?"
} > "$E/restore.log" 2>&1
