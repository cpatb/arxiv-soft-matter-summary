#!/bin/bash
# Pull the fabric campaign's RAW shards back from WSL. From PowerShell, by file path:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/pull_fabric.sh
#
# RAW RESULTS ONLY. `scripts/sync_back.sh` copies *.csv/*.txt as well and has already
# clobbered fifteen committed derived artifacts once (lessons.md 2026-09-08); it must
# never be used. Derived tables are RE-DERIVED on the Windows side, never copied.
#
# THE SHARDS LAND IN out/fabric/, NOT out/. `out/phase2*.jsonl` is a SHARED GLOB read by
# ten scripts (analyse_phase2, analyse_finite_size, grade_capfix, measure_z_iso, ...) and
# measure_z_iso.py already carries defect D1 about two campaigns writing into it. These
# cells are a SINGLE p_hat rung with seeds 4-12 that no earlier campaign ran, so dropping
# them into the shared glob would hand every phi_J ladder fit one-point "ladders" and
# would collide on keys with the seed 1-3 cells already there. A separate directory makes
# that impossible rather than merely unlikely.
set -u
DST=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/out/fabric
mkdir -p "$DST"
n=0
for f in ~/pj/out/phase2_fab*.jsonl; do
  [ -e "$f" ] || continue
  cp "$f" "$DST/$(basename "$f")" || exit 1
  n=$((n+1))
done
echo "pulled $n shard file(s) to $DST"
echo "--- line counts ---"
wc -l "$DST"/phase2_fab*.jsonl 2>/dev/null | tail -5
