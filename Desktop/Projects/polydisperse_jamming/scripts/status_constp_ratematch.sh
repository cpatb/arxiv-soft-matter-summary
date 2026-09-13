#!/bin/bash
# Status of the rate-matched ladder. Invoke from PowerShell by path:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/status_constp_ratematch.sh
echo "workers alive: $(pgrep -cf phase9_constp)"
for t in c00 c20 c50; do
  n=$(cat ~/pj/out/phase9_rm_${t}_shard*.jsonl 2>/dev/null | wc -l)
  echo "$t records: $n"
done
for f in ~/pj/out/constp_rm_c50_w*.log; do
  echo "== $f"; tail -n 2 "$f"
done
