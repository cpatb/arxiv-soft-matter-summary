#!/usr/bin/env bash
# resid.sh <log> -- last initial-residual per field, and the trend of the worst one.
# Exists because quoting a sed backreference through `wsl bash -lc '...'` mangles it.
LOG="$1"
echo "== last initial residual per field =="
for F in Ux Uy Uz p k omega; do
  V=$(grep "Solving for $F," "$LOG" | tail -1 | awk -F'Initial residual = ' '{print $2}' | awk -F',' '{print $1}')
  [ -n "$V" ] && printf "  %-6s %s\n" "$F" "$V"
done
echo "== k initial residual, every 200th solve =="
grep "Solving for k," "$LOG" | awk -F'Initial residual = ' 'NR%200==0{split($2,a,","); print NR, a[1]}'
echo "== iterations =="
grep -E "^Time = " "$LOG" | tail -1
grep -c "SIMPLE solution converged" "$LOG"
