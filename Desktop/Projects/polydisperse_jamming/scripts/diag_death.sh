#!/bin/bash
# Print the solver lines of err.log that explain a death (STOP / assert / abort / signal / nan).
#   wsl bash .../scripts/diag_death.sh <campaign_dir_name> <label>
d=/home/cpb/work/campaigns/$1/$2
grep -a -n -iE "STOP|assert|abort|signal|nan|segm|poisson|fail|killed|timeout" "$d/err.log" | grep -v -i "ffmpeg\|libav\|Error opening" | tail -8
echo "-- out.log tail"; tail -3 "$d/out.log" 2>/dev/null
