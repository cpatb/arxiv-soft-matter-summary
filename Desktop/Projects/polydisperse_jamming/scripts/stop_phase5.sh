#!/bin/bash
# Stop the phase-5 workers. In a SCRIPT FILE, because `pkill -f <pattern>` typed into
# `wsl.exe -e bash -lc '...'` matches the invoking shell's own command line and kills
# the launcher along with (or instead of) the targets. That happened once already.
for p in $(pgrep -f phase5_phase_diagram); do kill "$p" 2>/dev/null; done
sleep 3
echo "phase5 workers remaining: $(pgrep -c -f phase5_phase_diagram || echo 0)"
