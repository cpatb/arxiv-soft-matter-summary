#!/bin/bash
# shape_kill_stale.sh -- kill the 4-rung shape run and leave the registered 3-rung one.
#
#   wsl bash /mnt/c/.../scripts/shape_kill_stale.sh
#
# WHY THIS EXISTS. `pkill -f phase2_shape.py` kills the WORKERS and not the
# `run_shape.sh` / `launch_shape.sh` parents, and run_shape.sh's `for fam in ...` loop
# continues on failure -- so killing the workers merely advanced the old campaign to its
# next family, and two campaigns then ran side by side at 40 workers on 22 cores. The
# symptom was a 2x slowdown with nothing obviously wrong; `scripts/shape_status.sh`
# found it by printing the distinct `--phats` values in flight.
#
# Kill the PARENTS first, then the workers, or the parents start the next family while
# the workers are dying.
set -u

echo "--- before ---"
ps -eo args | grep -c '[p]hase2_shape.py' | xargs echo "  workers:"
ps -eo args | grep '[p]hase2_shape.py' | grep -o '\-\-phats [0-9e.,-]*' | sort | uniq -c

# The stale campaign is the one whose ladder still carries the 3e-5 rung.
STALE='--phats 1e-3,3e-4,1e-4,3e-5'

# 1. the launchers. Both scripts are shell, so match on the script name.
pkill -f 'run_shape.sh prod' 2>/dev/null
pkill -f 'launch_shape.sh' 2>/dev/null
sleep 1

# 2. the stale workers only, matched on their own ladder argument.
pkill -f "phase2_shape.py.*3e-05" 2>/dev/null
pkill -f "phase2_shape.py.*1e-3,3e-4,1e-4,3e-5" 2>/dev/null
sleep 3

echo "--- after ---"
ps -eo args | grep -c '[p]hase2_shape.py' | xargs echo "  workers:"
ps -eo args | grep '[p]hase2_shape.py' | grep -o '\-\-phats [0-9e.,-]*' | sort | uniq -c
ps -eo args | grep '[p]hase2_shape.py' | grep -o '\-\-family [a-z]*' | sort | uniq -c
echo "NOTE: killing the launchers means the remaining families must be re-invoked; the"
echo "runner is idempotent by cell key, so re-invoking RESUMES."
exit 0
