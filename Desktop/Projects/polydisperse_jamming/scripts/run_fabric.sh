#!/bin/bash
# Sync code into WSL, then run the fabric campaign. From PowerShell, BY FILE PATH:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_fabric.sh gate
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_fabric.sh prod
#   wsl bash /mnt/c/.../scripts/run_fabric.sh second
#
# Git Bash rewrites a /mnt/c/... argument into a Windows path before wsl.exe sees it
# (lessons.md 2026-09-12), so this must be invoked from PowerShell by file path, never
# inline from Git Bash.
#
# The file list is FIXED and explicit: a glob copies whatever happens to be in the
# directory, which is how a WSL tree silently drifts from the repo.
#
#   gate   -- F1, the bit-identity gate. Three cells that ALREADY EXIST on disk, one per
#             family, re-run under the fabric code. Tag `fabgate`. ~5 min.
#   prod   -- the primary rung, p_hat = 1e-3, 12 seeds, 216 cells. Tag `fab`.
#   second -- the secondary rung, p_hat = 1e-4, 6 seeds, 108 cells. Tag `fab`.
set -u
MODE="${1:?usage: run_fabric.sh gate|prod|second}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py phase2_phiJ.py scripts/launch_fabric.sh"
mkdir -p ~/pj/out ~/pj/scripts
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1

# md5 on the WSL side so the copy can be COMPARED with the Windows one rather than
# assumed equal (lessons.md 2026-09-09). CR is stripped above, so compare against
# `md5sum` of the CR-stripped Windows file, not of the file as committed.
echo "--- md5 in WSL (after CR strip) ---"
for f in $FILES; do echo "  $(md5sum "$f")"; done

CVS=0.0,0.1,0.2,0.3,0.4,0.5
rc=0
case "$MODE" in
  gate)
    # The three cells named in PREREG_FABRIC_2026-09-18.md F1. Each is a SINGLE cell, so
    # 20 workers would be 17 idle processes; NSH is dropped to 1 per family and the three
    # families run in parallel as three invocations.
    # The gate is always run FRESH. The runner is idempotent by key, so a leftover
    # fabgate shard from an earlier code version would be SKIPPED and the gate would
    # then report the old build's numbers as the new build's -- the stale-artifact
    # failure mode (lessons.md 2026-09-14: read the exit code, not the artifact).
    rm -f ~/pj/out/phase2_fabgate*.jsonl
    export NSH=1
    bash scripts/launch_fabric.sh fabgate bidisperse 1e-3 1 0.3 & p1=$!
    bash scripts/launch_fabric.sh fabgate uniform    1e-3 1 0.3 & p2=$!
    bash scripts/launch_fabric.sh fabgate lognormal  1e-3 1 0.4 & p3=$!
    for p in $p1 $p2 $p3; do wait "$p" || rc=1; done
    ;;
  prod)
    for fam in lognormal uniform bidisperse; do
      bash scripts/launch_fabric.sh fab "$fam" 1e-3 1,2,3,4,5,6,7,8,9,10,11,12 "$CVS" || rc=1
    done
    ;;
  second)
    for fam in lognormal uniform bidisperse; do
      bash scripts/launch_fabric.sh fab "$fam" 1e-4 1,2,3,4,5,6 "$CVS" || rc=1
    done
    ;;
  *)
    echo "unknown mode $MODE"; exit 2;;
esac
echo "RUN_FABRIC DONE mode=$MODE rc=$rc"
exit $rc
