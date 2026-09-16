#!/bin/bash
# Sync code into WSL, then run the shape campaign. From PowerShell, BY FILE PATH:
#   wsl bash /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/scripts/run_shape.sh s0
#
# Git Bash rewrites a /mnt/c/... argument into a Windows path before wsl.exe sees it
# (lessons.md 2026-09-12), so this must be invoked from PowerShell by file path, never
# inline from Git Bash.
#
# The file list is FIXED and explicit: a glob copies whatever happens to be in the
# directory, which is how a WSL tree silently drifts from the repo.
#
#   s0     -- gate S0. The three cells named in PREREG_FABRIC_2026-09-18.md F1, re-run
#             at alpha = 1 through the SHAPE runner. Tag `shpgate`. ~5 min.
#   time   -- a timing and clearance probe at the registered N. Writes shards under the
#             tag `shptime` which the grader refuses to admit.
#   prod   -- the registered production grid.
set -u
MODE="${1:?usage: run_shape.sh s0|time|prod}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
FILES="packing.py psd.py measure.py shape.py phase2_phiJ.py phase2_shape.py scripts/launch_shape.sh"
mkdir -p ~/pj/out ~/pj/scripts
for f in $FILES; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1

# md5 on the WSL side so the copy can be COMPARED with the Windows one rather than
# assumed equal (lessons.md 2026-09-09). CR is stripped above, so compare against the
# md5 of the CR-stripped Windows file, not of the file as committed.
echo "--- md5 in WSL (after CR strip) ---"
for f in $FILES; do echo "  $(md5sum "$f")"; done

# THE REGISTERED LADDER, amendment 5.6. The 3e-5 rung was dropped after the first
# production cells showed it does not jam at CV = 0.5 with N = 5400 -- z_all = 0.023 at
# the 200 000-iteration cap, against a converged 4.264 at CV = 0.0 -- so it censored the
# ladder preferentially at one end of the axis under test.
PHATS_FULL=1e-3,3e-4,1e-4
rc=0
case "$MODE" in
  s0)
    # Always run FRESH. The runner is idempotent by key, so a leftover shpgate shard
    # from an earlier code version would be SKIPPED and the gate would report the old
    # build's numbers as the new build's (lessons.md 2026-09-14: read the exit code,
    # not the artifact).
    rm -f ~/pj/out/phase2_shpgate*.jsonl
    export NSH=1
    bash scripts/launch_shape.sh shpgate bidisperse 1.0 1e-3 1 0.3 1500 & p1=$!
    bash scripts/launch_shape.sh shpgate uniform    1.0 1e-3 1 0.3 1500 & p2=$!
    bash scripts/launch_shape.sh shpgate lognormal  1.0 1e-3 1 0.4 1500 & p3=$!
    for p in $p1 $p2 $p3; do wait "$p" || rc=1; done
    ;;
  time)
    # The cost of one cell at each corner that matters, so the grid can be SIZED from a
    # measurement instead of an estimate. N is the value the clearance check demands at
    # the worst cell (lognormal CV 0.5, alpha 2), computed on the Windows side and
    # passed in; the cheap corners are run at the SAME N because N cannot vary across a
    # grid whose whole point is a comparison between its cells.
    NG="${2:?usage: run_shape.sh time N}"
    export NSH=12
    bash scripts/launch_shape.sh shptime lognormal 1.0,1.5,2.0 1e-3,3e-5 1 0.0,0.5 "$NG" || rc=1
    ;;
  prod)
    TAG="${2:?usage: run_shape.sh prod TAG ALPHAS CVS SEEDS N}"
    ALPHAS="${3:?}"; CVS="${4:?}"; SEEDS="${5:?}"; NG="${6:?}"
    for fam in lognormal uniform bidisperse; do
      bash scripts/launch_shape.sh "$TAG" "$fam" "$ALPHAS" "$PHATS_FULL" \
          "$SEEDS" "$CVS" "$NG" || rc=1
    done
    ;;
  *)
    echo "unknown mode $MODE"; exit 2;;
esac
echo "RUN_SHAPE DONE mode=$MODE rc=$rc"
exit $rc
