#!/bin/bash
# THE ITERATION-CAP REPAIR -- PREREG_DST_WINDOW_2026-09-09.md §5.
#
# 82 isotropic cells on disk are marked non-converged and EVERY one stopped at exactly
# MAX_ITER = 200000, i.e. slow rather than wrong. The censoring is NOT flat in the axis
# this project measures, and it is not symmetric between the two arms of the
# frictionless-minus-frictional window:
#
#     lognormal mu_p = 0.0 : converged 0.92 (CV 0) -> 0.33 (CV 0.5)
#     lognormal mu_p = 0.3 : converged 1.00 (CV 0) -> 0.75 (CV 0.5)
#
# and 18 of the 82 met the physical bar (unbalanced < 0.02) and were rejected ONLY on
# the iteration count, 12 of those in the starved frictionless arm.
#
# This re-runs the affected region at PJ_MAX_ITER = 1000000. `phase2_phiJ.cell_key`
# appends `_it1000000` when the cap is non-default, so these records CANNOT merge with
# the capped ones they exist to replace, and `packing` writes `max_iter` and `tol_met`
# into every record.
#
# ⚠ CV 0.0 IS THE CONTROL, not padding. Those cells converged at ~11k iters, three
# orders below either cap, so raising the cap cannot change them: phi must come back
# bit-identical. If it does not, the protocol is not deterministic under this knob and
# nothing else in the run is interpretable.
#
# Both frictions are repaired at the CV levels where either arm is censored. Repairing
# one arm of a difference and not the other is the error this run exists to fix, not one
# to repeat.
#
# Tag `capfix`; glob phase2_capfix*.jsonl matches no existing analysis.
set -u
# 10 per arm and TWO arms launched together = 20 workers on 22 cores.
NSHARDS=10
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
echo "PJ_MAX_ITER=$PJ_MAX_ITER"

pids=()

run_arm () {
  fam="$1"
  mus="$2"
  cvs="$3"
  nick="$4"
  for i in $(seq 0 $((NSHARDS - 1))); do
    LOG=~/pj/out/capfix_${nick}_w$(printf "%02d" "$i").log
    stdbuf -oL python3 -u phase2_phiJ.py \
        --shard "$i" --nshards "$NSHARDS" --family "$fam" \
        --mus "$mus" --cvs "$cvs" --seeds 1,2,3 --tag "capfix_${nick}" \
        > "$LOG" 2>&1 &
    pids+=($!)
  done
  echo "capfix $nick: $NSHARDS workers, family=$fam mus=$mus cvs=$cvs"
}

# frictionless arm: the control (0.0) plus every censored level
run_arm lognormal 0.0 0.0,0.3,0.4,0.5 lg0
run_arm uniform   0.0 0.0,0.3,0.4,0.5 un0

fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
pids=()

# Second wave. Three arms of 10 workers would be 30 on 22 cores, so NSHARDS drops to 7.
NSHARDS=7

# frictional arm, only where it is censored (CV 0.4, 0.5), run second so it does not
# compete with the frictionless one, which carries the claim.
run_arm lognormal 0.3 0.4,0.5 lg3
run_arm uniform   0.3 0.4,0.5 un3

# mu_p = 1.0 across CV -- PREREG_DST_WINDOW_2026-09-09.md §1, "the cheapest remaining
# gap". phi_J SATURATES in friction (the mu_p 1 -> 3 step at CV = 0 is inside seed
# scatter), so mu_p = 1 is Wyart-Cates' phi_mu as this protocol can measure it, and the
# window read there is 0.0741 against 0.0586 at the production friction -- a factor 1.27
# on the very quantity the window argument divides by. It exists ONLY at CV = 0 today, so
# `dphi_mu/dCV` at the saturated friction is unmeasured and the traverse currently pairs
# a saturated width with a mu_p = 0.3 slope. CV 0/0.2/0.5 matches the shear ladder's own
# levels. Run at the SAME raised cap as the frictionless arm, so the two legs of the
# window are never differenced across a protocol change.
run_arm lognormal 1.0 0.0,0.2,0.5 lg1

for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== capfix DONE, $fail nonzero exits ==="
