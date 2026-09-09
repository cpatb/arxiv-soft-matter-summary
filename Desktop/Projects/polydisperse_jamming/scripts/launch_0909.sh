#!/bin/bash
# Two work packages for 2026-09-09, run sequentially so neither oversubscribes the box.
#
#   A  LADDER EXTENSION (NEXT_STEPS_2026-09-09 task 0.2). Two more rungs, +0.10 and
#      +0.12, at cv = 0 and 0.2, so that applying the quality gates leaves more than
#      one seed with enough rungs for the four-parameter E2 fit. Tagged `prod2ext`,
#      which the default analysis glob `phase8_prod2*.jsonl` matches DELIBERATELY:
#      protocol, gamma_max, shear_steps and N are identical to prod2, and those four
#      are all in analyse_phase8.group_key, so the new rungs join the existing ladders
#      rather than forming a new group.
#
#   B  RATE-MATCHED CONTROL, found 2026-09-09 and not in any plan. `gdot` is NOT
#      constant across the CV ladder: gdot = gamma_max/(shear_steps*dt) and dt is
#      DT_SAFETY*PWaveTimeStep(), which is set by the SMALLEST particle. Measured on
#      prod2: gdot = 3.000e-05 at cv = 0 and 1.437e-04 at cv = 0.5, a factor 4.8.
#      So the CV contrast is confounded with a shear-rate contrast. This re-runs the
#      cv = 0 and cv = 0.2 ladders with shear_steps chosen per cell so that gdot lands
#      on the cv = 0.5 value, 1.430e-04. Every record stores its realised gdot, so the
#      match is CHECKED afterwards, never assumed.
#
#      steps = gamma_max / (gdot_target * dt), dt read from the prod2 records:
#        cv 0.0, all seeds  dt 0.250000 -> 41958
#        cv 0.2, seed 1     dt 0.121394 -> 86411
#        cv 0.2, seed 2     dt 0.133128 -> 78795
#        cv 0.2, seed 3     dt 0.126717 -> 82782
#      Tag `ratem`, which does NOT match phase8_prod2*, because these cells are a
#      DIFFERENT protocol in the only sense that matters and must not pool.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
L="$SRC/scripts/launch_phase8.sh"

echo "=== A: ladder extension, cv 0.0 and 0.2, offsets +0.10 +0.12 ==="
bash "$L" 12 --cvs 0.0,0.2 --seeds 1,2,3 --offsets 0.10,0.12 --tag prod2ext

echo "=== B: rate-matched control, cv 0.0 (steps 41958) ==="
bash "$L" 15 --cvs 0.0 --seeds 1,2,3 --shear-steps 41958 --tag ratem_c0

echo "=== B: rate-matched control, cv 0.2 seed 1 (steps 86411) ==="
bash "$L" 5 --cvs 0.2 --seeds 1 --shear-steps 86411 --tag ratem_c2s1
echo "=== B: rate-matched control, cv 0.2 seed 2 (steps 78795) ==="
bash "$L" 5 --cvs 0.2 --seeds 2 --shear-steps 78795 --tag ratem_c2s2
echo "=== B: rate-matched control, cv 0.2 seed 3 (steps 82782) ==="
bash "$L" 5 --cvs 0.2 --seeds 3 --shear-steps 82782 --tag ratem_c2s3

echo "=== ALL DONE 0909 ==="
