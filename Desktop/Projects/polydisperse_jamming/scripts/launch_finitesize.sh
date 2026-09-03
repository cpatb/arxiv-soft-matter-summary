#!/bin/bash
# Matched-N family comparison + a lognormal convergence rung.
#
# WHY THIS RUN EXISTS
# -------------------
# The report's actionable claim (REPORT 6a) is that the lognormal family keeps climbing
# to phi_J ~ 0.637 on its tail while the zero-skew families saturate at ~0.596. Two facts
# discovered 2026-09-03 put that claim at risk:
#
#   1. The code's own finite-size flag (d_max/L > 0.25) fires on 100% of lognormal cells
#      at CV = 0.4 and 0.5, and on 0% of uniform and bidisperse cells. The risk is
#      confined to exactly the family that carries the effect -- the long tail that
#      raises phi_J is the same tail that puts the largest grain against the box.
#   2. Lognormal CV = 0.5 moves phi 0.6668 (N = 1500) -> 0.6436 (N = 3000), about 15 seed
#      scatters, and the flag STILL fires at N = 3000 (d_max/L = 0.268).
#
# And the comparison as it stands is CROSS-N: lognormal has N = 3000 data, the zero-skew
# families only N = 1500. Comparing them is comparing two different system sizes.
#
# So this run does two things:
#   A. MATCHED N = 3000 for all three families at CV 0.4 and 0.5, so the contrast is read
#      at one system size.
#   B. N = 6000 for the lognormal, where d_max/L should fall to ~0.21 and clear the flag,
#      to see whether phi_J is still falling or has converged.
#
# PRE-REGISTERED READING (fixed before the numbers were seen):
#   Let D(N) = phi_J(lognormal) - phi_J(uniform) at matched N and CV = 0.5.
#   D(1500) = 0.0428.
#   * If D(3000) and D(6000) stay within ~30% of D(1500), the skewness effect is real and
#     the finite-size flag is a caution rather than a refutation.
#   * If D falls monotonically toward zero, the "lognormal climbs on its tail" claim is a
#     finite-size artefact and 6a must be retracted to "no resolved family difference".
#   * If D shrinks but plateaus well above zero, the effect is real and its MAGNITUDE was
#     overstated at N = 1500; report the converged value.
#
# Cost note: N = 6000 cells are ~4x an N = 1500 cell. Run last so the matched-N block
# lands first and is usable on its own.
#
# SCOPE TRIMMED 2026-09-03 after measuring the cost: N = 3000 cells run 40-60 min each, so
# the original 5-block plan (2 CVs, 3 families at N=3000, 2 at N=6000) needed ~10 h. The
# decisive comparison is the MATCHED-N contrast at the CV where the effect is largest, so
# this now runs CV = 0.5 only at N = 3000 for all three families. CV = 0.4 keeps its
# existing N = 1500 data, and the N = 6000 convergence rung is NOT run -- recorded as not
# done rather than quietly dropped.
set -u
NSHARDS="${1:-11}"
CVS="${2:-0.5}"
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out
cp "$SRC"/*.py ~/pj/ 2>/dev/null    # NOT scripts/: would rewrite this file mid-run
find ~/pj -maxdepth 2 -name "*.py" -o -maxdepth 2 -name "*.sh" | xargs -r sed -i 's/\r$//'
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

run_block () {           # family  N  tag
  local fam="$1" n="$2" tag="$3"
  echo "=== block: family=$fam N=$n tag=$tag ==="
  local pids=()
  for i in $(seq 0 $((NSHARDS - 1))); do
    LOG=~/pj/out/fs_${tag}_w$(printf "%02d" "$i").log
    stdbuf -oL python3 -u phase2_phiJ.py \
        --family "$fam" --n "$n" --cvs "$CVS" --seeds 1,2,3 --mus 0.3 \
        --tag "$tag" --shard "$i" --nshards "$NSHARDS" > "$LOG" 2>&1 &
    pids+=($!)
  done
  local fail=0
  for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
  echo "=== block $tag done, $fail nonzero exits ==="
}

# A. matched N = 3000, all three families
run_block uniform    3000 n3k_uniform
run_block bidisperse 3000 n3k_bidisperse
run_block lognormal  3000 n3k_lognormal

# B. NOT RUN this session: the N = 6000 convergence rung. See the scope note above.
# run_block lognormal  6000 n6k_lognormal

echo "=== ALL FINITE-SIZE BLOCKS FINISHED ==="
