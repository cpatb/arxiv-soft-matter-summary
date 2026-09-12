#!/bin/bash
# THE COMMON ABSOLUTE PHI LADDER -- PREREG_COMMON_LADDER_2026-09-09.md.
#
# Every phi_c ladder in this project is centred on that cell's published isotropic
# phi_J, so the fitting window moves with CV by exactly delta phi_J_iso and a contrast
# computed on it can be inherited rather than measured. This run puts all three CV
# levels on the SAME SIX ABSOLUTE phi:
#
#     phi = 0.545, 0.570, 0.595, 0.620, 0.645, 0.670
#
# realised by giving each CV its own --offsets = grid - PHI_J_ISO[cv]. Spacing 0.025 is
# uniform, as --offsets requires, and 5.3x the largest seed scatter in PHI_J_ISO.
#
# Every other knob is held at the centred ladders' setting on purpose (N=1500,
# mu_p=0.3, lognormal, gamma_max=1.5, shear_steps=200000, NOT rate-matched), because
# the comparison is common-window vs centred-window and the two arms must differ in one
# knob only. The rate confound is measured separately in out/rate_matched_verdict.txt.
#
# 3 CV x 3 seeds x 6 rungs = 54 cells. Runs against the copy already in ~/pj (no
# concurrent cp), same pattern as launch_0909b.sh.
#
# ⚠ Each CV gets its OWN tag (commonlad_c00 / _c20 / _c50). The three invocations run
# concurrently and phase8 APPENDS to out/phase8<tag>_shard<i>of<n>.jsonl, so a shared
# tag would have three processes appending to one file; a JSON record here is ~30 kB,
# far above PIPE_BUF, so those appends are not atomic. The analysis glob
# `phase8_commonlad*.jsonl` matches all three.
set -u
NSHARDS=6
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1

pids=()

run_cv () {
  cv="$1"
  offs="$2"
  nick="$3"
  for i in $(seq 0 $((NSHARDS - 1))); do
    LOG=~/pj/out/commonlad_${nick}_w$(printf "%02d" "$i").log
    stdbuf -oL python3 -u phase8_criticalstate.py \
        --shard "$i" --nshards "$NSHARDS" --cvs "$cv" --seeds 1,2,3 \
        --offsets "$offs" --tag "commonlad_${nick}" > "$LOG" 2>&1 &
    pids+=($!)
    echo "commonlad cv=$cv worker $i -> $LOG (pid ${pids[-1]})"
  done
}

run_cv 0.0 "-0.0295333,-0.0045333,0.0204667,0.0454667,0.0704667,0.0954667" c00
run_cv 0.2 "-0.0430504,-0.0180504,0.0069496,0.0319496,0.0569496,0.0819496" c20
run_cv 0.5 "-0.0918540,-0.0668540,-0.0418540,-0.0168540,0.0081460,0.0331460" c50

fail=0
for p in "${pids[@]}"; do wait "$p" || fail=$((fail+1)); done
echo "=== commonlad DONE, $fail nonzero exits ==="
