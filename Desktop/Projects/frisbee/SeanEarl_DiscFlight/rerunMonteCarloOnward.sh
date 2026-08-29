#!/bin/bash
# Re-run runDownstream.sh steps 3-8 only: the Monte Carlo and everything derived
# from it. Steps 1-2 (CFD validation, field re-sampling) are deliberately skipped
# so ../P3_CFD/flow_fields.npz stays byte-identical to the field the current
# manuscript and GATE G1 were both computed on. The Monte Carlo is then the only
# variable that moves, which is the point: landing_stats.npz predates monte_carlo.py
# by eleven hours, carries none of the schema P4_Coupling/README.md documents, and
# does not record which gust model built it.
#
# Usage: bash rerunMonteCarloOnward.sh [n_jobs] [n_throws]
cd "$(dirname "$0")"
JOBS=${1:-18}
NTHROWS=${2:-150}
STAMP=$(date +%Y-%m-%d)
set -e

echo "############ 0. back up the artifacts this run will replace"
BK="_pre_rerun_${STAMP}"
mkdir -p "$BK/P4_Coupling" "$BK/P5_Synthesis"
for f in landing_stats.npz; do
    [ -f "P4_Coupling/$f" ] && cp -p "P4_Coupling/$f" "$BK/P4_Coupling/$f"
done
for f in Q.npz q_significance.json hypothesis_check.json \
         results_generated.tex significance_generated.tex; do
    [ -f "P5_Synthesis/$f" ] && cp -p "P5_Synthesis/$f" "$BK/P5_Synthesis/$f"
done
[ -d P5_Synthesis/figs ] && cp -rp P5_Synthesis/figs "$BK/P5_Synthesis/figs"
echo "backed up to $BK"
ls -l "$BK/P4_Coupling" "$BK/P5_Synthesis" | head -20

echo
echo "############ 3. Monte Carlo on the real field  (jobs=$JOBS throws=$NTHROWS)"
cd P4_Coupling
python monte_carlo.py \
    --field ../P3_CFD/flow_fields.npz \
    --launch-grid -25 40 13 -5 5 10 \
    --headings 90 210 270 \
    --throws right_backhand right_forehand hammer \
    --throw-bearings 0 60 120 180 240 300 \
    --n-throws "$NTHROWS" \
    --jobs "$JOBS" \
    --out landing_stats.npz 2>&1 | tail -30

echo
echo "############ 4. throwability metric"
cd ../P5_Synthesis
# --exclude-throws hammer is REQUIRED to reproduce the committed Q.npz; see the
# note in runDownstream.sh.
python throwability.py \
    --fields ../P3_CFD/flow_fields.npz \
    --landings ../P4_Coupling/landing_stats.npz \
    --exclude-throws hammer \
    --u-ref 5.0 --out Q.npz 2>&1 | tail -20

echo
echo "############ 5. is the direction ordering resolved?"
python q_significance.py \
    --Q Q.npz --landing ../P4_Coupling/landing_stats.npz \
    --tex significance_generated.tex --json q_significance.json 2>&1 | tail -25

echo
echo "############ 6. hypothesis check"
python check_hypotheses.py \
    --Q Q.npz --fields ../P3_CFD/flow_fields.npz \
    --landings ../P4_Coupling/landing_stats.npz \
    --significance q_significance.json \
    --out-json hypothesis_check.json 2>&1 | tail -40 || \
    echo "!! HYPOTHESIS CHECK FAILED -- continuing"

echo
echo "############ 7. figures"
mkdir -p figs
python plot_throwability_map.py Q.npz --outdir figs 2>&1 | tail -5
python paper_figures.py --Q Q.npz --fields ../P3_CFD/flow_fields.npz \
    --outdir figs 2>&1 | tail -5

echo
echo "############ 8. generated results macros"
python make_results_tex.py \
    --Q Q.npz --fields ../P3_CFD/flow_fields.npz \
    --landings ../P4_Coupling/landing_stats.npz \
    --validation ../P3_CFD/cfd_validation.json \
    --hypotheses hypothesis_check.json \
    --out results_generated.tex 2>&1 | tail -5

echo
echo "############ RERUN COMPLETE ############"
