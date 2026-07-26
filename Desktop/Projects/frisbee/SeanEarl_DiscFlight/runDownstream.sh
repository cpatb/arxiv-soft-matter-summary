#!/bin/bash
# Stages 4d-6: sample the averaged CFD fields, validate, run the real Monte
# Carlo, compute throwability, make the figures, and check the hypotheses.
#
# Every step writes an artifact that the manuscript's numbers are generated
# from (P5_Synthesis/make_results_tex.py), so nothing is transcribed by hand.
#
# Usage: bash runDownstream.sh [n_jobs] [n_throws]
cd "$(dirname "$0")"
JOBS=${1:-20}
NTHROWS=${2:-150}
set -e

echo "############ 1. CFD validation (ABL, grid independence, convergence)"
cd P3_CFD
python3 validate_cfd.py \
    --empty cases/dir_090_abl \
    --fine cases/dir_090 \
    --coarse cases/dir_090_coarse \
    --out-json cfd_validation.json 2>&1 | tail -45

echo
echo "############ 2. sample the averaged fields to flow_fields.npz"
python3 sample_disc_altitude.py \
    cases/dir_090 cases/dir_210 cases/dir_270 \
    --out flow_fields.npz --dxy 1.0 2>&1 | tail -30

echo
echo "############ 3. Monte Carlo on the real field"
cd ../P4_Coupling
# Launch grid over the playing surface (OSM polygon x[-31,47], y[-10,15]).
python3 monte_carlo.py \
    --field ../P3_CFD/flow_fields.npz \
    --launch-grid -25 40 13 -8 12 10 \
    --headings 90 210 270 \
    --throws right_backhand right_forehand hammer \
    --throw-bearings 0 60 120 180 240 300 \
    --n-throws "$NTHROWS" \
    --jobs "$JOBS" \
    --out landing_stats.npz 2>&1 | tail -30

echo
echo "############ 4. throwability metric"
cd ../P5_Synthesis
python3 throwability.py \
    --fields ../P3_CFD/flow_fields.npz \
    --landings ../P4_Coupling/landing_stats.npz \
    --u-ref 5.0 --out Q.npz 2>&1 | tail -20

echo
echo "############ 5. hypothesis check"
python3 check_hypotheses.py \
    --Q Q.npz --fields ../P3_CFD/flow_fields.npz \
    --landings ../P4_Coupling/landing_stats.npz \
    --out-json hypothesis_check.json 2>&1 | tail -40

echo
echo "############ 6. figures"
mkdir -p figs
python3 plot_throwability_map.py Q.npz --outdir figs 2>&1 | tail -5
python3 paper_figures.py --Q Q.npz --fields ../P3_CFD/flow_fields.npz \
    --outdir figs 2>&1 | tail -5

echo
echo "############ 7. generated results macros"
python3 make_results_tex.py \
    --Q Q.npz --fields ../P3_CFD/flow_fields.npz \
    --landings ../P4_Coupling/landing_stats.npz \
    --validation ../P3_CFD/cfd_validation.json \
    --hypotheses hypothesis_check.json \
    --out results_generated.tex 2>&1 | tail -5

echo
echo "############ DOWNSTREAM COMPLETE ############"
