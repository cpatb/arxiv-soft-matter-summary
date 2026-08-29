#!/usr/bin/env bash
# Post-processing. Run after run_pipeline.py has produced results/*_meta.json
# for every clip. The five stages are independent and are run concurrently;
# thread pools are capped so they do not oversubscribe the cores.
set -u
cd "$(dirname "$0")"

export OMP_NUM_THREADS=3 OPENBLAS_NUM_THREADS=3 MKL_NUM_THREADS=3

mkdir -p logs

echo "=== post-processing started $(date) ==="

# Body length must be re-measured before analysis: everything is divided by it,
# and the batch pass used a scale-dependent acceptance window that biased it low.
python -W ignore src/measure_bodylength.py > logs/bodylength.log 2>&1
echo "measure_bodylength: exit $?"

python -W ignore src/analyse_results.py > logs/analyse.log 2>&1
echo "analyse_results: exit $?"

python -W ignore src/camera_qa_figure.py    > logs/camera_qa.log     2>&1 &
python -W ignore src/collective_state.py    > logs/collective.log    2>&1 &
python -W ignore src/correlation_null.py    > logs/corr_null.log     2>&1 &
python -W ignore src/fringe_analysis.py     > logs/fringe.log        2>&1 &
python -W ignore src/make_visual_figures.py > logs/visual.log        2>&1 &
python -W ignore src/fission_analysis.py    > logs/fission.log       2>&1 &
wait

python -W ignore src/digest.py > logs/digest.log 2>&1
echo "=== post-processing finished $(date) ==="
# digest is in this list deliberately: it was omitted once, and its log ended in
# a traceback for a day without anything reporting it.
for f in bodylength analyse camera_qa collective corr_null fringe visual fission digest; do
  printf '%-12s : ' "$f"
  if grep -qE "Traceback|Error" "logs/$f.log" 2>/dev/null; then
    echo "ERRORS -- see logs/$f.log"
  else
    echo "ok"
  fi
done
ls -la figures/*.png 2>/dev/null | wc -l | xargs echo "figures produced:"
