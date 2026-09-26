#!/bin/bash
# PREREG_ROLLCHI execution order: the mu-0 control queue (rchi0A/B, 2 workers) is restarted so the eta = 1.0 partners
# of already-finished eta = 0 cells run first (R0's inert half needs PAIRS). Grid unchanged; finished keys skipped.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
for pid in $(pgrep -f "[p]hase2_shape.py --tag rchi0[AB] "); do
  pp=$(ps -o ppid= -p "$pid" | tr -d ' ')
  [ -n "$pp" ] && [ "$pp" != 1 ] && kill "$pp" 2>/dev/null
done
sleep 1
pkill -f "[p]hase2_shape.py --tag rchi0[AB] "
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
RD="/home/cpb/pj/out/phase2_rchi*_lognormal*.jsonl"
for i in 0 1; do
  ( stdbuf -oL python3 -u phase2_shape.py --tag rchi0C --family lognormal --n 1700 --seeds 1,2,3 --phats 1e-3,3e-5 \
        --alphas 1.0 --mus 0.0 --cvs 0.0 --eta-rolls 1.0,0 --also-done "$RD" --shard $i --nshards 2 >> ~/pj/out/rchi0C_w0$i.log 2>&1
    stdbuf -oL python3 -u phase2_shape.py --tag rchi0D --family lognormal --n 1700 --seeds 1,2,3 --phats 1e-3,3e-5 \
        --alphas 1.0 --mus 0.0 --cvs 0.5,0.2 --eta-rolls 0,1.0 --also-done "$RD" --shard $i --nshards 2 >> ~/pj/out/rchi0D_w0$i.log 2>&1 ) &
done
echo "rchi0 reordered $(date -Is)"
wait
bash "$SRC/scripts/sync_rollchi.sh"
