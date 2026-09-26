#!/bin/bash
# PREREG_ROLLCHI Amendment 2: switch the sphere queue (G2) to a reverse-order clump pool.
# 1. kill the PARENT subshells of every rchiA / rchiC..F python (so no further sphere waves start);
# 2. kill rchiC..F pythons (later waves); rchiA pythons finish their in-flight shard as orphans;
# 3. start NW workers of achiR (arm-A CV-0 grid, reverse order) then achiBR (CV 0.5, reverse), both
#    skipping any key already written by any achi* pool.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/phase2_shape.py" > ~/pj/phase2_shape.py
md5sum "$SRC/phase2_shape.py" ~/pj/phase2_shape.py
for pid in $(pgrep -f "[p]hase2_shape.py --tag rchi[ACDEF] "); do
  pp=$(ps -o ppid= -p "$pid" | tr -d ' ')
  [ -n "$pp" ] && [ "$pp" != 1 ] && kill "$pp" 2>/dev/null
done
sleep 1
pkill -f "[p]hase2_shape.py --tag rchi[CDEF] "
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
NW=${1:-10}
AD="/home/cpb/pj/out/phase2_achi*_lognormal*.jsonl"
for i in $(seq 0 $((NW - 1))); do
  ( stdbuf -oL python3 -u phase2_shape.py --tag achiR --family lognormal --n 1700 \
        --alphas 1.75,1.5,1.3,1.2,1.1 --mus 0.3,0.0 --cvs 0.0 --seeds 3,2,1 --phats 3e-5,1e-3 \
        --dump-dir /home/cpb/pj/out/achi_dump --also-done "$AD" --shard "$i" --nshards "$NW" \
        >> ~/pj/out/achiR_w$(printf %02d "$i").log 2>&1
    stdbuf -oL python3 -u phase2_shape.py --tag achiBR --family lognormal --n 1700 \
        --alphas 1.75,1.5,1.3,1.2,1.1 --mus 0.3,0.0 --cvs 0.5 --seeds 3,2,1 --phats 3e-5,1e-3 \
        --dump-dir /home/cpb/pj/out/achi_dump --also-done "$AD" --shard "$i" --nshards "$NW" \
        >> ~/pj/out/achiBR_w$(printf %02d "$i").log 2>&1 ) &
done
echo "achiR pool: $NW workers $(date -Is)"
wait
echo "achiR pool done $(date -Is)"
bash "$SRC/scripts/sync_rollchi.sh"
