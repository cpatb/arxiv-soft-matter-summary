#!/bin/bash
# PREREG_ROLLCHI Amendment 2 (execution only): replace the sphere queue G2 by two pools.
#   X (5 workers): the mu 0.3 rolling ladder with the decisive etas first (1.0, 0.3, 0.1, then 0.03, 0.01),
#                  then the remaining sphere waves;
#   Y (5 workers): the arm-A CV-0 clump grid in REVERSE order, then CV 0.5 reverse.
# Every pool skips keys already written by any pool of the same arm (--also-done). In-flight G2 cells are lost.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/phase2_shape.py" > ~/pj/phase2_shape.py
md5sum "$SRC/phase2_shape.py" ~/pj/phase2_shape.py
for pid in $(pgrep -f "[p]hase2_shape.py --tag rchi[ACDEF] "); do
  pp=$(ps -o ppid= -p "$pid" | tr -d ' ')
  [ -n "$pp" ] && [ "$pp" != 1 ] && kill "$pp" 2>/dev/null
done
sleep 1
pkill -f "[p]hase2_shape.py --tag rchi[ACDEF] "
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
RD="/home/cpb/pj/out/phase2_rchi*_lognormal*.jsonl"
AD="/home/cpb/pj/out/phase2_achi*_lognormal*.jsonl"
S="--family lognormal --n 1700 --seeds 1,2,3 --phats 1e-3,3e-5"
x () {  # tag i ns args
  local tag=$1 i=$2 ns=$3; shift 3
  stdbuf -oL python3 -u phase2_shape.py --tag "$tag" $S "$@" --also-done "$RD" --shard "$i" --nshards "$ns" \
      >> ~/pj/out/${tag}_w$(printf %02d "$i").log 2>&1
}
for i in 0 1 2 3 4; do
  ( x rchiX $i 5 --alphas 1.0 --mus 0.3 --cvs 0.0,0.5 --eta-rolls 1.0,0.3,0.1
    x rchiX2 $i 5 --alphas 1.0 --mus 0.3 --cvs 0.0,0.5 --eta-rolls 0.03,0.01,0
    x rchiC $i 5 --alphas 1.0 --mus 0.1,1.0 --cvs 0.0 --eta-rolls 0,1.0,0.1,0.3,0.03,0.01
    x rchiD $i 5 --alphas 1.0 --mus 3.0 --cvs 0.5 --eta-rolls 0,1.0,0.1,0.3,0.03,0.01
    x rchiE $i 5 --alphas 1.0 --mus 0.1,0.3,1.0,3.0 --cvs 0.2 --eta-rolls 0,1.0,0.1,0.3,0.03,0.01
    x rchiF $i 5 --alphas 1.0 --mus 0.1,1.0 --cvs 0.5 --eta-rolls 0,1.0,0.1,0.3,0.03,0.01 ) &
  ( stdbuf -oL python3 -u phase2_shape.py --tag achiR --family lognormal --n 1700 \
        --alphas 1.75,1.5,1.3,1.2,1.1 --mus 0.3,0.0 --cvs 0.0 --seeds 3,2,1 --phats 3e-5,1e-3 \
        --dump-dir /home/cpb/pj/out/achi_dump --also-done "$AD" --shard "$i" --nshards 5 \
        >> ~/pj/out/achiR_w$(printf %02d "$i").log 2>&1
    stdbuf -oL python3 -u phase2_shape.py --tag achiBR --family lognormal --n 1700 \
        --alphas 1.75,1.5,1.3,1.2,1.1 --mus 0.3,0.0 --cvs 0.5 --seeds 3,2,1 --phats 3e-5,1e-3 \
        --dump-dir /home/cpb/pj/out/achi_dump --also-done "$AD" --shard "$i" --nshards 5 \
        >> ~/pj/out/achiBR_w$(printf %02d "$i").log 2>&1 ) &
done
echo "switched: pools X and Y, 10 workers $(date -Is)"
wait
echo "pools done $(date -Is)"
bash "$SRC/scripts/sync_rollchi.sh"
