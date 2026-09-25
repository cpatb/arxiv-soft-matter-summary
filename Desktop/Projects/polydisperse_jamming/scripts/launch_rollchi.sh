#!/bin/bash
# rollchi -- PREREG_ROLLCHI_2026-09-25.md. Arm R (468 rolling-resistant sphere cells) + arm A (150 clump
# cells) + one G0 wiring cell. Foreground; hold wsl.exe open from PowerShell. Restartable (keys skipped).
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out ~/pj/out/achi_dump ~/pj/scripts
for f in packing.py phase2_shape.py phase2_phiJ.py shape.py psd.py measure.py; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
md5sum "$SRC/packing.py" ~/pj/packing.py "$SRC/phase2_shape.py" ~/pj/phase2_shape.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
COMMON="--family lognormal --seeds 1,2,3 --phats 1e-3,3e-5 --n 1700"
pids=()
go () {  # tag nshards args...
  local tag=$1 ns=$2; shift 2
  for i in $(seq 0 $((ns - 1))); do
    stdbuf -oL python3 -u phase2_shape.py --tag "$tag" $COMMON "$@" --shard "$i" --nshards "$ns" \
        > ~/pj/out/${tag}_w$(printf %02d "$i").log 2>&1 &
    pids+=($!)
  done
}
go rchi 12 --alphas 1.0 --mus 0.1,0.3,1.0,3.0 --cvs 0.0,0.2,0.5 --eta-rolls 0,0.01,0.03,0.1,0.3,1.0
go rchi0 1 --alphas 1.0 --mus 0.0 --cvs 0.0,0.2,0.5 --eta-rolls 0,1.0
go achi 4 --alphas 1.1,1.2,1.3,1.5,1.75 --mus 0.0,0.3,3.0 --cvs 0.0 --dump-dir /home/cpb/pj/out/achi_dump
go achi5 3 --alphas 1.1,1.2,1.3,1.5,1.75 --mus 0.0,0.3 --cvs 0.5 --dump-dir /home/cpb/pj/out/achi_dump
stdbuf -oL python3 -u phase2_shape.py --tag achig --family lognormal --seeds 1 --phats 1e-3 --n 1700 \
    --alphas 2.0 --mus 0.3 --cvs 0.0 > ~/pj/out/achig_w00.log 2>&1 &
pids+=($!)
echo "rollchi: ${#pids[@]} workers"
date -u +"start %Y-%m-%dT%H:%M:%SZ"
for p in "${pids[@]}"; do wait "$p"; done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
bash "$SRC/scripts/sync_rollchi.sh"
