#!/bin/bash
# rollchi -- PREREG_ROLLCHI_2026-09-25.md, Amendment 1 (EXECUTION ORDER ONLY; the registered grid is unchanged).
# 21 workers in 4 queues; each worker runs its shard of each wave in turn. Waves carry their own tag so that
# resume files do not collide; the grader dedups by cell key. Foreground; hold wsl.exe open from PowerShell.
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
C="--family lognormal --seeds 1,2,3 --phats 1e-3,3e-5 --n 1700"
D="--dump-dir /home/cpb/pj/out/achi_dump"
ETA="--eta-rolls 0,0.01,0.03,0.1,0.3,1.0"
run () {  # tag shard nshards args...
  local tag=$1 i=$2 ns=$3; shift 3
  stdbuf -oL python3 -u phase2_shape.py --tag "$tag" $C "$@" --shard "$i" --nshards "$ns" \
      >> ~/pj/out/${tag}_w$(printf %02d "$i").log 2>&1
}
pids=()
for i in $(seq 0 6); do   # G1: clump ladder (7 workers)
  ( run achiA $i 7 --alphas 1.1,1.2,1.3,1.5,1.75 --mus 0.0,0.3 --cvs 0.0 $D
    run achiB $i 7 --alphas 1.1,1.2,1.3,1.5,1.75 --mus 0.0,0.3 --cvs 0.5 $D
    run achiC $i 7 --alphas 1.1,1.2,1.3,1.5,1.75 --mus 3.0 --cvs 0.0 $D ) &
  pids+=($!)
done
for i in $(seq 0 9); do   # G2: rolling spheres (10 workers)
  ( run rchiA $i 10 --alphas 1.0 --mus 0.3 --cvs 0.0,0.5 $ETA
    run rchiC $i 10 --alphas 1.0 --mus 0.1,1.0 --cvs 0.0 $ETA
    run rchiD $i 10 --alphas 1.0 --mus 3.0 --cvs 0.5 $ETA
    run rchiE $i 10 --alphas 1.0 --mus 0.1,0.3,1.0,3.0 --cvs 0.2 $ETA
    run rchiF $i 10 --alphas 1.0 --mus 0.1,1.0 --cvs 0.5 $ETA ) &
  pids+=($!)
done
for i in 0 1; do          # G3: mu 0 controls (2 workers)
  ( run rchi0A $i 2 --alphas 1.0 --mus 0.0 --cvs 0.0,0.5 --eta-rolls 0,1.0
    run rchi0B $i 2 --alphas 1.0 --mus 0.0 --cvs 0.2 --eta-rolls 0,1.0 ) &
  pids+=($!)
done
for i in 0 1; do          # G4: mu 3 rolling, CV 0 (2 workers)
  ( run rchiB $i 2 --alphas 1.0 --mus 3.0 --cvs 0.0 $ETA ) &
  pids+=($!)
done
echo "rollchi: ${#pids[@]} workers"
date -u +"start %Y-%m-%dT%H:%M:%SZ"
for p in "${pids[@]}"; do wait "$p"; done
date -u +"end %Y-%m-%dT%H:%M:%SZ"
bash "$SRC/scripts/sync_rollchi.sh"
