#!/bin/bash
# PREREG_ROLLCHI Amendment 3: alpha 1.9 / 1.95 at mu_p 0.3, CV 0 (12 cells), 3 workers. Foreground.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
cd ~/pj || exit 1
md5sum "$SRC/phase2_shape.py" ~/pj/phase2_shape.py "$SRC/packing.py" ~/pj/packing.py
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000
for i in 0 1 2; do
  stdbuf -oL python3 -u phase2_shape.py --tag achiN --family lognormal --n 1700 --alphas 1.9,1.95 --mus 0.3 \
      --cvs 0.0 --seeds 1,2,3 --phats 1e-3,3e-5 --dump-dir /home/cpb/pj/out/achi_dump \
      --shard "$i" --nshards 3 >> ~/pj/out/achiN_w0$i.log 2>&1 &
done
wait
echo "achiN done $(date -Is)"
bash "$SRC/scripts/sync_rollchi.sh"
