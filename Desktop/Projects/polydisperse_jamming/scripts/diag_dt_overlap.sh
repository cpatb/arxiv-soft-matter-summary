#!/bin/bash
# Is the deep inter-grain overlap in frictional clump packs a time-step artefact?
# One bad cell (mu 3, w2 = 1, p 1e-3, seed 1; 3.66 % of contacts beyond u/r 0.05 at the default
# DT_SAFETY 0.5) re-run at DT_SAFETY 0.5 / 0.25 / 0.1. Separate dir; the running campaigns are untouched.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
D=~/pjdt
mkdir -p $D/out $D/scripts
cp ~/pj/*.py $D/ 2>/dev/null
for f in packing.py phase2_shape.py; do tr -d '\r' < "$SRC/$f" > $D/"$f"; done
cp ~/pj/scripts/env.sh $D/scripts/
cd $D || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1 PJ_MAX_ITER=1000000
for dts in 0.5 0.25 0.1; do
  ( PJ_DT_SAFETY=$dts python3 -u phase2_shape.py --tag dt${dts} --alphas 1.0 --mus 3.0 --cvs 0.0 \
      --seeds 1 --phats 1e-3 --n 1700 --kmers "1.0:0" --dump-dir $D/out/dump_$dts \
      > out/dt_$dts.log 2>&1 ) &
done
wait
python3 - <<'EOF'
import glob, json, numpy as np, os
for dts in ("0.5", "0.25", "0.1"):
    for f in glob.glob(f"out/phase2_dt{dts}_*.jsonl"):
        for l in open(f):
            r = json.loads(l)
            d = np.load(glob.glob(f"out/dump_{dts}/*.npz")[0])
            u = d["fn"] / 0.25
            print(f"DT_SAFETY {dts}: phi {r['phi']:.5f} iters {r['iters']} conv {r['converged']} "
                  f"dt {r['dt']:.3e} frac u/r>0.05 {(u > 0.05).mean():.4f} max u/r {u.max():.3f} "
                  f"Z_nr {r['z_nr']:.3f}")
EOF
