#!/bin/bash
# Smoke for the per-contact dump (PREREG_ACTSHAPE_2026-09-22.md). NOT a measurement:
# N = 600, 4000 steps. Checks only that the Yade API returns (fn, ft, r1, r2) and that the
# npz lands with N_DUMP snapshots, for mu_p = 0 and 0.3.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out ~/pj/scripts
for f in shear.py phase9_constp.py packing.py measure.py psd.py; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
rm -rf ~/pj/out/actshape_smoke ~/pj/out/phase9_actsmoke.jsonl
python3 -u phase9_constp.py --cvs 0.2 --mus 0.0,0.3 --p-stars 3e-4 --phi-starts 0.55 \
    --gains 0.15 --seeds 1 --n 1000 --tag actsmoke --shear-steps 4000 --gamma-max 0.3 \
    --dump-dir ~/pj/out/actshape_smoke 2>&1 | tail -5
python3 - <<'EOF'
import glob, numpy as np, os
for f in sorted(glob.glob(os.path.expanduser('~/pj/out/actshape_smoke/*.npz'))):
    d = np.load(f)
    print(os.path.basename(f)[:60], 'snaps', len(d['n']), 'n', d['n'][:3],
          'fn>0', int((d['fn'] > 0).sum()), 'ft max', float(d['ft'].max()),
          'r range', float(d['r1'].min()), float(d['r1'].max()), 'mu', float(d['mu_p']))
EOF
