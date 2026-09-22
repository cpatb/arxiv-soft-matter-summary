#!/bin/bash
# PREREG_KMER gate K0 (identity: the historical path through the new code) and K1 (smoke:
# one w2 = 0.5 and one w3 = 0.5 cell). Foreground; hold wsl.exe open from PowerShell.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out/shape ~/pj/out/kmer_dump ~/pj/scripts
for f in packing.py phase2_shape.py shape.py psd.py measure.py; do
  tr -d '\r' < "$SRC/$f" > ~/pj/"$f" || exit 1
done
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
export PJ_MAX_ITER=1000000

( python3 -u phase2_shape.py --tag kgate --alphas 1.0 --mus 0.3 --cvs 0.0 --seeds 1 \
    --phats 1e-3 --n 1700 > ~/pj/out/kgate_a1.log 2>&1 ) &
( python3 -u phase2_shape.py --tag kgate --alphas 2.0 --mus 0.3 --cvs 0.0 --seeds 1 \
    --phats 1e-3 --n 1700 > ~/pj/out/kgate_a2.log 2>&1 ) &
( python3 -u phase2_shape.py --tag kgate2 --alphas 1.0 --mus 0.3 --cvs 0.0 --seeds 1 \
    --phats 1e-3 --n 1700 --kmers "0.5:0.0" --dump-dir ~/pj/out/kmer_dump \
    > ~/pj/out/kgate_w2.log 2>&1 ) &
( python3 -u phase2_shape.py --tag kgate3 --alphas 1.0 --mus 0.3 --cvs 0.0 --seeds 1 \
    --phats 1e-3 --n 1700 --kmers "0.0:0.5" --dump-dir ~/pj/out/kmer_dump \
    > ~/pj/out/kgate_w3.log 2>&1 ) &
wait
python3 - <<'EOF'
import json, glob, os
ref = {}
for f in glob.glob(os.path.expanduser('~/pj/out/shape/phase2_shapechi*.jsonl')) + \
         glob.glob('/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/out/shape/phase2_shapechi*.jsonl'):
    for l in open(f):
        r = json.loads(l); ref[r['key']] = r
for f in sorted(glob.glob(os.path.expanduser('~/pj/out/phase2_kgate*.jsonl'))):
    for l in open(f):
        r = json.loads(l)
        if r.get('failed'):
            print('FAILED', r['key'], r.get('error')); continue
        base = ref.get(r['key'])
        if base is not None:
            same = all(r[k] == base[k] for k in ('phi', 'iters', 'n_contacts', 'z_nr', 'force'))
            print(f"K0 {r['key']}: phi {r['phi']!r} vs {base['phi']!r}  identical={same}")
        else:
            print(f"K1 {r['key']}: phi {r['phi']:.5f} z_nr {r['z_nr']:.3f} mass/vol "
                  f"{r['mass_over_vol_solid']:.6f} n2 {r['n_doublets']} n3 {r['n_triplets']} "
                  f"conv {r['converged']} internal {r['n_contacts_internal']} notes {r['notes']}")
EOF
ls -la ~/pj/out/kmer_dump
cp ~/pj/out/phase2_kgate*.jsonl "$SRC/out/shape/"
