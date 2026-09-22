#!/bin/bash
# Smoke for the straight-chain trimer (one cell). Syncs into a SEPARATE dir so the running
# kmer campaign's ~/pj copy is not touched.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
D=~/pjchain
mkdir -p $D/out $D/scripts
for f in packing.py phase2_shape.py phase2_phiJ.py shape.py psd.py measure.py \
         analyse_phase2.py mixture_model.py; do
  [ -f "$SRC/$f" ] && tr -d '\r' < "$SRC/$f" > $D/"$f"
done
cp ~/pj/*.py $D/ -n 2>/dev/null
for f in packing.py phase2_shape.py; do tr -d '\r' < "$SRC/$f" > $D/"$f"; done
cp ~/pj/scripts/env.sh $D/scripts/
cd $D || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1 PJ_MAX_ITER=1000000
python3 -u phase2_shape.py --tag cgate --alphas 1.0 --mus 0.3 --cvs 0.0 --seeds 1 \
    --phats 1e-3 --n 1700 --kmers "0:0.5:c" > out/cgate.log 2>&1
python3 - <<'EOF'
import json, glob
for f in glob.glob('out/phase2_cgate*.jsonl'):
    for l in open(f):
        r = json.loads(l)
        print(r['key'], r.get('failed'), r.get('error'), 'phi', r.get('phi'), 'z_nr', r.get('z_nr'),
              'mass', r.get('mass_over_vol_solid'), 'n3', r.get('n_triplets'),
              'chain', r.get('triplet_chain'), 'span', r.get('d_max_over_L'),
              'internal', r.get('n_contacts_internal'), 'conv', r.get('converged'))
EOF
tail -2 out/cgate.log
