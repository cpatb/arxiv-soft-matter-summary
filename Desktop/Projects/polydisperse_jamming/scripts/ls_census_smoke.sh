#!/bin/bash
# Smoke + identity check for the post-run contact census: same cell with and without --census must give identical phi.
set -u
SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
tr -d '\r' < "$SRC/ls_pack.py" > ~/pj/ls_pack.py
md5sum "$SRC/ls_pack.py" ~/pj/ls_pack.py
cd ~/pj || exit 1
source scripts/env.sh
export OMP_NUM_THREADS=1
for e in 0.5 1.0; do
python3 -u ls_pack.py --eps $e --mu 0.0 --phat 1e-3 --n 60 --seed 7 --tag csmoke_a --nodes 1000 --max-iter 4000 > /dev/null 2>&1 &
python3 -u ls_pack.py --eps $e --mu 0.0 --phat 1e-3 --n 60 --seed 7 --tag csmoke_b --nodes 1000 --max-iter 4000 --census > /dev/null 2>&1 &
done
wait
python3 - <<'PY'
import json, numpy as np, glob
a=[json.loads(l) for l in open('out/lsdem/ls_csmoke_a.jsonl')]
b=[json.loads(l) for l in open('out/lsdem/ls_csmoke_b.jsonl')]
for x in a:
    y=[r for r in b if r['eps']==x['eps']][0]
    print('eps',x['eps'],'phi',x['phi'],y['phi'],'identical',x['phi']==y['phi'] and x['z_all']==y['z_all'],'iters',x['iters'],'wall',round(x['wall_s']))
    d=np.load('out/lsdem/'+y['census'])
    print('  contacts',len(d['m1']),'2nc/n',2*len(d['m1'])/60,'z_all',y['z_all'],'m1 range',d['m1'].min().round(3),d['m1'].max().round(3),'arm max',max(d['arm1'].max(),d['arm2'].max()).round(4),'fn>0',(d['fn']>0).mean())
PY
