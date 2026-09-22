"""W2'' for PREREG_ACTSHAPE: recompute force_moments from each dump's LAST snapshot ON THE
PLATFORM THAT WROTE THE RECORD (WSL), and require bit-equality with the record's `force`.

Why: the registered W2 compares on Windows, where np.polyfit (LAPACK) differs from WSL's at
~1e-15 in `tail_slope` while n/mean/cv are bit-identical -- a platform difference, not a
network one. Run in WSL from ~/pj:  python3 scripts/w2_same_platform.py > out/actshape_w2wsl.txt
"""
import glob
import json
import os

import numpy as np

from measure import force_moments

n_ok = n_all = 0
for f in sorted(glob.glob("out/phase9_actshape_*.jsonl")):
    for line in open(f):
        r = json.loads(line)
        p = os.path.join("out", "actshape_dump", r["key"] + ".npz")
        if "error" in r or not os.path.exists(p):
            continue
        d = np.load(p)
        e = np.concatenate([[0], np.cumsum(d["n"])])
        ok = force_moments(d["fn"][e[-2]:e[-1]]) == r["force"]
        n_all += 1
        n_ok += ok
        print(("OK   " if ok else "DIFF ") + r["key"])
print(f"W2'' same-platform identical: {n_ok}/{n_all}")
