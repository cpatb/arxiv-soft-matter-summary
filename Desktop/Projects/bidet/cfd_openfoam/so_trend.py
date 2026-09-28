"""so_trend.py -- worst initial residual per iteration for the PREREG_so_triplet cases,
sampled every N iterations, plus a log-linear projection of when it crosses 1e-5.

Run inside WSL: python3 so_trend.py [every]. Read-only.
"""
import glob
import math
import os
import re
import sys

EVERY = int(sys.argv[1]) if len(sys.argv) > 1 else 200
TARGET = 1e-5
pat = re.compile(r"Solving for (\w+),\s+Initial residual = ([0-9.eE+-]+)")
for d in sorted(glob.glob(os.path.expanduser("~/bidet_so/case*"))):
    log = os.path.join(d, "log.simpleFoam")
    if not os.path.exists(log):
        continue
    it, cur, hist = None, {}, []
    with open(log, errors="replace") as fh:
        for line in fh:
            if line.startswith("Time = "):
                if it is not None and cur:
                    # Amendment 2: Uz is round-off on the wedge; H1 excludes it
                    nz = {k: v for k, v in cur.items() if k != "Uz"}
                    hist.append((it, max(nz.values()), max(nz, key=nz.get)))
                it, cur = int(float(line.split("=")[1])), {}
            else:
                m = pat.search(line)
                if m and m.group(1) not in cur:
                    cur[m.group(1)] = float(m.group(2))
    if not hist:
        continue
    print(os.path.basename(d))
    for i, r, f in hist[::EVERY]:
        print(f"   it {i:6d}  worst {r:.3e} ({f})")
    i1, r1, f1 = hist[-1]
    print(f"   it {i1:6d}  worst {r1:.3e} ({f1})  <- last")
    # project from the last 40% of the history
    k = max(0, int(0.6 * len(hist)))
    i0, r0, _ = hist[k]
    if i1 > i0 and r1 < r0:
        rate = (math.log(r1) - math.log(r0)) / (i1 - i0)
        need = (math.log(TARGET) - math.log(r1)) / rate if r1 > TARGET else 0
        print(f"   log-linear projection: {need:.0f} more iterations to {TARGET:g}")
    else:
        print("   not decreasing over the last 40% -- no projection")
