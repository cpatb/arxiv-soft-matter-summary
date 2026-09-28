"""binload.py -- POST HOC (added 2026-09-28 05:20, after the Am. 7 first-crossing ratio was seen): the load-size exponent
in the BINARY arm, measured on an observable independent of activation.

    python lfdem/binload.py   ->  out/binload_lfsig2bi.txt

Per cell (out/lfsig2bi/), interaction snapshots in strain [2, 4]: mean total normal load (lubrication normal col 8 +
contact normal col 12 + repulsion col 16) on COMPRESSIVE pairs that are NOT YET TOUCHING (state 0), for SS and LL pair
types. With R*_LL/R*_SS = 2, beta = ln(F_LL/F_SS)/ln 2 (mean field: 2). Under a Derjaguin threshold (~R*) a pair
activates when sigma R*^(beta-1) exceeds a constant, so the predicted activation exponent is m = beta - 1, to be compared
with the first-crossing m_eff of `pairtype.py` (a size-independent critical load would predict m = beta instead).
"""
from __future__ import annotations

import glob
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mech_lfsig import snapshots  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "out")


def main() -> int:
    lines = ["binload -- POST HOC: LL/SS load on not-yet-touching compressive pairs in the binary arm; beta = ln ratio / ln 2"]
    rows = []
    for d in sorted(glob.glob(os.path.join(OUT, "lfsig2bi", "bi_*"))):
        kv = dict(re.findall(r"(\w+)=(\S+)", open(os.path.join(d, "META"), encoding="utf-8").read()))
        r = np.loadtxt(os.path.join(d, "conf.dat"), comments="#", ndmin=2)[:, 3]
        small = r < 1.5
        FS, FL = [], []
        for strain, _q, s in snapshots(glob.glob(os.path.join(d, "int_*.dat"))[0]):
            if not (2.0 <= strain <= 4.0):
                continue
            i, j = s[:, 0].astype(int), s[:, 1].astype(int)
            ns = small[i].astype(int) + small[j].astype(int)
            f = s[:, 7] + s[:, 11] + s[:, 15]
            nt = (s[:, 2] < 1) & (f > 0)
            FS.extend(f[nt & (ns == 2)]); FL.extend(f[nt & (ns == 0)])
        ratio = float(np.mean(FL) / np.mean(FS))
        rows.append((float(kv["SIG"]), int(kv["SEED"]), ratio, math.log(ratio) / math.log(2)))
    for sig in sorted({x[0] for x in rows}):
        g = [x for x in rows if x[0] == sig]
        lines.append(f"  sigma {sig:5.1f}: LL/SS " + " ".join(f"{x[2]:.2f}" for x in g) + "   beta " + " ".join(f"{x[3]:.2f}" for x in g))
    b = np.array([x[3] for x in rows if x[0] <= 4.5])
    lines.append(f"  sigma <= 4.5 (the activation range): beta {b.mean():.3f} (sd {b.std(ddof=1):.3f}, n {b.size}) "
                 f"-> Derjaguin-predicted m = {b.mean() - 1:.3f}; constant-force-predicted m = {b.mean():.3f}")
    txt = "\n".join(lines) + "\n"
    print(txt)
    with open(os.path.join(OUT, "binload_lfsig2bi.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
