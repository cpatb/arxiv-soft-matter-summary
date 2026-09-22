"""grade_kmerdrop.py -- PREREG_KMERDROP_2026-09-22.md. Writes out/grade_kmerdrop.txt.
Input: out/kmerdrop/<label>/ from `scripts/pull_jamwin.sh kmerdrop`. Written before the data.
"""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "kmerdrop")
LADDER = (("w0.00", 0.00, 0.54932, 1.0), ("w0.25", 0.25, 0.53258, 1.0),
          ("w0.50", 0.50, 0.51708, 0.9036), ("w0.75", 0.75, 0.50282, 0.7937),
          ("w1.00", 1.00, 0.48856, 0.7075))


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("grade_kmerdrop -- PREREG_KMERDROP_2026-09-22.md (x = 0.85, sigma* 300, delta 1.3169)")
    rows = []
    for lab, w, pm, jf in LADDER:
        d = os.path.join(ROOT, lab)
        if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            say(f"  {lab}: not landed")
            continue
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        m = re.search(r"phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+) sigma\*=([\d.e+]+) Pa "
                      r"delta=([\d.]+)", ban)
        jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
        prov = bool(m and jb) and all(math.isclose(float(v), t, rel_tol=1e-6) for v, t in
                                      zip(m.groups(), (0.527, 0.62, pm, 300.0, 1.3169))) \
            and math.isclose(float(jb.group(1)), jf, rel_tol=1e-6)
        a = np.loadtxt(os.path.join(d, "spread.dat"), comments="#", ndmin=2) \
            if os.path.exists(os.path.join(d, "spread.dat")) else np.zeros((0, 21))
        healthy = st["rc"] == 0 and a.size and a[-1, 1] >= 5.9
        r = {"lab": lab, "w": w, "prov": prov, "healthy": bool(healthy), "rc": st["rc"],
             "t_end": st["last_t_over_tau"]}
        if a.size:
            r.update(beta=float(a[:, 2].max()), J=float(a[:, 7].max()),
                     fmean=float(a[:, 6].max()),
                     urms_late=float(a[a[:, 1] >= 5.0, 19].mean()) if (a[:, 1] >= 5.0).any()
                     else math.nan)
        rows.append(r)
        say(f"  {lab}: prov {prov} healthy {bool(healthy)} (rc {st['rc']}, t/tau {st['last_t_over_tau']})"
            + (f"  beta_max {r['beta']:.3f}  J(f_c={jf}) {r['J']:.4f}  <fhat>max {r['fmean']:.3f}"
               f"  u_rms(t/tau>=5) {r['urms_late']:.4f}" if a.size else ""))
    ok = [r for r in rows if r["prov"] and r["healthy"]]
    say(f"D3 health: {len(ok)}/{len(rows)} (of 5)")
    if len(ok) == 5:
        b = [r["beta"] for r in ok]
        d1 = all(b[i] > b[i + 1] for i in range(4))
        say(f"D1 beta_max monotone falling in w: {[round(v, 3) for v in b]}  {'PASS' if d1 else 'FAIL'}")
        j = [r["J"] for r in ok[2:]]
        d2 = ok[0]["J"] == 0 and ok[1]["J"] == 0 and j[0] > 0 and j[0] <= j[1] <= j[2]
        say(f"D2 J at w 0.5/0.75/1 = {[round(v, 4) for v in j]} (w 0, 0.25: "
            f"{ok[0]['J']}, {ok[1]['J']})  {'PASS' if d2 else 'FAIL'}")
    else:
        say("D1/D2: not graded -- fewer than 5 healthy cells; rows above are descriptive")
    with open(os.path.join(HERE, "out", "grade_kmerdrop.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
