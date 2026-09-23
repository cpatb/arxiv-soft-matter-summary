"""grade_epssig.py -- PREREG_EPSSIG_2026-09-23.md. Writes out/grade_epssig.txt.
Input: out/epssig/<label>/ from `scripts/pull_jamwin.sh epssig 2026-09-23` (PowerShell) and the default-
EPS_MIN references out/kmerdrop/w0.25, w0.75. S2 as defined in grade_jamsig.sig. Written before the data."""
import json
import os
import re
import sys

import numpy as np

from grade_jamsig import sig

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
EPS = {"e003_w0.25": 0.003, "e003_w0.75": 0.003, "e03_w0.25": 0.03, "e03_w0.75": 0.03}
DEFAULT_GAP = None  # computed from the references below


def read(d, eps_want):
    st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
    a = np.loadtxt(os.path.join(d, "spread.dat"), comments="#", ndmin=2) \
        if os.path.exists(os.path.join(d, "spread.dat")) else np.zeros((0, 21))
    ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
    m = re.search(r"eps_min=([\d.e+-]+)", ban)
    eps_ok = bool(m) and abs(float(m.group(1)) - eps_want) <= 1e-9 * max(1.0, eps_want)
    healthy = bool(st["rc"] == 0 and a.size and a[-1, 1] >= 5.9 and eps_ok)
    return dict(healthy=healthy, eps_ok=eps_ok, rc=st["rc"], t=float(a[-1, 1]) if a.size else 0.0,
                J=float(a[:, 7].max()) if a.size else float("nan"),
                s2=sig(a, 2.0)[1] if a.size else float("nan"),
                rg=float(a[:, 11].max()) if a.size else float("nan"))


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("grade_epssig -- PREREG_EPSSIG_2026-09-23.md")
    ref = {w: read(os.path.join(OUT, "kmerdrop", w), 0.01) for w in ("w0.25", "w0.75")}
    for w, v in ref.items():
        say(f"  default EPS_MIN 0.01 {w}: healthy {v['healthy']} J {v['J']:.4f} S2 {v['s2']:+.4f} Rg {v['rg']:.4f}")
    gap0 = ref["w0.25"]["s2"] - ref["w0.75"]["s2"]
    say(f"  default S2 gap (w0.25 - w0.75) = {gap0:+.4f}; G2 bar = half = {gap0 / 2:.4f}")
    c = {}
    for lab, e in EPS.items():
        d = os.path.join(OUT, "epssig", lab)
        if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            say(f"  {lab}: not run yet")
            continue
        c[lab] = v = read(d, e)
        say(f"  {lab}: healthy {v['healthy']} (eps tag ok {v['eps_ok']}, rc {v['rc']}, t/tau {v['t']})  "
            f"J {v['J']:.4f}  S2 {v['s2']:+.4f}  Rg {v['rg']:.4f}")
    say("")
    a, b = c.get("e003_w0.25"), c.get("e003_w0.75")
    if a and b and a["healthy"] and b["healthy"]:
        g1 = b["s2"] < 0.1 <= a["s2"]
        say(f"G1 EPS_MIN 0.003: S2(w0.75) {b['s2']:.4f} < 0.1 <= S2(w0.25) {a['s2']:.4f}  {'PASS' if g1 else 'FAIL'}")
    else:
        g1 = None
        say("G1: INCOMPLETE (a cell is missing or dead; a dead cell's S2 is partial)")
    a, b = c.get("e03_w0.25"), c.get("e03_w0.75")
    if a and b and a["healthy"] and b["healthy"]:
        gap = a["s2"] - b["s2"]
        g2 = abs(gap) < abs(gap0) / 2
        say(f"G2 EPS_MIN 0.03: gap {gap:+.4f} vs half the default {abs(gap0) / 2:.4f}  {'PASS' if g2 else 'FAIL'}")
    else:
        g2 = None
        say("G2: INCOMPLETE")
    reading = {(True, True): "S2 is a signature of the jammed-vs-thickened VISCOSITY CONTRAST",
               (False, True): "S2 is not robust to the cap: WITHDRAWN as a candidate",
               (False, False): "S2 is not robust to the cap: WITHDRAWN as a candidate",
               (True, False): "S2 does not depend on the viscosity contrast: mechanism unexplained"}
    if g1 is not None and g2 is not None:
        say(f"READING (fixed in the prereg): {reading[(g1, g2)]}")
    elif g1 is False:
        say(f"READING (fixed in the prereg): {reading[(False, True)]}")
    # ---- Amendment 1 (registered before its cells): kmermap R1 under the cap -------------
    say("")
    say("Amendment 1: kmermap R1 (w1, x 0.88) under EPS_MIN 0.003 / 0.03; default J 0.0579 (U0 2), 0.1958 (U0 4, dead)")
    for e, tag in ((0.003, "e003"), (0.03, "e03")):
        js = []
        for u in ("2.0", "4.0"):
            d = os.path.join(OUT, "epssig", f"r1_{tag}_U{u}")
            if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
                say(f"  r1_{tag}_U{u}: not run yet")
                continue
            v = read(d, e)
            usable = v["eps_ok"] and (v["healthy"] or v["J"] > 0.01)
            say(f"  r1_{tag}_U{u}: healthy {v['healthy']} (eps tag {v['eps_ok']}, rc {v['rc']}, t/tau {v['t']})  "
                f"J {v['J']:.4f}{'' if v['healthy'] else ' (lower bound)'}")
            if usable:
                js.append(v["J"])
        say(f"H1 EPS_MIN {e}: J > 0.01 in >= 1 cell: "
            + ("PASS" if any(j > 0.01 for j in js) else "FAIL" if len(js) == 2 else "INCOMPLETE"))
    with open(os.path.join(OUT, "grade_epssig.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
