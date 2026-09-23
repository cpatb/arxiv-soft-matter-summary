"""grade_kmermap.py -- PREREG_KMERMAP_2026-09-23.md. Writes out/grade_kmermap.txt.
Input: out/kmermap/<label>/ from `scripts/pull_jamwin.sh kmermap` (PowerShell). Written before data."""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "kmermap")
XS = (0.85, 0.88, 0.92)
US = ("1.0", "2.0", "4.0")
PHI = {0.85: 0.527, 0.88: 0.5456, 0.92: 0.5704}
PAR = {"w0": dict(pm=0.56315, ss=300.0, d=1.3169, jf={0.85: 1.0, 0.88: 1.0, 0.92: 0.8724}),
       "w1": dict(pm=0.53512, ss=235.8, d=1.1163, jf={0.85: 1.0, 0.88: 0.8766, 0.92: 0.5844})}
Q = 1 / 64


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("grade_kmermap -- PREREG_KMERMAP_2026-09-23.md")
    c = {}
    for w, pr in PAR.items():
        for x in XS:
            for u in US:
                lab = f"{w}_x{x:.2f}_U{u}"
                d = os.path.join(ROOT, lab)
                if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
                    continue
                st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
                ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
                m = re.search(r"U0=([\d.]+) .*phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+) sigma\*=([\d.e+]+) "
                              r"Pa delta=([\d.]+)", ban)
                jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
                prov = bool(m and jb) and all(math.isclose(float(v), t, rel_tol=1e-6) for v, t in zip(
                    m.groups(), (float(u), PHI[x], 0.62, pr["pm"], pr["ss"], pr["d"]))) and \
                    math.isclose(float(jb.group(1)), pr["jf"][x], rel_tol=1e-6)
                a = np.loadtxt(os.path.join(d, "spread.dat"), comments="#", ndmin=2) \
                    if os.path.exists(os.path.join(d, "spread.dat")) else np.zeros((0, 21))
                healthy = bool(st["rc"] == 0 and a.size > 0 and a[-1, 1] >= 5.9)
                c[(w, x, u)] = dict(prov=prov, healthy=healthy, rc=st["rc"],
                                   t=st["last_t_over_tau"],
                                   J=float(a[:, 7].max()) if a.size else math.nan,
                                   beta=float(a[:, 2].max()) if a.size else math.nan,
                                   rg=float(a[:, 11].max()) if a.size else math.nan,
                                   fm=float(a[:, 6].max()) if a.size else math.nan)
    say(f"cells: {len(c)} of 18; provenance {sum(v['prov'] for v in c.values())}/{len(c)}")
    say("  w   x     U0   healthy  J(f_c)   beta_max  Rg_max  <fhat>max")
    for k in sorted(c):
        v = c[k]
        say(f"  {k[0]}  {k[1]:.2f}  {k[2]}  {str(v['healthy']):5s}   {v['J']:.4f}   {v['beta']:.4f}"
            f"   {v['rg']:.4f}  {v['fm']:.3f}" + ("" if v["healthy"] and v["prov"] else
                                               f"   [EXCLUDED rc {v['rc']} t {v['t']}]"))
    ok = {k: v for k, v in c.items() if v["healthy"] and v["prov"]}
    say("")
    r1 = [ok[("w1", 0.88, u)]["J"] for u in US if ("w1", 0.88, u) in ok]
    say(f"R1 x 0.88: doublets J = {[round(j, 4) for j in r1]}  "
        + ("PASS" if any(j > 0.01 for j in r1) else "FAIL" if len(r1) == 3 else "INCOMPLETE"))
    r2 = [(u, ok[("w1", 0.92, u)]["J"], ok[("w0", 0.92, u)]["J"]) for u in US
          if ("w1", 0.92, u) in ok and ("w0", 0.92, u) in ok]
    say(f"R2 x 0.92: (U0, J doublet, J singlet) {[(u, round(a, 4), round(b, 4)) for u, a, b in r2]}  "
        + ("PASS" if len(r2) == 3 and all(a > b for _, a, b in r2) else
           "FAIL" if r2 and not all(a > b for _, a, b in r2) else "INCOMPLETE (deaths)"))
    r3 = []
    for w, x in (("w0", 0.92), ("w1", 0.88), ("w1", 0.92)):
        js = [ok[(w, x, u)]["J"] for u in US if (w, x, u) in ok]
        if len(js) == 3:
            r3.append(js[0] <= js[1] <= js[2])
            say(f"  R3 {w} x {x}: J over U0 1/2/4 = {[round(j, 4) for j in js]}")
    say(f"R3: {sum(r3)}/{len(r3)} non-decreasing" + ("  PASS" if r3 and all(r3) else
                                                   "  FAIL" if r3 else "  INCOMPLETE"))
    r4 = [ok[("w1", x, u)]["beta"] <= ok[("w0", x, u)]["beta"] + Q
          for x in XS for u in US if ("w1", x, u) in ok and ("w0", x, u) in ok]
    say(f"R4 beta(doublets) <= beta(singlets) + one quantum: {sum(r4)}/{len(r4)}"
        + ("  PASS" if r4 and all(r4) else "  FAIL" if r4 else "  INCOMPLETE"))
    say(f"R5 health: {len(ok)}/{len(c)}")

    # ---- POST-HOC (added after the data; does NOT change any verdict above) ------------------
    # A dead cell's J_max is a max over a TRUNCATED record, so it is a LOWER BOUND on the
    # full-run max (a running max can only rise). A healthy cell's J_max is exact. So an
    # ordering "A > B" is decided when lower(A) > exact-or-upper(B); B dead has no upper bound.
    say("")
    say("POST-HOC censored reading (dead cells as lower bounds; registered verdicts unchanged):")
    for u in US:
        a, b = c.get(("w1", 0.92, u)), c.get(("w0", 0.92, u))
        if not (a and b and a["prov"] and b["prov"]):
            continue
        a_ok, b_ok = a["healthy"], b["healthy"]
        if b_ok and a["J"] > b["J"]:
            v = "doublet > singlet DECIDED" + ("" if a_ok else " (doublet J is a lower bound)")
        elif a_ok and b_ok:
            v = "doublet <= singlet"
        else:
            v = "UNDECIDABLE (singlet dead: no upper bound)"
        say(f"  R2' x 0.92 U0 {u}: J doublet {a['J']:.4f}{'' if a_ok else ' (>=)'}"
            f" vs singlet {b['J']:.4f}{'' if b_ok else ' (>=)'}: {v}")
    for w, x in (("w0", 0.92), ("w1", 0.88), ("w1", 0.92)):
        cs = [c.get((w, x, u)) for u in US]
        if all(v and v["prov"] for v in cs):
            say(f"  R3' {w} x {x}: J over U0 1/2/4 = "
                + ", ".join(f"{v['J']:.4f}{'' if v['healthy'] else '(>=)'}" for v in cs))
    with open(os.path.join(HERE, "out", "kmermap_cells.json"), "w", encoding="utf-8") as fh:
        json.dump([dict(w=k[0], x=k[1], U0=float(k[2]), **v) for k, v in sorted(c.items())], fh,
                  indent=1)
    with open(os.path.join(HERE, "out", "grade_kmermap.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
