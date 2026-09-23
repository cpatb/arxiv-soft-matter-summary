"""grade_kmermapE03.py -- PREREG_KMERMAPE03_2026-09-23.md. Writes out/grade_kmermapE03.txt.
Input: out/kmermapE03/<label>/ from `scripts/pull_jamwin.sh kmermapE03 2026-09-23` (PowerShell); the
default-cap cells in out/kmermap_cells.json. Written before the data."""
import json
import math
import os
import re
import sys

import numpy as np

from grade_kmermap import PAR, PHI, Q, US, XS

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "kmermapE03")
EPS = 0.03


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("grade_kmermapE03 -- PREREG_KMERMAPE03_2026-09-23.md (EPS_MIN 0.03)")
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
                em = re.search(r"eps_min=([\d.e+-]+)", ban)
                prov = bool(m and jb and em) and all(math.isclose(float(v), t, rel_tol=1e-6) for v, t in zip(
                    m.groups(), (float(u), PHI[x], 0.62, pr["pm"], pr["ss"], pr["d"]))) \
                    and math.isclose(float(jb.group(1)), pr["jf"][x], rel_tol=1e-6) \
                    and math.isclose(float(em.group(1)), EPS, rel_tol=1e-9)
                p = os.path.join(d, "spread.dat")
                a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
                healthy = bool(st["rc"] == 0 and a.size > 0 and a[-1, 1] >= 5.9)
                c[(w, x, u)] = dict(prov=prov, healthy=healthy, rc=st["rc"],
                                   t=float(a[-1, 1]) if a.size else 0.0,
                                   J=float(a[:, 7].max()) if a.size else math.nan,
                                   beta=float(a[:, 2].max()) if a.size else math.nan,
                                   rg=float(a[:, 11].max()) if a.size else math.nan)
    say(f"cells: {len(c)} of 18; provenance (incl. eps_min 0.03) {sum(v['prov'] for v in c.values())}/{len(c)}")
    for k in sorted(c):
        v = c[k]
        say(f"  {k[0]}  {k[1]:.2f}  {k[2]}  healthy {str(v['healthy']):5s} J {v['J']:.4f}  beta {v['beta']:.4f}"
            f"  Rg {v['rg']:.4f}" + ("" if v["healthy"] and v["prov"] else f"   [EXCLUDED rc {v['rc']} t {v['t']}]"))
    ok = {k: v for k, v in c.items() if v["healthy"] and v["prov"]}
    say("")
    r1 = [ok[("w1", 0.88, u)]["J"] for u in US if ("w1", 0.88, u) in ok]
    say(f"R1' x 0.88 doublets J {[round(j, 4) for j in r1]}  "
        + ("PASS" if any(j > 0.01 for j in r1) else "FAIL" if len(r1) == 3 else "INCOMPLETE"))
    r2 = [(u, ok[("w1", 0.92, u)]["J"], ok[("w0", 0.92, u)]["J"]) for u in US
          if ("w1", 0.92, u) in ok and ("w0", 0.92, u) in ok]
    say(f"R2' x 0.92 (U0, J d, J s) {[(u, round(a, 4), round(b, 4)) for u, a, b in r2]}  "
        + ("PASS" if len(r2) == 3 and all(a > b for _, a, b in r2) else
           "FAIL" if r2 and not all(a > b for _, a, b in r2) else "INCOMPLETE"))
    r3 = []
    for w, x in (("w0", 0.92), ("w1", 0.88), ("w1", 0.92)):
        js = [ok[(w, x, u)]["J"] for u in US if (w, x, u) in ok]
        if len(js) == 3:
            r3.append(js[0] <= js[1] <= js[2])
            say(f"  R3' {w} x {x}: J over U0 1/2/4 = {[round(j, 4) for j in js]}")
    say(f"R3': {sum(r3)}/{len(r3)} non-decreasing" + ("  PASS" if len(r3) == 3 and all(r3) else
                                                    "  FAIL" if r3 and not all(r3) else "  INCOMPLETE"))
    r4 = [ok[("w1", x, u)]["beta"] <= ok[("w0", x, u)]["beta"] + Q
          for x in XS for u in US if ("w1", x, u) in ok and ("w0", x, u) in ok]
    say(f"R4' {sum(r4)}/{len(r4)}" + ("  PASS" if r4 and all(r4) else "  FAIL" if r4 else "  INCOMPLETE"))
    say(f"R5' health: {len(ok)}/{len(c)}")
    say("")
    ref = {(r["w"], r["x"], str(r["U0"])): r for r in
           json.load(open(os.path.join(HERE, "out", "kmermap_cells.json"), encoding="utf-8"))}

    def cls(j):
        return "J+" if j > 0.01 else "J0" if j == 0.0 else "amb"

    x1 = []
    for k, v in sorted(ok.items()):
        r = ref.get(k)
        if r and r["healthy"] and r["prov"]:
            same = cls(v["J"]) == cls(r["J"])
            x1.append(same)
            say(f"  X1 {k}: default {cls(r['J'])} ({r['J']:.4f})  e03 {cls(v['J'])} ({v['J']:.4f})"
                + ("" if same else "   CLASS CHANGED"))
    say(f"X1 cross-cap class: {sum(x1)}/{len(x1)} unchanged" + ("  PASS" if x1 and all(x1) else
                                                              "  FAIL" if x1 else "  INCOMPLETE"))
    # ---- Amendment 1 (registered before its cells): R1 at LEVEL 10 --------------------------
    say("")
    say("Amendment 1 (R1 at LEVEL 10, EPS_MIN 0.03):")
    r1l = []
    for u in ("2.0", "4.0"):   # the two registered cells (not US: U0 1 was not registered)
        d = os.path.join(HERE, "out", "kmermapE03L10", f"w1_x0.88_U{u}")
        if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            say(f"  L10 w1 x 0.88 U0 {u}: not run yet")
            continue
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        em = re.search(r"eps_min=([\d.e+-]+)", ban)
        prov = "_L10_" in ban and bool(em) and math.isclose(float(em.group(1)), EPS, rel_tol=1e-9)
        p = os.path.join(d, "spread.dat")
        a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
        healthy = bool(st["rc"] == 0 and a.size and a[-1, 1] >= 5.9 and prov)
        J = float(a[:, 7].max()) if a.size else math.nan
        l9 = c.get(("w1", 0.88, u), {}).get("J", math.nan)
        say(f"  L10 w1 x 0.88 U0 {u}: healthy {healthy} (tag/eps ok {prov}, rc {st['rc']}, t/tau "
            f"{a[-1, 1] if a.size else 0})  J {J:.4f}  (L9 {l9:.4f}, ratio {J / l9 if l9 else math.nan:.3f})")
        if healthy:
            r1l.append(J)
    say("R1L10: " + ("PASS" if any(j > 0.01 for j in r1l) else "FAIL" if len(r1l) == 2 else "INCOMPLETE")
        + f"  healthy J {[round(j, 4) for j in r1l]}")
    # ---- Amendment 2 (registered before its cells): the p_hat 1e-3 rung, doublets at x 0.91 --
    say("")
    say("Amendment 2 (chi from the p_hat 1e-3 rung; doublets at x 0.91, f_c 0.8806):")
    rg = []
    for lab, e in (("rg_d_x0.91_U2.0_e03", 0.03), ("rg_d_x0.91_U4.0_e03", 0.03),
                   ("rg_d_x0.91_U2.0_e01", 0.01), ("rg_d_x0.91_U4.0_e01", 0.01)):
        d = os.path.join(HERE, "out", "kmermapRG", lab)
        if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            say(f"  {lab}: not run yet")
            continue
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        em = re.search(r"eps_min=([\d.e+-]+)", ban)
        pm = re.search(r"phiM=([\d.]+)", ban)
        prov = bool(em and pm) and math.isclose(float(em.group(1)), e, rel_tol=1e-9)             and math.isclose(float(pm.group(1)), 0.55664, rel_tol=1e-6)
        p = os.path.join(d, "spread.dat")
        a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
        healthy = bool(st["rc"] == 0 and a.size and a[-1, 1] >= 5.9 and prov)
        J = float(a[:, 7].max()) if a.size else math.nan
        say(f"  {lab}: healthy {healthy} (prov {prov}, rc {st['rc']}, t/tau {a[-1, 1] if a.size else 0})  "
            f"J {J:.4f}{'' if healthy else ' (lower bound if dead)'}")
        if e == 0.03 and healthy:
            rg.append(J)
    say("RG1: " + ("PASS" if any(j > 0.01 for j in rg) else "FAIL" if len(rg) == 2 else "INCOMPLETE")
        + f"  healthy e03 J {[round(j, 4) for j in rg]}")
    with open(os.path.join(HERE, "out", "kmermapE03_cells.json"), "w", encoding="utf-8") as fh:
        json.dump([dict(w=k[0], x=k[1], U0=float(k[2]), **v) for k, v in sorted(c.items())], fh, indent=1)
    with open(os.path.join(HERE, "out", "grade_kmermapE03.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
