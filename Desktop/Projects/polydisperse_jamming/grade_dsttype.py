"""grade_dsttype.py -- PREREG_DSTTYPE_2026-09-26.md. Writes out/grade_dsttype.txt.

Input: out/dsttype/<label>/ from `scripts/pull_jamwin.sh dsttype 2026-09-26` (run from PowerShell), and the registered
cell list out/dsttype_cells.txt (provenance checked per cell against the solver banner). Written before any cell
landed. Health graded SEPARATELY: a dead cell is excluded and printed, never a FAIL.
spread.dat columns: 1 t/tau, 4 h/D0, 11 Rg_r/D0, 19 u_r (as grade_cvmap / grade_jamsig).
"""
import json
import math
import os
import re

import numpy as np

from grade_jamsig import sig

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "dsttype")
CELLS = os.path.join(HERE, "out", "dsttype_cells.txt")
PHI0, N_EXP, EPS = 0.62, 2.0, 0.03
XDST = {"C5": 0.89854, "C0": 0.87919}          # registered (PREREG table), from thickening_type x_dst
CHI = {"C5": 0.56855 / PHI0, "C0": 0.56459 / PHI0}
L = []


def say(s=""):
    print(s)
    L.append(s)


def load():
    reg = {}
    for line in open(CELLS, encoding="utf-8"):
        f = line.split()
        if f:
            reg[f[0]] = dict(phi=float(f[1]), d=float(f[2]), ss=float(f[3]), jf=float(f[4]), pm=float(f[5]),
                             u0=float(f[6]), set=f[0][:2], x=float(f[0].split("_x")[1].split("_")[0]))
    c = {}
    for lab, r in reg.items():
        d = os.path.join(ROOT, lab)
        if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            c[lab] = dict(r, missing=True, healthy=False, prov=False)
            continue
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        m = re.search(r"U0=([\d.]+) .*phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+) sigma\*=([\d.e+]+) "
                      r"Pa delta=([\d.]+)", ban)
        jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
        em = re.search(r"eps_min=([\d.e+-]+)", ban)
        prov = bool(m and jb and em) and math.isclose(float(em.group(1)), EPS, rel_tol=1e-9) and all(
            math.isclose(float(v), t, rel_tol=1e-6)
            for v, t in zip(m.groups(), (r["u0"], r["phi"], PHI0, r["pm"], r["ss"], r["d"]))) \
            and math.isclose(float(jb.group(1)), r["jf"], rel_tol=1e-6)
        p = os.path.join(d, "spread.dat")
        a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
        healthy = bool(st["rc"] == 0 and a.size > 0 and a[-1, 1] >= 5.9)
        chi = CHI[r["set"]]
        c[lab] = dict(r, missing=False, prov=prov, healthy=healthy and prov, rc=st["rc"],
                      t=float(a[-1, 1]) if a.size else 0.0,
                      rg=float(a[:, 11].max()) if a.size else math.nan,
                      s2=sig(a, r["u0"])[1] if healthy else math.nan,
                      theta=((1 - r["x"]) / (1 - r["x"] / chi)) ** N_EXP,
                      type="DST" if r["x"] >= XDST[r["set"]] else "CST")
    return c


def main():
    c = load()
    say("grade_dsttype -- PREREG_DSTTYPE_2026-09-26.md")
    n_land = sum(not v["missing"] for v in c.values())
    n_ok = sum(v["healthy"] for v in c.values())
    say(f"cells registered {len(c)}, landed {n_land}, healthy+provenance {n_ok}")
    for lab, v in sorted(c.items()):
        if not v["missing"] and not v["healthy"]:
            say(f"   EXCLUDED {lab}: rc {v.get('rc')} t {v['t']:.2f} prov {v['prov']}")
    say("\n  label              type  Theta   Rg_max   S2")
    for lab, v in sorted(c.items(), key=lambda kv: (kv[1]["u0"], kv[1]["set"], kv[1]["x"])):
        if v["healthy"]:
            say(f"  {lab:18s} {v['type']:4s}  {v['theta']:6.1f}  {v['rg']:.4f}  {v['s2']:+.3f}")

    def rg(s, x, u):
        v = c.get(f"{s}_x{x:.4f}_U{u:.1f}")
        return v["rg"] if v and v["healthy"] else math.nan

    say("\n=== D1 (primary): matched Theta, different type ===")
    d1 = []
    for u in (2.0, 4.0):
        d16 = math.log(rg("C5", 0.8925, u)) - math.log(rg("C0", 0.8842, u))
        d10 = math.log(rg("C5", 0.8835, u)) - math.log(rg("C0", 0.8732, u))
        if math.isnan(d16) or math.isnan(d10):
            d1.append(None)
            say(f"   U0 {u}: INCOMPLETE")
            continue
        v = "TYPE" if (d16 > 0.10 and d16 > 3 * abs(d10)) else ("THETA" if abs(d16) < 0.05 else "UNDETERMINED")
        d1.append(v)
        say(f"   U0 {u}: Delta16 (C5 CST - C0 DST) = {d16:+.4f}; Delta10 (both CST) = {d10:+.4f} -> {v}")
    say(f"D1: {'INCOMPLETE' if None in d1 else (d1[0] if d1[0] == d1[1] else 'SPLIT ' + '/'.join(d1))}")

    say("\n=== D2: a jump at x_DST (step across > 2x the largest same-type step) ===")
    passed, n = 0, 0
    xs = {"C5": (0.8735, 0.8835, 0.8925, 0.9035, 0.9105), "C0": (0.8542, 0.8642, 0.8732, 0.8842, 0.8912)}
    for s in ("C5", "C0"):
        for u in (2.0, 4.0):
            lr = [math.log(rg(s, x, u)) for x in xs[s]]
            if any(math.isnan(y) for y in lr):
                say(f"   {s} U0 {u}: INCOMPLETE")
                continue
            steps = np.diff(lr)
            across = steps[2]                       # x3 (CST) -> x4 (DST)
            same = max(abs(steps[0]), abs(steps[1]), abs(steps[3]))
            ok = abs(across) > 2 * same
            n += 1
            passed += ok
            say(f"   {s} U0 {u}: steps {np.round(steps, 4).tolist()}; across {across:+.4f} vs 2x max same-type {2*same:.4f} -> {'jump' if ok else 'no jump'}")
    say(f"D2: {'INCOMPLETE' if n < 4 else ('PASS' if passed >= 3 else 'FAIL')} ({passed}/{n})")

    say("\n=== D3: recoil S2 (DST <= 0.2, CST >= 0.3) ===")
    ok = [v for v in c.values() if v["healthy"] and ((v["type"] == "DST" and v["s2"] <= 0.2) or (v["type"] == "CST" and v["s2"] >= 0.3))]
    say(f"D3: {len(ok)}/{n_ok} healthy cells on their registered side")

    say("\n=== D4 (report): Rg monotone decreasing in Theta across all cells, per U0? ===")
    for u in (2.0, 4.0):
        v = sorted([x for x in c.values() if x["healthy"] and x["u0"] == u], key=lambda x: x["theta"])
        r = [x["rg"] for x in v]
        say(f"   U0 {u}: n {len(v)}  monotone {all(np.diff(r) <= 0)}  (Theta, Rg): "
            + ", ".join(f"({x['theta']:.1f}, {x['rg']:.3f})" for x in v))
    with open(os.path.join(HERE, "out", "grade_dsttype.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
