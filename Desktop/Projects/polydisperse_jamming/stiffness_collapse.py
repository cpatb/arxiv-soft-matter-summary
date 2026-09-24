"""stiffness_collapse.py -- PREREG_STIFFNESS_2026-09-24.md. Writes out/stiffness_collapse.txt.

Does the drop respond to the PEAK VISCOSITY it reaches (H), or to whether the closure calls it jammed (H0)?
Post-processing of every admissible L9 / MP_EXPO 2 / EPS_SOFT 0 / theta0 5 / clock 1 cs_fit cell on disk.
The viscosity law is the solver's own (wc_drop_impact.c:487-500):
    eta/eta_s = max(eps, EPS_MIN)^-2,  eps = 1 - x / (1 - (1 - chi) fhat),  x = phi/phi0, chi = phiM/phi0.
spread.dat columns (0-based): 1 t/tau, 4 h/D0, 5 fhat_max, 6 fhat_mean, 7 J, 8 umax/U0, 10 zc/D0, 11 Rg_r/D0.
"""
import glob
import json
import math
import os
import re

import numpy as np
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

from grade_jamsig import sig

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
T_IMPACT = 0.5          # impact phase for S (prereg)
J_ON = 0.01
TOL = 0.05              # |median residual| bar for H2 (prereg)
RHO_BAR = -0.7


def eta_ratio(x, chi, eps_min, f):
    """Solver law, MP_EXPO 2, hard floor (EPS_SOFT 0 cells only)."""
    f = np.clip(f, 0.0, 1.0)
    eps = 1.0 - x / (1.0 - (1.0 - chi) * f)
    return np.maximum(eps, eps_min) ** -2.0


def cells(say):
    n = {"dirs": 0, "banner": 0, "config": 0, "healthy": 0, "admissible": 0, "classed": 0}
    rows = []
    for d in sorted(glob.glob(os.path.join(OUT, "*", "*"))):
        if not os.path.isfile(os.path.join(d, "spread.dat")) or not os.path.isfile(os.path.join(d, "RUN_STATUS.json")):
            continue
        n["dirs"] += 1
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read() \
            if os.path.exists(os.path.join(d, "banner.txt")) else ""
        m = re.search(r"U0=([\d.]+) .*phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+)", ban)
        tag = re.search(r"# tag=(\S+)", ban)
        em = re.search(r"eps_min=([\d.e+-]+)", ban)
        es = re.search(r"EPS_SOFT=([\d.e+-]+)", ban)
        jf = re.search(r"JAM_FHAT=([\d.]+)", ban)
        if not (m and tag and em and es):
            continue
        n["banner"] += 1
        t_ = tag.group(1)
        if not ("preset=cs_fit" in ban and "_L9_" in t_ and "_mpe" not in t_ and float(es.group(1)) == 0
                and "_th5_" in t_ and "_clk1" in t_):
            continue
        n["config"] += 1
        u0, phi, phi0, phim = (float(g) for g in m.groups())
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        a = np.loadtxt(os.path.join(d, "spread.dat"), comments="#", ndmin=2)
        if not (st["rc"] == 0 and a.size and a[-1, 1] >= 5.9):
            continue
        n["healthy"] += 1
        late = a[:, 1] >= T_IMPACT
        if a[late, 8].max() > 1.0 or a[late, 10].max() > 0.5 or a[late, 4].max() > 1.0:
            continue
        n["admissible"] += 1
        x, chi, eps_min = phi / phi0, phim / phi0, float(em.group(1))
        inside = x >= chi
        fc = (1 - x) / (1 - chi) if inside else math.nan
        J = float(a[:, 7].max())
        if inside:
            if not (jf and abs(float(jf.group(1)) - fc) < 1e-3):
                continue
            cls = "J+" if J > J_ON else "J0"
        else:
            cls = "J0"
        n["classed"] += 1
        imp = a[:, 1] <= T_IMPACT
        rows.append(dict(camp=os.path.basename(os.path.dirname(d)), lab=os.path.basename(d), u0=u0, x=x,
                         chi=chi, eps=eps_min, cls=cls, J=J,
                         S=float(np.log(eta_ratio(x, chi, eps_min, a[imp, 5]).max())),
                         Sbar=float(np.log(eta_ratio(x, chi, eps_min, a[imp, 6]).max())),
                         lrg=float(np.log(a[:, 11].max())), s2=sig(a, u0)[1]))
    say("N per filter: " + ", ".join(f"{k} {v}" for k, v in n.items()))
    return rows


def grade(rows, key, say, graded):
    verdicts = []
    for u0 in sorted({r["u0"] for r in rows}):
        g = [r for r in rows if r["u0"] == u0]
        if len(g) < 6:
            say(f"  U0 {u0:g}: n {len(g)} < 6, not graded")
            continue
        rho = spearmanr([r[key] for r in g], [r["lrg"] for r in g]).correlation
        j0 = [r for r in g if r["cls"] == "J0"]
        jp = [r for r in g if r["cls"] == "J+"]
        lo, hi = min(r[key] for r in j0), max(r[key] for r in j0)
        iso = IsotonicRegression(increasing=False).fit([r[key] for r in j0], [r["lrg"] for r in j0])
        inner = [r for r in jp if lo <= r[key] <= hi]
        res = [r["lrg"] - float(iso.predict([r[key]])[0]) for r in inner]
        med = float(np.median(res)) if res else math.nan
        say(f"  U0 {u0:g}: n {len(g)} (J0 {len(j0)}, J+ {len(jp)}); {key} range J0 [{lo:.2f}, {hi:.2f}], "
            f"J+ [{min((r[key] for r in jp), default=math.nan):.2f}, {max((r[key] for r in jp), default=math.nan):.2f}]")
        say(f"     H1 Spearman rho({key}, ln Rg_max) = {rho:+.3f}  " + ("PASS" if rho <= RHO_BAR else "FAIL"))
        say(f"     H2 interior J+ {len(inner)}: median residual {med:+.4f}"
            + (f"  [{min(res):+.4f}, {max(res):+.4f}]" if res else "")
            + ("  -> " + ("PASS" if abs(med) <= TOL else "FAIL (H0 side)" if med < -TOL else "FAIL (above)")
               if len(inner) >= 3 else "  not informative (< 3)"))
        for r in inner:
            say(f"       {r['camp']}/{r['lab']}: {key} {r[key]:.2f}  ln Rg {r['lrg']:+.4f}  "
                f"fit {float(iso.predict([r[key]])[0]):+.4f}  J {r['J']:.3f}  eps_min {r['eps']:g}")
        verdicts.append((rho <= RHO_BAR, (abs(med) <= TOL) if len(inner) >= 3 else None))
    if graded:
        h1 = [v[0] for v in verdicts]
        h2 = [v[1] for v in verdicts if v[1] is not None]
        say(f"H1: {sum(h1)}/{len(h1)} groups  " + ("PASS" if h1 and all(h1) else "FAIL"))
        say(f"H2: {sum(h2)}/{len(h2)} informative groups  " + ("PASS" if h2 and all(h2) else
                                                              "FAIL" if h2 else "NOT INFORMATIVE"))


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("stiffness_collapse -- PREREG_STIFFNESS_2026-09-24.md")
    rows = cells(say)
    say(f"campaigns used: {sorted({r['camp'] for r in rows})}")
    say("")
    say("=== S = ln(eta_peak/eta_s) at fhat_max, t/tau <= 0.5 (GRADED) ===")
    grade(rows, "S", say, True)
    say("")
    say("=== Sbar (fhat_mean), reported only ===")
    grade(rows, "Sbar", say, False)
    say("")
    say("=== every cell: camp/label U0 x chi eps_min class J S Sbar lnRg S2 ===")
    for r in sorted(rows, key=lambda r: (r["u0"], r["S"])):
        say(f"  {r['camp']}/{r['lab']:22s} {r['u0']:g}  {r['x']:.4f} {r['chi']:.4f} {r['eps']:<6g} {r['cls']} "
            f"{r['J']:.3f}  {r['S']:6.3f} {r['Sbar']:6.3f}  {r['lrg']:+.4f}  {r['s2']:+.3f}")
    with open(os.path.join(OUT, "stiffness_collapse.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
