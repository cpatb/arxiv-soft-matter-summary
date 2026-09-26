"""cvmap_by_type_2026-09-26.py -- POST HOC: classify the 24 cvmap DNS drops (EPS_MIN 0.03) as CST / DST / SJ.

Inputs (read, not typed): out/cvmap_cells.txt (label, PHI, DELTA, SIGMA_STAR, JAM_FHAT, PHIM, U0 per cell) and
out/grade_cvmapE03.txt (J, beta_max, Rg_max, S2, Rg_term per cell). PHI0 = 0.62 and n = MP_EXPO = 2 are the cvmap
campaign's fixed solver inputs (PREREG_CVMAP_2026-09-24; kmerdecomp launcher PHI0=0.62; MP_EXPO default 2).
The type classifier is derivations/thickening_type_2026-09-25_check.x_dst (imported, not re-implemented).
Question: does the steady-shear TYPE order the drop's outcome, and in particular are the "thickened but unjammed"
drops that no DNS observable separates from jammed ones DST-band drops?
Writes derivations/cvmap_by_type_2026-09-26.txt.
"""
from __future__ import annotations

import importlib.util
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location("tt", os.path.join(HERE, "thickening_type_2026-09-25_check.py"))
PHI0, N_EXP = 0.62, 2.0
L = []


def say(s=""):
    print(s)
    L.append(s)


def x_dst_fn():
    import numpy as np
    from scipy.optimize import brentq

    def Hmax(x, chi, nd):
        uu = np.logspace(-6, 3, 40001)
        J = 1 - (1 - chi) * np.exp(-uu)
        return float((nd * x * (1 - chi) * uu * np.exp(-uu) / (J * (J - x))).max())

    def x_dst(chi, nd):
        f = lambda x: Hmax(x, chi, nd) - 1.0  # noqa: E731
        if f(chi * (1 - 1e-9)) < 0:
            return float("nan")
        return brentq(f, 1e-6, chi * (1 - 1e-9), xtol=1e-10)
    return x_dst


x_dst = x_dst_fn()
cells = {}
for line in open(os.path.join(ROOT, "out", "cvmap_cells.txt"), encoding="utf-8"):
    p = line.split()
    if len(p) == 7:
        cells[p[0]] = dict(phi=float(p[1]), delta=float(p[2]), phim=float(p[5]), u0=float(p[6]))
grade = {}
pat = re.compile(r"^\s+(\S+)\s+(in|out)\s+(True|False)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([+-][\d.]+)\s+([\d.]+)")
for line in open(os.path.join(ROOT, "out", "grade_cvmapE03.txt"), encoding="utf-8"):
    m = pat.match(line)
    if m:
        grade[m.group(1)] = dict(J=float(m.group(4)), beta=float(m.group(5)), rg=float(m.group(6)),
                                 s2=float(m.group(7)), rgt=float(m.group(8)), healthy=m.group(3) == "True")
say(f"N cells in cvmap_cells.txt: {len(cells)}; N graded (EPS_MIN 0.03): {len(grade)}; joined: {len(set(cells) & set(grade))}")
say("label               x       chi     n*delta  x_DST    type  | J       Rg_max  S2     Rg_term")
rows = []
for lab in sorted(set(cells) & set(grade)):
    c, g = cells[lab], grade[lab]
    x, chi, nd = c["phi"] / PHI0, c["phim"] / PHI0, N_EXP * c["delta"]
    xd = x_dst(chi, nd)
    t = "SJ" if x >= chi else ("DST" if x >= xd else "CST")
    rows.append((lab, x, chi, nd, xd, t, g))
    say(f"{lab:19s} {x:.4f}  {chi:.4f}  {nd:.3f}    {xd:.4f}   {t:4s}  | {g['J']:.4f}  {g['rg']:.4f}  {g['s2']:+.3f} {g['rgt']:.4f}")
say("\nby type (U0 pooled): n, Rg_max range, S2 range, J range")
for t in ("SJ", "DST", "CST"):
    v = [r for r in rows if r[5] == t]
    if v:
        say(f"   {t}: n {len(v)}  Rg_max [{min(r[6]['rg'] for r in v):.4f}, {max(r[6]['rg'] for r in v):.4f}]  "
            f"S2 [{min(r[6]['s2'] for r in v):+.3f}, {max(r[6]['s2'] for r in v):+.3f}]  J [{min(r[6]['J'] for r in v):.3f}, {max(r[6]['J'] for r in v):.3f}]")
with open(os.path.join(HERE, "cvmap_by_type_2026-09-26.txt"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
