"""grade_cvmap.py -- PREREG_CVMAP_2026-09-24.md. Writes out/grade_cvmap.txt.

Input: out/cvmap/<label>/ from `scripts/pull_jamwin.sh cvmap 2026-09-24` (PowerShell), and the
registered cell list out/cvmap_cells.txt (provenance is checked against it, per cell, from the
solver's own banner). Written before any cell landed. Health (C5) is graded SEPARATELY from the
bars: a dead cell is excluded and printed, never counted as a FAIL or as a jam.
spread.dat columns: 1 t/tau, 2 beta, 4 h/D0, 7 jam_frac (J at JAM_FHAT), 11 Rg_r/D0.
"""
import json
import math
import os
import re

import numpy as np

from grade_jamsig import sig

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "cvmap")
CELLS = os.path.join(HERE, "out", "cvmap_cells.txt")
Q = 1 / 64              # one beta quantum at L9
S2_STAR = 0.1           # pre-fixed in PREREG_JAMSIG Amendment 1, unchanged
J_ON = 0.01             # "jams" (kmermap R1 convention)
X0S = ("0.920", "0.932", "0.950")
CVS = ("0.0", "0.2", "0.5")
US = ("2.0", "4.0")


def load():
    reg = {}
    for line in open(CELLS, encoding="utf-8"):
        f = line.split()
        if f:
            reg[f[0]] = dict(phi=float(f[1]), d=float(f[2]), ss=float(f[3]), jf=float(f[4]),
                             pm=float(f[5]), u0=float(f[6]))
    c = {}
    for lab, r in reg.items():
        d = os.path.join(ROOT, lab)
        if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            c[lab] = dict(missing=True, healthy=False, prov=False)
            continue
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        m = re.search(r"U0=([\d.]+) .*phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+) sigma\*=([\d.e+]+) "
                      r"Pa delta=([\d.]+)", ban)
        jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
        prov = bool(m and jb) and all(
            math.isclose(float(v), t, rel_tol=1e-6)
            for v, t in zip(m.groups(), (r["u0"], r["phi"], 0.62, r["pm"], r["ss"], r["d"]))) \
            and math.isclose(float(jb.group(1)), r["jf"], rel_tol=1e-6)
        p = os.path.join(d, "spread.dat")
        a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
        healthy = bool(st["rc"] == 0 and a.size > 0 and a[-1, 1] >= 5.9)
        c[lab] = dict(missing=False, prov=prov, healthy=healthy, rc=st["rc"],
                      t=float(a[-1, 1]) if a.size else 0.0, inside=r["jf"] < 1.0,
                      J=float(a[:, 7].max()) if a.size else math.nan,
                      beta=float(a[:, 2].max()) if a.size else math.nan,
                      rg=float(a[:, 11].max()) if a.size else math.nan,
                      s2=sig(a, r["u0"])[1] if healthy else math.nan)
    return c


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    c = load()
    ok = {k: v for k, v in c.items() if v["healthy"] and v["prov"]}
    say("grade_cvmap -- PREREG_CVMAP_2026-09-24.md")
    say(f"cells landed {sum(not v['missing'] for v in c.values())}/{len(c)}; provenance "
        f"{sum(v['prov'] for v in c.values())}; healthy+provenance {len(ok)}")
    say("  label               window  healthy  J        beta_max  Rg_max   S2")
    for lab in c:
        v = c[lab]
        if v["missing"]:
            say(f"  {lab:18s}  MISSING")
            continue
        tag = "" if lab in ok else f"   [EXCLUDED rc {v['rc']} t/tau {v['t']:.3g}; prov {v['prov']}]"
        say(f"  {lab:18s}  {'in ' if v['inside'] else 'out'}     {str(v['healthy']):5s}    "
            f"{v['J']:.4f}   {v['beta']:.4f}    {v['rg']:.4f}   {v['s2']:+.3f}{tag}")

    # C1 -- by construction, recorded
    c05 = [c[f"c0.5_x{x}_U{u}"] for x in X0S for u in US if not c[f"c0.5_x{x}_U{u}"]["missing"]]
    say("")
    say(f"C1 (by construction): CV 0.5 max J = {max((v['J'] for v in c05), default=math.nan):.4g} over {len(c05)} landed cells "
        "(outside the window, JAM_FHAT 1.0)")
    c0j = [(x, u, ok[f"c0.0_x{x}_U{u}"]["J"]) for x in X0S for u in US if f"c0.0_x{x}_U{u}" in ok]
    say(f"   CV 0 healthy cells and their J: {[(x, u, round(j, 4)) for x, u, j in c0j]}")

    # C2 -- spreading
    say("")
    dec, good = 0, 0
    for x in X0S:
        for u in US:
            a, b, m = (f"c0.0_x{x}_U{u}", f"c0.5_x{x}_U{u}", f"c0.2_x{x}_U{u}")
            if a not in ok or b not in ok:
                say(f"  C2 x0 {x} U0 {u}: undecided (CV 0 healthy {a in ok}, CV 0.5 healthy {b in ok})")
                continue
            dec += 1
            s = ok[b]["rg"] > ok[a]["rg"] and ok[b]["beta"] >= ok[a]["beta"] - Q
            chain = ""
            if m in ok:
                mono = ok[a]["rg"] <= ok[m]["rg"] <= ok[b]["rg"]
                s = s and mono
                chain = f"  Rg(0.2) {ok[m]['rg']:.4f}  monotone {mono}"
            good += s
            say(f"  C2 x0 {x} U0 {u}: Rg {ok[a]['rg']:.4f} -> {ok[b]['rg']:.4f} "
                f"(ln {math.log(ok[b]['rg'] / ok[a]['rg']):+.4f}); beta {ok[a]['beta']:.4f} -> "
                f"{ok[b]['beta']:.4f}{chain}  -> {'ok' if s else 'VIOLATED'}")
    say(f"C2: {good}/{dec} decided pairs  " + ("PASS" if dec >= 2 and good == dec else
                                             "FAIL" if good < dec else "NOT INFORMATIVE (< 2 decided)"))

    # C3 -- channel split
    say("")
    c3 = []
    for x in X0S:
        need = [f"c0.0_x{x}_U4.0", f"c0.5_x{x}_U4.0", f"W_x{x}_U4.0", f"A_x{x}_U4.0"]
        if not all(k in ok for k in need):
            say(f"  C3 x0 {x}: undecided (healthy: {[k for k in need if k in ok]})")
            continue
        r0, r5, rw, ra = (ok[k]["rg"] for k in need)
        D = math.log(r5 / r0)
        share = math.log(rw / r0) / D
        act = abs(math.log(ra / r0)) / abs(D)
        pa, pb = 0.8 <= share <= 1.2, act <= 0.2
        conv = None
        if f"D_x{x}_U4.0" in ok:
            conv = abs(math.log(ok[f"D_x{x}_U4.0"]["rg"] / r5)) / abs(D)
        pc = conv is not None and conv <= 0.2
        c3.append(pa and pb and pc)
        say(f"  C3 x0 {x}: Delta ln Rg {D:+.4f}; (a) window share {share:.3f} [{pa}]; (b) activation "
            f"|share| {act:.3f} [{pb}]; (c) d32 convention |share| "
            + (f"{conv:.3f} [{pc}]" if conv is not None else "undecided (D dead)"))
    say(f"C3: {sum(c3)}/{len(c3)} decided x0  " + ("PASS" if c3 and all(c3) else
                                                  "FAIL" if c3 else "INCOMPLETE"))

    # C4 -- apex recoil
    say("")
    jp = [(k, v["s2"]) for k, v in ok.items() if v["J"] > J_ON]
    j0 = [(k, v["s2"]) for k, v in ok.items() if v["J"] <= J_ON]
    right = [s < S2_STAR for _, s in jp] + [s >= S2_STAR for _, s in j0]
    say(f"C4 J+ (n {len(jp)}): {[(k, round(s, 3)) for k, s in jp]}")
    say(f"   J0 (n {len(j0)}): S2 range [{min(s for _, s in j0):.3f}, {max(s for _, s in j0):.3f}]"
        if j0 else "   J0: none")
    wrong = [k for (k, s) in jp if s >= S2_STAR] + [k for (k, s) in j0 if s < S2_STAR]
    say(f"C4: {sum(right)}/{len(right)} on the predicted side of S2* = {S2_STAR}  "
        + ("NOT INFORMATIVE (a class is empty)" if not jp or not j0 else
           "PASS" if all(right) else f"FAIL (wrong side: {wrong})"))
    # near-miss J0 cells (0 < J <= J_ON inside a window) are printed, not re-classified
    grey = [(k, round(v["J"], 4)) for k, v in ok.items() if v["inside"] and 0 < v["J"] <= J_ON]
    if grey:
        say(f"   note: in-window cells with 0 < J <= {J_ON} (counted J0 as registered): {grey}")

    say("")
    say(f"C5 health: {len(ok)}/{len(c)} healthy with provenance; excluded: "
        f"{[k for k in c if k not in ok]}")
    with open(os.path.join(HERE, "out", "grade_cvmap.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
