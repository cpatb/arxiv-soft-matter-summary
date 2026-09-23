"""grade_jamsig.py -- PREREG_JAMSIG_2026-09-23.md. Writes out/grade_jamsig.txt.

Post-processing only: camera-visible candidate signatures of jamming (vs thickening) on existing
healthy L9 DNS cells. spread.dat columns: 1 t/tau, 4 h/D0, 7 jam_frac (J), 11 Rg_r/D0, 19 u_rms.
Written with the prereg, before any observable was computed.
"""
import glob
import json
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
BAR = 1.6


def cells():
    """Yield (campaign, label, U0, x, spread array) for every healthy cell with a readable banner."""
    for camp in ("kmermap", "kmerdecomp", "kmerdrop", "kmerdrop2"):
        for d in sorted(glob.glob(os.path.join(OUT, camp, "*"))):
            lab = os.path.basename(d)
            if lab.startswith("_ref") or not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
                continue
            st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
            p = os.path.join(d, "spread.dat")
            if not os.path.exists(p):
                continue
            a = np.loadtxt(p, comments="#", ndmin=2)
            if not (st["rc"] == 0 and a.size and a[-1, 1] >= 5.9):
                continue
            ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
            m = re.search(r"U0=([\d.]+) .*phi=([\d.]+) phi0=([\d.]+)", ban)
            if not m:
                continue
            u0, phi, phi0 = (float(g) for g in m.groups())
            yield camp, lab, u0, round(phi / phi0, 4), a


def sig(a, u0):
    t, h, rg, ur = a[:, 1], a[:, 4], a[:, 11], a[:, 19]
    late = (t >= 5.0) & (t <= 6.0)
    s1 = float(ur[late].mean() / u0)
    i0 = int(np.argmin(h))
    s2 = float(h[i0:].max() - h[i0])
    r3 = float(np.interp(3.0, t, rg))
    s3 = float((np.interp(6.0, t, rg) - r3) / r3)
    return s1, s2, s3


def separates(jp, j0, name, say):
    """Non-overlapping ranges AND median ratio >= BAR (either direction)."""
    if not jp or not j0:
        return None
    mp, m0 = float(np.median(jp)), float(np.median(j0))
    overlap = not (max(jp) < min(j0) or min(jp) > max(j0))
    ratio = (max(abs(mp), abs(m0)) / min(abs(mp), abs(m0))) if min(abs(mp), abs(m0)) > 0 else float("inf")
    ok = (not overlap) and ratio >= BAR
    say(f"    {name}: J+ median {mp:+.4g} [{min(jp):+.4g}, {max(jp):+.4g}] (n {len(jp)})  "
        f"J0 median {m0:+.4g} [{min(j0):+.4g}, {max(j0):+.4g}] (n {len(j0)})  ratio {ratio:.2f}  "
        f"overlap {overlap}  -> {'SEPARATES (J+ ' + ('lower' if mp < m0 else 'higher') + ')' if ok else 'not informative'}")
    return ok, mp < m0


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("grade_jamsig -- PREREG_JAMSIG_2026-09-23.md")
    rows = []
    for camp, lab, u0, x, a in cells():
        J = float(a[:, 7].max())
        s = sig(a, u0)
        rows.append(dict(camp=camp, lab=lab, u0=u0, x=x, J=J, s=s, rg=float(a[:, 11].max())))
        say(f"  {camp:10s} {lab:16s} U0 {u0:.0f} x {x:.4f} J {J:.4f} Rg {rows[-1]['rg']:.4f}  "
            f"S1 {s[0]:.4g}  S2 {s[1]:+.4g}  S3 {s[2]:+.4g}")
    names = ("S1 late motion", "S2 recoil", "S3 late creep")
    verdict = {n: [] for n in names}
    say("")
    for u0 in sorted({r["u0"] for r in rows}):
        g = [r for r in rows if r["u0"] == u0]
        jp = [r for r in g if r["J"] > 0.01]
        j0 = [r for r in g if r["J"] == 0.0]
        if not jp or not j0:
            continue
        say(f"  U0 {u0:.0f}: {len(jp)} J+ vs {len(j0)} J0 (all campaigns; x differs between groups)")
        for i, n in enumerate(names):
            r = separates([x["s"][i] for x in jp], [x["s"][i] for x in j0], n, say)
            if r is not None:
                verdict[n].append(r)
    say("")
    say("  SHARED-x subset (kmerdrop ladder, x 0.85, U0 2):")
    kd = [r for r in rows if r["camp"] == "kmerdrop"]
    jp, j0 = [r for r in kd if r["J"] > 0.01], [r for r in kd if r["J"] == 0.0]
    shared = {}
    for i, n in enumerate(names):
        shared[n] = separates([x["s"][i] for x in jp], [x["s"][i] for x in j0], n, say)
    say("")
    s1 = verdict[names[0]]
    say("J1 S1 separates at every U0 with J+ LOWER: "
        + ("PASS" if s1 and all(ok and low for ok, low in s1) else "FAIL" if s1 else "INCOMPLETE"))
    say("J2 S2 and S3 do NOT separate at any U0: "
        + ("PASS" if all(not ok for n in names[1:] for ok, _ in verdict[n]) else "FAIL"))
    for n in names:
        sep_all = bool(verdict[n]) and all(ok for ok, _ in verdict[n])
        sh = shared[n]
        tag = ("JAMMING SIGNATURE (separates across U0 and at shared x)" if sep_all and sh and sh[0]
               else "x-confounded (separates across x only)" if sep_all else "no signature")
        say(f"J3 {n}: {tag}")
    # ---- Amendment 1 (registered before its cells): held-out K1 and grid K2, S2* = 0.1 -------
    say("")
    say("Amendment 1 (S2* = 0.1 fixed before the cells):")
    am = {}
    for d in sorted(glob.glob(os.path.join(OUT, "jamsig", "*"))):
        lab = os.path.basename(d)
        if lab.startswith("_ref") or not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
            continue
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        p = os.path.join(d, "spread.dat")
        a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        lvl_ok = ("_L10_" in ban) if lab.startswith("L10") else ("_L9_" in ban)
        healthy = bool(st["rc"] == 0 and a.size and a[-1, 1] >= 5.9 and lvl_ok)
        u0 = 4.0 if "U4.0" in lab else 2.0
        J = float(a[:, 7].max()) if a.size else float("nan")
        s2 = sig(a, u0)[1] if a.size else float("nan")
        am[lab] = dict(healthy=healthy, J=J, s2=s2)
        say(f"  {lab}: healthy {healthy} (level tag ok {lvl_ok})  J {J:.4f}  S2 {s2:+.4f}  "
            f"class {'J+' if J > 0.01 else 'J0' if J == 0 else 'ambiguous'}")
    k1 = [(v["J"] > 0.01, v["s2"]) for k, v in am.items() if k.startswith("kd_") and v["healthy"]]
    if k1:
        right = [(s2 < 0.1) if jp else (s2 >= 0.1) for jp, s2 in k1]
        classes = {jp for jp, _ in k1}
        say(f"K1 held-out: {sum(right)}/{len(right)} on the predicted side of 0.1"
            + ("  PASS" if all(right) and len(classes) == 2 else
               "  UNINFORMATIVE (one class only)" if all(right) else "  FAIL"))
    else:
        say("K1 held-out: INCOMPLETE")
    a25, a50 = am.get("L10_w0.25_U2.0"), am.get("L10_w0.50_U2.0")
    if a25 and a50 and a25["healthy"] and a50["healthy"]:
        cls = (a25["J"] == 0.0) and (a50["J"] > 0.01)
        side = (a25["s2"] >= 0.1) and (a50["s2"] < 0.1)
        say(f"K2 L10: classes kept {cls}, S2 sides kept {side}  "
            + ("PASS" if cls and side else "UNDECIDABLE (class changed)" if not cls else "FAIL"))
    else:
        say("K2 L10: INCOMPLETE")
    # ---- Amendment 2 (registered before its cells): singlets across x, K3 ------------------
    say("")
    say("Amendment 2 (singlets across x, default cap; S2* = 0.1):")
    k3 = []
    for d in sorted(glob.glob(os.path.join(OUT, "jamsig2", "s_*"))):
        lab = os.path.basename(d)
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        p = os.path.join(d, "spread.dat")
        a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        em = re.search(r"eps_min=([\d.e+-]+)", ban)
        prov = bool(em) and abs(float(em.group(1)) - 0.01) < 1e-12 and "_L9_" in ban
        healthy = bool(st["rc"] == 0 and a.size and a[-1, 1] >= 5.9 and prov)
        u0 = 4.0 if "U4.0" in lab else 2.0
        J = float(a[:, 7].max()) if a.size else float("nan")
        s2 = sig(a, u0)[1] if a.size and a[-1, 1] >= 5.0 else float("nan")
        c_ = "J+" if J > 0.01 else "J0" if J == 0 else "ambiguous"
        say(f"  {lab}: healthy {healthy} (rc {st['rc']}, t/tau {a[-1, 1] if a.size else 0})  J {J:.4f}  "
            f"S2 {s2:+.4f}  class {c_}")
        if healthy and c_ != "ambiguous":
            k3.append((c_, s2, (s2 < 0.1) if c_ == "J+" else (s2 >= 0.1)))
    classes = {c for c, _, _ in k3}
    say(f"K3: {sum(r for *_, r in k3)}/{len(k3)} on the predicted side; classes present {sorted(classes)}  "
        + ("PASS" if k3 and all(r for *_, r in k3) and len(classes) == 2 else
           "FAIL" if any(not r for *_, r in k3) else "UNINFORMATIVE/INCOMPLETE"))
    # ---- Amendment 3 (registered before its cells): EPS_SOFT 0.002 rescue -------------------
    say("")
    say("Amendment 3 (EPS_SOFT 0.002; S2* = 0.1):")
    resc = ("s_x0.93_U2.0", "s_x0.94_U2.0", "s_x0.94_U4.0")
    t1, t2, k3b = [], [], []
    ctrl = {"s_x0.93_U4.0": ("J+", 0.0), "s_x0.86_U4.0": ("J0", 0.5)}
    for d in sorted(glob.glob(os.path.join(OUT, "jamsig3", "s_*"))):
        lab = os.path.basename(d)
        st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
        p = os.path.join(d, "spread.dat")
        a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
        ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
        prov = "_esft0.002" in ban and "_L9_" in ban
        healthy = bool(st["rc"] == 0 and a.size and a[-1, 1] >= 5.9 and prov)
        u0 = 4.0 if "U4.0" in lab else 2.0
        J = float(a[:, 7].max()) if a.size else float("nan")
        s2 = sig(a, u0)[1] if a.size and a[-1, 1] >= 5.0 else float("nan")
        c_ = "J+" if J > 0.01 else "J0" if J == 0 else "ambiguous"
        say(f"  {lab}: healthy {healthy} (esft tag {prov}, rc {st['rc']}, t/tau {a[-1, 1] if a.size else 0})"
            f"  J {J:.4f}  S2 {s2:+.4f}  class {c_}")
        if lab in resc:
            t1.append(healthy)
            if healthy and c_ != "ambiguous":
                k3b.append((s2 < 0.1) if c_ == "J+" else (s2 >= 0.1))
        if lab in ctrl and healthy:
            want_c, _ = ctrl[lab]
            t2.append(c_ == want_c and ((s2 < 0.1) if want_c == "J+" else (s2 >= 0.1)))
    say(f"T1 rescue: {sum(t1)}/{len(t1)} healthy (bar >= 2 of 3)  " + ("PASS" if sum(t1) >= 2 else
                                                                      "FAIL" if len(t1) == 3 else "INCOMPLETE"))
    say(f"T2 controls keep class and S2 side: {sum(t2)}/{len(t2)}  " + ("PASS" if len(t2) == 2 and all(t2)
                                                                    else "FAIL" if t2 and not all(t2) else "INCOMPLETE"))
    say(f"K3' rescued cells on the predicted side: {sum(k3b)}/{len(k3b)}  "
        + ("PASS" if k3b and all(k3b) else "FAIL" if k3b else "INCOMPLETE"))
    with open(os.path.join(OUT, "grade_jamsig.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
