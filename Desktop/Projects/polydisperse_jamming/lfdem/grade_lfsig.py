"""grade_lfsig.py -- grades PREREG_LFSIG_2026-09-24.md (+ Amendment 1). Written before any cell landed.

    bash lfdem/sync_lfsig.sh        (in WSL: pulls raw outputs to out/lfsig/)
    python lfdem/grade_lfsig.py     -> out/grade_lfsig.txt

Per cell (out/lfsig/<label>/): META (label, CV, PHI, SIG, rc), conf.json (realised CV), data_*.dat
(LF_DEM columns: 2 strain, 4 eta_r, 16 min gap, 19 contact number per particle).
Window: strain [1.5, 4] (Amendment 1); 4 equal blocks -> mean and block s.e.
Health (registered): rc == 0, strain reached >= 3, max overlap (-min gap in window) < 0.05.
A cell with rc 124 (wall cap) is CENSORED under the registered rule; a post-hoc variant that admits
rc 124 cells whose window is complete (strain >= 4) is printed separately and labelled post hoc.
"""
from __future__ import annotations

import glob
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")
W0, W1 = 1.5, 4.0
SIGMAS = (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
CVS = (0.2, 0.5)
XMS = (0.90, 0.95)
R_YADE = 0.735 / 0.957             # P1 point prediction (predict_cvmap.txt, d10 convention)
P1_TOL = math.log(1.3)
P2_TOL = 0.02
P4_TOL = math.log(1.2)
MAX_OVERLAP = 0.05
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def load_cell(d: str) -> dict | None:
    meta = open(os.path.join(d, "META"), encoding="utf-8").read()
    kv = dict(re.findall(r"(\w+)=(\S+)", meta))
    files = glob.glob(os.path.join(d, "data_*.dat"))
    if not files:
        return None
    a = np.loadtxt(files[0], comments="#", ndmin=2)
    conf = json.load(open(os.path.join(d, "conf.json"), encoding="utf-8"))
    rc = int(kv["rc"]) if "rc" in kv else None
    if rc is None:  # Amendment 4: runner killed; read LF_DEM's own completion line
        logp = os.path.join(d, "run.log")
        done = os.path.exists(logp) and "Job done ok" in open(logp, encoding="utf-8", errors="replace").read()[-400:]
        rc = 0 if done else 124
    strain = a[:, 1]
    w = (strain >= W0) & (strain <= W1)
    c = {"label": kv["LABEL"], "cv": float(kv["CV"]), "phi": float(kv["PHI"]), "sig": float(kv["SIG"]),
         "rc": rc, "smax": float(strain.max()), "cv_real": conf["cv_real"], "n_win": int(w.sum())}
    c["xm"] = 0.90 if "_x90_" in c["label"] else 0.95
    if w.sum() >= 8:
        blocks = np.array_split(np.where(w)[0], 4)
        for key, col in (("eta", 3), ("zc", 18)):
            bm = np.array([a[b, col].mean() for b in blocks])
            c[key], c[key + "_se"] = float(bm.mean()), float(bm.std(ddof=1) / 2.0)
        c["overlap"] = float(-a[w, 15].min())
    else:
        c["eta"] = c["zc"] = c["overlap"] = float("nan")
    c["complete"] = c["smax"] >= W1 - 0.05
    # Amendment 3: window complete, overlap bounded, rc 0 or rc 124 (wall cap after the window)
    c["healthy_reg"] = rc in (0, 124) and c["complete"] and c["overlap"] < MAX_OVERLAP
    # the original rule, printed for the record: rc 0 and strain >= 3
    c["healthy_ph"] = (rc == 0) and c["smax"] >= 3.0 and c["overlap"] < MAX_OVERLAP
    return c


def sigma_eta(cells: dict, cv: float, xm: float, key: str) -> tuple[float, str]:
    pts = [(s, cells.get((cv, xm, s))) for s in SIGMAS]
    pts = [(s, c) for s, c in pts if c is not None and c[key]]
    if len(pts) < len(SIGMAS):
        return float("nan"), f"INCOMPLETE ({len(pts)}/{len(SIGMAS)} healthy)"
    ls = np.log([s for s, _ in pts])
    le = np.log([c["eta"] for _, c in pts])
    se = np.array([c["eta_se"] / c["eta"] for _, c in pts])
    for i in range(len(le) - 1):
        if le[i + 1] < le[i] - 2 * math.hypot(se[i], se[i + 1]):
            return float("nan"), "UNDEFINED (non-monotone ln eta)"
    mid = 0.5 * (le[0] + le[-1])
    for i in range(len(le) - 1):
        if (le[i] - mid) * (le[i + 1] - mid) <= 0 and le[i + 1] != le[i]:
            t = (mid - le[i]) / (le[i + 1] - le[i])
            return float(math.exp(ls[i] + t * (ls[i + 1] - ls[i]))), "ok"
    return float("nan"), "UNDEFINED (no crossing)"


def yade_cq_ratio() -> float:
    """C and 1/q ratios for CV 0.5 and 0.2 vs 0 at mu 0.3, p 3e-4, parsed from the geometry output."""
    txt = open(os.path.join(OUT, "sigma_star_geometry.txt"), encoding="utf-8").read()
    r = {}
    for cv in ("0.2", "0.5"):
        m = re.search(rf"^\s+0\.3 3e-04\s+{cv}\s+\d+\s+\S+\s+=\s+\S+\s+x\s+(\S+)\s+x\s+(\S+)", txt, re.M)
        r[cv] = float(m.group(1)) * float(m.group(2))
    return r["0.5"] / r["0.2"]


def main() -> int:
    say("grade_lfsig -- PREREG_LFSIG_2026-09-24.md (+ Amendment 1)")
    cells_all = [c for c in (load_cell(d) for d in sorted(glob.glob(os.path.join(OUT, "lfsig", "c*"))) if ".killed_" not in d)
                 if c is not None]
    say(f"cells with data: {len(cells_all)} of {len(CVS) * len(XMS) * len(SIGMAS)} registered")
    say("  label            CVreal  rc    strain  n_win  eta (se)            Zc (se)          overlap  reg  ph")
    for c in sorted(cells_all, key=lambda c: (c["cv"], c["xm"], c["sig"])):
        say(f"  {c['label']:16s} {c['cv_real']:.3f}  {str(c['rc']):5s} {c['smax']:6.2f}  {c['n_win']:5d}  "
            f"{c['eta']:9.3f} ({c.get('eta_se', float('nan')):.3f})  {c['zc']:6.3f} ({c.get('zc_se', float('nan')):.3f})  "
            f"{c['overlap']:.4f}  {int(c['healthy_reg'])}    {int(c['healthy_ph'])}")
    for key, name in (("healthy_reg", "REGISTERED health rule (Amendment 3)"), ("healthy_ph", "ORIGINAL rule (rc 0, strain >= 3), for the record")):
        cells = {(c["cv"], c["xm"], c["sig"]): c for c in cells_all if c[key]}
        say(f"\n==== {name}: {len(cells)} healthy ====")
        ss = {}
        for cv in CVS:
            for xm in XMS:
                ss[(cv, xm)] = sigma_eta(cells, cv, xm, key)
                say(f"  sigma*_eta CV {cv} x_m {xm}: {ss[(cv, xm)][0]:.4g}  [{ss[(cv, xm)][1]}]")
        say(f"\n  P1: R = sigma*(0.5)/sigma*(0.2), Yade point {R_YADE:.3f}, PASS iff R < 1 and |ln R - ln {R_YADE:.3f}| < ln 1.3 at BOTH x_m")
        ok1 = []
        for xm in XMS:
            a, b = ss[(0.5, xm)][0], ss[(0.2, xm)][0]
            if math.isfinite(a) and math.isfinite(b):
                R = a / b
                p = R < 1 and abs(math.log(R / R_YADE)) < P1_TOL
                ok1.append(p)
                say(f"    x_m {xm}: R = {R:.3f}  ln(R/{R_YADE:.3f}) = {math.log(R / R_YADE):+.3f}  -> {'pass' if p else 'fail'}")
            else:
                ok1.append(None)
                say(f"    x_m {xm}: not computable")
        say(f"  P1 verdict: {'INCOMPLETE' if None in ok1 else ('PASS' if all(ok1) else 'FAIL')}")
        # P2: WC inversion
        say("\n  P2: WC inversion phi_J = phi/(1 - eta^-1/2); phi_0 from sigma 0.1, phi_m from sigma 100")
        chi = {}
        for cv in CVS:
            p0, pm = [], []
            for xm in XMS:
                lo, hi = cells.get((cv, xm, 0.1)), cells.get((cv, xm, 100.0))
                if lo and hi:
                    p0.append(lo["phi"] / (1 - lo["eta"] ** -0.5))
                    pm.append(hi["phi"] / (1 - hi["eta"] ** -0.5))
            if len(p0) == 2:
                agree = abs(p0[0] - p0[1]) < 0.01 and abs(pm[0] - pm[1]) < 0.01
                chi[cv] = (np.mean(pm) / np.mean(p0), agree)
                say(f"    CV {cv}: phi_0 {p0[0]:.4f}/{p0[1]:.4f}  phi_m {pm[0]:.4f}/{pm[1]:.4f}  "
                    f"two-phi agreement {'yes' if agree else 'NO'}  chi_LF {chi[cv][0]:.4f}")
            else:
                say(f"    CV {cv}: incomplete")
        if len(chi) == 2 and all(v[1] for v in chi.values()):
            d = chi[0.5][0] - chi[0.2][0]
            say(f"  P2 verdict: dchi = {d:+.4f} -> {'PASS' if abs(d) < P2_TOL else 'FAIL'} (|dchi| < {P2_TOL})")
        elif len(chi) == 2:
            say("  P2 verdict: NOT USABLE (the two phi disagree: the prefactor-1 WC form fails, inversion not used)")
        else:
            say("  P2 verdict: INCOMPLETE")
        # P3
        say("\n  P3: Zc(sigma 0.1) < 1 at x_m 0.90, and eta rises sigma 0.1 -> 100 at every (CV, x_m)")
        ok3 = []
        for cv in CVS:
            lo = cells.get((cv, 0.90, 0.1))
            if lo:
                ok3.append(lo["zc"] < 1)
                say(f"    CV {cv} x_m 0.90: Zc(0.1) = {lo['zc']:.3f}")
            else:
                ok3.append(None)
            for xm in XMS:
                lo, hi = cells.get((cv, xm, 0.1)), cells.get((cv, xm, 100.0))
                if lo and hi:
                    ok3.append(hi["eta"] > lo["eta"])
                    say(f"    CV {cv} x_m {xm}: eta {lo['eta']:.2f} -> {hi['eta']:.2f}")
                else:
                    ok3.append(None)
        say(f"  P3 verdict: {'INCOMPLETE' if None in ok3 else ('PASS' if all(ok3) else 'FAIL')}")
        # P4
        cq = yade_cq_ratio()
        say(f"\n  P4: R_LF vs [Zc phi]_0.5/[Zc phi]_0.2 x (C/q ratio, Yade mu .3 p 3e-4 = {cq:.4f}); PASS iff within x1.2 at both x_m")
        say("      (the registered text rounded the C/q ratio to 1.05; it also said 'recomputed exactly in the grader', which is what is used)")
        ok4 = []
        for xm in XMS:
            z = {}
            for cv in CVS:
                hs = [cells.get((cv, xm, s)) for s in (30.0, 100.0)]
                if all(hs):
                    z[cv] = np.mean([h["zc"] for h in hs]) * hs[0]["phi"]
            a, b = ss[(0.5, xm)][0], ss[(0.2, xm)][0]
            if len(z) == 2 and math.isfinite(a) and math.isfinite(b):
                pred = z[0.5] / z[0.2] * cq
                R = a / b
                p = abs(math.log(R / pred)) < P4_TOL
                ok4.append(p)
                say(f"    x_m {xm}: predicted {pred:.3f} (Zc phi ratio {z[0.5] / z[0.2]:.3f})  measured {R:.3f}  -> {'pass' if p else 'fail'}")
            else:
                ok4.append(None)
                say(f"    x_m {xm}: not computable")
        say(f"  P4 verdict: {'INCOMPLETE' if None in ok4 else ('PASS' if all(ok4) else 'FAIL')}")
        # P5 (Amendment 2): literature anchor, printed not gated
        say("\n  P5 anchor (descriptive; Mari 2014 onset ~0.3 sigma_0; Singh 2018 phi_J^0 0.646, phi_J^inf 0.562)")
        for xm in XMS:
            base = cells.get((0.2, xm, 0.1))
            if not base:
                say(f"    CV 0.2 x_m {xm}: no sigma 0.1 cell")
                continue
            onset = next((s for s in SIGMAS[1:] if (0.2, xm, s) in cells
                          and cells[(0.2, xm, s)]["eta"] > 1.1 * base["eta"]), None)
            p0 = base["phi"] / (1 - base["eta"] ** -0.5)
            hi = cells.get((0.2, xm, 100.0))
            pm = hi["phi"] / (1 - hi["eta"] ** -0.5) if hi else float("nan")
            say(f"    CV 0.2 x_m {xm}: onset sigma {onset} (in [0.1, 1]: {onset is not None and onset <= 1})  "
                f"phi_0 {p0:.4f} (|-0.646| {abs(p0 - 0.646):.4f})  phi_m {pm:.4f} (in [0.562, 0.646]: {0.562 <= pm <= 0.646})")
    with open(os.path.join(OUT, "grade_lfsig.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
