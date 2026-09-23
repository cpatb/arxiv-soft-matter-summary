"""grade_kmerdecomp.py -- PREREG_KMERDECOMP_2026-09-23.md. Writes out/grade_kmerdecomp.txt.
Input: out/kmerdecomp/<label>/ from `scripts/pull_jamwin.sh kmerdecomp 2026-09-23` (PowerShell), and the
kmermap reference arms in out/kmermap_cells.json (written by grade_kmermap.py). Written before the data."""
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "kmerdecomp")
# label: (phi, phim, sigma*, delta, JAM_FHAT, U0)
CELLS = {"win_x0.88_U2.0": (0.5456, 0.53512, 300.0, 1.3169, 0.8766, 2.0),
         "win_x0.88_U4.0": (0.5456, 0.53512, 300.0, 1.3169, 0.8766, 4.0),
         "win_x0.85_U4.0": (0.527, 0.53512, 300.0, 1.3169, 1.0, 4.0),
         "act_x0.85_U4.0": (0.527, 0.56315, 235.8, 1.1163, 1.0, 4.0),
         # Amendment 1
         "act_x0.85_U1.0": (0.527, 0.56315, 235.8, 1.1163, 1.0, 1.0),
         "act_x0.85_U2.0": (0.527, 0.56315, 235.8, 1.1163, 1.0, 2.0)}


def read(lab, par):
    d = os.path.join(ROOT, lab)
    if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
        return None
    st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
    ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read()
    m = re.search(r"U0=([\d.]+) .*phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+) sigma\*=([\d.e+]+) "
                  r"Pa delta=([\d.]+)", ban)
    jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
    phi, pm, ss, dl, jf, u0 = par
    prov = bool(m and jb) and all(math.isclose(float(v), t, rel_tol=1e-6) for v, t in zip(
        m.groups(), (u0, phi, 0.62, pm, ss, dl))) and math.isclose(float(jb.group(1)), jf, rel_tol=1e-6)
    p = os.path.join(d, "spread.dat")
    a = np.loadtxt(p, comments="#", ndmin=2) if os.path.exists(p) else np.zeros((0, 21))
    healthy = bool(st["rc"] == 0 and a.size > 0 and a[-1, 1] >= 5.9)
    return dict(prov=prov, healthy=healthy, rc=st["rc"], t=st["last_t_over_tau"],
                J=float(a[:, 7].max()) if a.size else math.nan,
                rg=float(a[:, 11].max()) if a.size else math.nan,
                beta=float(a[:, 2].max()) if a.size else math.nan)


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    say("grade_kmerdecomp -- PREREG_KMERDECOMP_2026-09-23.md")
    ref = {(r["w"], r["x"], r["U0"]): r for r in
           json.load(open(os.path.join(HERE, "out", "kmermap_cells.json"), encoding="utf-8"))}
    c = {lab: read(lab, par) for lab, par in CELLS.items()}
    for lab, v in c.items():
        if v is None:
            say(f"  {lab}: not run yet")
            continue
        flag = "" if v["healthy"] and v["prov"] else f"   [EXCLUDED rc {v['rc']} t {v['t']} prov {v['prov']}]"
        say(f"  {lab}: J {v['J']:.4f}  Rg_max {v['rg']:.4f}  beta_max {v['beta']:.4f}{flag}")
    ok = {k: v for k, v in c.items() if v and v["healthy"] and v["prov"]}
    bound = {k: v for k, v in c.items() if v and v["prov"] and not v["healthy"]}

    say("")
    d1 = [ok[k]["J"] for k in ("win_x0.88_U2.0", "win_x0.88_U4.0") if k in ok]
    d1b = [bound[k]["J"] for k in ("win_x0.88_U2.0", "win_x0.88_U4.0") if k in bound]
    say(f"D1 window alone jams at x 0.88 (J > 0.01 in >= 1 of U0 2/4): healthy J {[round(j, 4) for j in d1]}"
        f"{'' if not d1b else f', dead-cell lower bounds {[round(j, 4) for j in d1b]}'}  "
        + ("PASS" if any(j > 0.01 for j in d1 + d1b) else "FAIL" if len(d1) == 2 else "INCOMPLETE"))
    for u in (2.0, 4.0):
        k = f"win_x0.88_U{u}"
        a = ref.get(("w1", 0.88, u))
        if k in ok and a and a["healthy"] and a["prov"]:
            say(f"D2 x 0.88 U0 {u}: J(window only) - J(all channels) = {ok[k]['J'] - a['J']:+.4f}"
                f"  ({ok[k]['J']:.4f} vs {a['J']:.4f})")
        else:
            say(f"D2 x 0.88 U0 {u}: not decidable (a cell is missing, dead or off-provenance)")
    w0 = ref.get(("w0", 0.85, 4.0))
    al = ref.get(("w1", 0.85, 4.0))
    win, act = ok.get("win_x0.85_U4.0"), ok.get("act_x0.85_U4.0")
    if w0 and al and win and act and w0["healthy"] and al["healthy"]:
        lw, la, lall = (math.log(win["rg"] / w0["rg"]), math.log(act["rg"] / w0["rg"]),
                        math.log(al["rg"] / w0["rg"]))
        say(f"D3 x 0.85 U0 4: ln cut window {lw:+.4f}, activation {la:+.4f}, all {lall:+.4f}  "
            + ("PASS (window cut larger)" if lw < la else "FAIL (activation cut >= window cut)"))
        say(f"D4 additivity: ln(all) {lall:+.4f} vs ln(win)+ln(act) {lw + la:+.4f}"
            f" (difference {lall - (lw + la):+.4f}; no bar)")
    else:
        say("D3/D4: INCOMPLETE (a needed cell is missing, dead or off-provenance)")
    say("")
    say("Amendment 1 (low U0, activation only, x 0.85):")
    e1 = None
    for u in (1.0, 2.0):
        k, r = f"act_x0.85_U{u}", ref.get(("w0", 0.85, u))
        if k in ok and r and r["healthy"]:
            cut = math.log(ok[k]["rg"] / r["rg"])
            say(f"  U0 {u}: ln(R_act/R_w0) = {cut:+.4f} ({ok[k]['rg']:.4f} vs {r['rg']:.4f})")
            if u == 1.0:
                e1 = abs(cut) > 0.0008
        else:
            say(f"  U0 {u}: not decidable (missing, dead or reference unhealthy)")
    say("E1 |cut| at U0 1 > 0.0008: " + ("INCOMPLETE" if e1 is None else "PASS" if e1 else "FAIL"))
    say(f"D5 health: {len(ok)}/{sum(v is not None for v in c.values())}")
    with open(os.path.join(HERE, "out", "grade_kmerdecomp.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
