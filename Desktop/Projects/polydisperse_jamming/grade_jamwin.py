"""grade_jamwin.py -- grades PREREG_JAMWIN_2026-09-22.md. Writes out/grade_jamwin.txt.

Written before the grid cells landed. Input: out/jamwin/<label>/{spread.dat, RUN_STATUS.json,
banner.txt}, pulled by scripts/pull_jamwin.sh.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "out", "jamwin")
ARMS = {"A": dict(phi=0.589, fc=0.25, jf=0.25), "B": dict(phi=0.527, fc=0.75, jf=0.75),
        "N": dict(phi=0.465, fc=None, jf=0.5)}
PHI0, PHIM = 0.62, 0.496
SIGS = (30, 300, 3000)
DELTAS = ("0.7", "1.3169", "2.6")
T_END_MIN = 5.9
LO, HI = 1e-3, 0.99                 # informative-pair bounds on J, registered
COLS = ("t t_tau beta beta_wet h fhat_max fhat_mean jam_frac umax vol zc Rg_r Rg_z AR "
        "A_if zc_if Rgs_r Rgs_z AR_if u_rms u_rms_jam").split()


def md5(p):
    with open(p, "rb") as fh:
        return hashlib.md5(fh.read()).hexdigest()


def load(lab):
    d = os.path.join(ROOT, lab)
    st = json.load(open(os.path.join(d, "RUN_STATUS.json"), encoding="utf-8"))
    dat = os.path.join(d, "spread.dat")
    arr = np.loadtxt(dat, comments="#", ndmin=2) if os.path.exists(dat) else None
    if arr is not None and arr.size == 0:
        arr = None                       # header-only file: no rows were written
    ban = open(os.path.join(d, "banner.txt"), encoding="utf-8", errors="replace").read() \
        if os.path.exists(os.path.join(d, "banner.txt")) else ""
    return st, arr, ban


def main() -> int:
    out = []

    def say(s=""):
        print(s)
        out.append(s)

    say("grade_jamwin -- PREREG_JAMWIN_2026-09-22.md")
    # ---- G0 ---------------------------------------------------------------------
    ref = os.path.join(ROOT, "_ref_DS_d1.3169_L9", "spread.dat")
    for g in ("G0a_unset", "G0b_half"):
        p = os.path.join(ROOT, g, "spread.dat")
        if os.path.exists(p) and os.path.exists(ref):
            same = md5(p) == md5(ref)
            say(f"G0 {g}: md5 {md5(p)[:12]} vs deltastiff {md5(ref)[:12]}  "
                f"{'PASS' if same else 'FAIL'}")
        else:
            say(f"G0 {g}: missing")

    # ---- load grid ----------------------------------------------------------------
    cells = {}
    for arm, a in ARMS.items():
        for s in SIGS:
            for d in DELTAS:
                lab = f"{arm}_s{s}_d{d}"
                if not os.path.exists(os.path.join(ROOT, lab, "RUN_STATUS.json")):
                    continue
                st, arr, ban = load(lab)
                # G1 provenance
                m = re.search(r"phi=([\d.]+) phi0=([\d.]+) phiM=([\d.]+) sigma\*=([\d.e+]+) "
                              r"Pa delta=([\d.]+)", ban)
                prov = bool(m) and all(math.isclose(float(v), w, rel_tol=1e-9) for v, w in
                                       zip(m.groups(), (a["phi"], PHI0, PHIM, s, float(d))))
                jb = re.search(r"JAM_FHAT=([\d.]+)", ban)
                prov &= bool(jb) and math.isclose(float(jb.group(1)), a["jf"])
                healthy = st["rc"] == 0 and arr is not None and arr[-1, 1] >= T_END_MIN
                c = {"arm": arm, "s": s, "d": float(d), "prov": prov, "healthy": healthy,
                     "st": st}
                if arr is not None and len(arr):
                    col = {k: arr[:, i] for i, k in enumerate(COLS) if i < arr.shape[1]}
                    c["J"] = float(col["jam_frac"].max())
                    c["t_J"] = float(col["t_tau"][np.argmax(col["jam_frac"])])
                    c["fmean"] = float(col["fhat_mean"].max())
                    c["beta"] = float(col["beta"].max())
                    late = col["t_tau"] >= 5.0
                    c["jam_late"] = float(col["jam_frac"][late].mean()) if late.any() else math.nan
                    c["h_end"] = float(col["h"][-1])
                cells[lab] = c
    say(f"cells landed: {len(cells)} of 27")
    say(f"G1 PROVENANCE: {sum(c['prov'] for c in cells.values())}/{len(cells)} banners match "
        f"phi/phi0/phiM/sigma*/delta/JAM_FHAT")
    for k, c in cells.items():
        if not c["prov"]:
            say(f"  PROVENANCE FAIL {k}")
    bad = [k for k, c in cells.items() if not c["healthy"]]
    say(f"G2 HEALTH: {len(cells) - len(bad)}/{len(cells)} rc=0 and t/tau >= {T_END_MIN}")
    for k in bad:
        c = cells[k]["st"]
        say(f"  EXCLUDED {k}: rc={c['rc']} last t/tau={c['last_t_over_tau']}")

    say("")
    say("arm  sigma*  delta    J=max jam_frac  t/tau@J  <fhat>max  beta_max  jam(t/tau>=5)")
    for k, c in cells.items():
        if "J" in c:
            say(f"  {c['arm']}  {c['s']:5d}  {c['d']:6.4f}   {c['J']:.4f}         "
                f"{c['t_J']:.2f}     {c['fmean']:.3f}     {c['beta']:.3f}    "
                f"{c['jam_late']:.4f}{'' if c['healthy'] and c['prov'] else '   [EXCLUDED]'}")

    # ---- predictions ------------------------------------------------------------------
    def series(arm, s):
        row = [cells.get(f"{arm}_s{s}_d{d}") for d in DELTAS]
        if any(r is None or not (r["healthy"] and r["prov"]) or "J" not in r for r in row):
            return None
        return [r["J"] for r in row]

    say("")
    verdict = {}
    for arm, want in (("A", "down"), ("B", "up")):
        inf_pairs, passes = 0, 0
        for s in SIGS:
            J = series(arm, s)
            if J is None:
                say(f"  {arm} sigma* {s}: incomplete")
                continue
            informative = not (all(j <= LO for j in J) or all(j >= HI for j in J))
            if not informative:
                say(f"  {arm} sigma* {s}: J = {[round(j, 4) for j in J]}  NOT INFORMATIVE")
                continue
            inf_pairs += 1
            if want == "down":
                ok = J[0] >= J[1] >= J[2] and J[0] > J[2]
            else:
                ok = J[0] <= J[1] <= J[2] and J[2] > J[0]
            passes += ok
            say(f"  {arm} sigma* {s}: J(delta 0.7/1.32/2.6) = {[round(j, 4) for j in J]}  "
                f"predicted {'non-increasing' if want == 'down' else 'non-decreasing'}  "
                f"{'PASS' if ok else 'FAIL'}")
        verdict[arm] = (inf_pairs, passes)
        tag = "P1" if arm == "A" else "P2"
        say(f"{tag} arm {arm}: {passes}/{inf_pairs} informative sigma* pass"
            + ("  PASS" if inf_pairs and passes == inf_pairs else
               "  UNDETERMINED (observable unmeasurable: no informative pair among cells "
               "passing G1/G2)" if not inf_pairs else "  FAIL"))
    a, b = verdict.get("A", (0, 0)), verdict.get("B", (0, 0))
    p3 = a[0] and b[0] and a[1] == a[0] and b[1] == b[0]
    say(f"P3 reversal: {'PASS' if p3 else 'NOT ESTABLISHED'}")
    # P4
    n30 = [cells.get(f"N_s30_d{d}") for d in DELTAS]
    if all(c is not None and "J" in c for c in n30):
        # only cells that pass G1 and G2 count (reviewer 2026-09-22: a dead cell's pre-death
        # rows must not demonstrate anything)
        ok4 = all(c["healthy"] and c["prov"] and c["J"] > 0.1 for c in n30)
        say(f"P4 arm N at sigma* 30 (x < chi, cannot jam): J(fhat>0.5) = "
            f"{[round(c['J'], 3) for c in n30]}  {'DEMONSTRATED' if ok4 else 'NOT DEMONSTRATED'}"
            f"   gap floor 1 - x/chi = {1 - 0.75 / 0.8:.4f} > EPS_MIN 0.01")
    with open(os.path.join(HERE, "out", "grade_jamwin.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
