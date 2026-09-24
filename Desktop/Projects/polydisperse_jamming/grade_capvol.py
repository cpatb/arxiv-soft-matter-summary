"""grade_capvol.py -- PREREG_CAPVOL_2026-09-24.md. Writes out/grade_capvol.txt.

Input: out/capvol/<label>/ (scripts/pull_jamwin.sh capvol 2026-09-24), the registered list
out/capvol_cells.txt, and the original runs out/cvmapE03/, out/kmermapE03/. Written before any cell.
spread.dat columns (0-based): 1 t/tau, 4 h/D0, 7 jam_frac, 8 umax/U0, 10 zc/D0, 11 Rg_r/D0, 20 u_rms_jam.
"""
import json
import math
import os

import numpy as np
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
T_IMPACT = 0.5
J_ON = 0.01
TOL = 0.05
DIAG_COLS = (7, 20)      # the only columns JAM_FHAT may change


def load(path, u0):
    st = json.load(open(os.path.join(path, "RUN_STATUS.json"), encoding="utf-8"))
    a = np.loadtxt(os.path.join(path, "spread.dat"), comments="#", ndmin=2)
    healthy = bool(st["rc"] == 0 and a.size and a[-1, 1] >= 5.9)
    late = a[:, 1] >= T_IMPACT
    adm = healthy and a[late, 8].max() <= 1.0 and a[late, 10].max() <= 0.5 and a[late, 4].max() <= 1.0
    return a, healthy and adm


def orig_dir(lab):
    camp, rest = lab.split("_", 1)
    return os.path.join(OUT, {"cv": "cvmapE03", "km": "kmermapE03"}[camp], rest) if camp in ("cv", "km") else None


def main():
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    reg = [ln.split() for ln in open(os.path.join(OUT, "capvol_cells.txt"), encoding="utf-8")
           if ln.strip() and not ln.startswith("#")]
    say("grade_capvol -- PREREG_CAPVOL_2026-09-24.md")
    rows, g0, g1 = [], [], []
    for f in reg:
        lab, phi, phim, u0, jf, jf_orig = f[0], float(f[1]), float(f[5]), float(f[6]), float(f[4]), float(f[9])
        x, chi = phi / 0.62, phim / 0.62
        inside = x >= chi
        od = orig_dir(lab)
        if jf >= 1.0:                       # never reaches the cap: original run, V = 0
            a, ok = load(od, u0)
            V = 0.0
            src = "orig"
        else:
            d = os.path.join(OUT, "capvol", lab)
            if not os.path.exists(os.path.join(d, "RUN_STATUS.json")):
                say(f"  {lab}: MISSING")
                continue
            a, ok = load(d, u0)
            src = "rerun"
            imp = a[:, 1] <= T_IMPACT
            V = float(a[imp, 7].max())
            if od:
                o, _ = load(od, u0)
                keep = [c for c in range(a.shape[1]) if c not in DIAG_COLS]
                same = a.shape == o.shape and np.array_equal(a[:, keep], o[:, keep])
                g0.append((lab, same))
        if od:
            ao, _ = load(od, u0)
            J = float(ao[:, 7].max()) if abs(jf_orig - ((1 - x) / (1 - chi) if inside else 1.0)) < 1e-3 else math.nan
        else:
            J = 0.0
        cls = "J+" if (inside and J > J_ON) else "J0" if not inside else "J?"
        if cls == "J+" and src == "rerun":
            g1.append((lab, V >= J - 1e-12, V, J))
        rows.append(dict(lab=lab, u0=u0, cls=cls, V=V, J=J, ok=ok, lrg=float(np.log(a[:, 11].max())), src=src))

    say(f"G0 wiring (non-diagnostic columns identical to the original): {sum(s for _, s in g0)}/{len(g0)}  "
        + ("PASS" if g0 and all(s for _, s in g0) else f"FAIL {[l for l, s in g0 if not s]}"))
    say(f"G1 V >= J on J+ reruns: {sum(s for _, s, _, _ in g1)}/{len(g1)}  "
        + ("PASS" if g1 and all(s for _, s, _, _ in g1) else f"FAIL {[(l, v, j) for l, s, v, j in g1 if not s]}"))
    say("")
    say("  label                   U0  class  src    healthy  V (capped)  J(orig)  ln Rg_max")
    for r in sorted(rows, key=lambda r: (r["u0"], r["V"])):
        say(f"  {r['lab']:22s}  {r['u0']:g}   {r['cls']}     {r['src']:5s}  {str(r['ok']):5s}    "
            f"{r['V']:.4f}      {r['J']:.3f}    {r['lrg']:+.4f}")
    ok = [r for r in rows if r["ok"] and r["cls"] in ("J+", "J0")]
    say(f"\nhealthy+admissible+classed: {len(ok)}/{len(rows)}")
    h3, res = [], []
    for u0 in sorted({r["u0"] for r in ok}):
        g = [r for r in ok if r["u0"] == u0]
        if len(g) < 6:
            say(f"  U0 {u0:g}: n {len(g)} < 6, H3 not graded")
        else:
            rho = spearmanr([r["V"] for r in g], [r["lrg"] for r in g]).statistic
            h3.append(rho <= -0.7)
            say(f"  U0 {u0:g}: H3 Spearman rho(V, ln Rg_max) = {rho:+.3f} (n {len(g)})  "
                + ("PASS" if rho <= -0.7 else "FAIL"))
        jp = [r for r in g if r["cls"] == "J+"]
        j0 = [r for r in g if r["cls"] == "J0" and r["V"] > 0]
        if len(jp) < 2:
            continue
        lo, hi = min(r["V"] for r in jp), max(r["V"] for r in jp)
        iso = IsotonicRegression(increasing=False).fit([r["V"] for r in jp], [r["lrg"] for r in jp])
        for r in j0:
            if lo <= r["V"] <= hi:
                rr = r["lrg"] - float(iso.predict([r["V"]])[0])
                res.append(rr)
                say(f"     H4 J0-at-cap {r['lab']}: V {r['V']:.4f} in J+ range [{lo:.4f}, {hi:.4f}]  residual {rr:+.4f}")
            else:
                say(f"     H4 J0-at-cap {r['lab']}: V {r['V']:.4f} OUTSIDE J+ range [{lo:.4f}, {hi:.4f}] (not graded)")
    say(f"H3: {sum(h3)}/{len(h3)} groups  " + ("PASS" if h3 and all(h3) else "FAIL" if h3 else "NOT GRADED"))
    med = float(np.median(res)) if res else math.nan
    say(f"H4: {len(res)} interior J0-at-cap cells, median residual {med:+.4f}  "
        + ("NOT INFORMATIVE (< 3)" if len(res) < 3 else "PASS" if abs(med) <= TOL else
           "FAIL (HJ side: unjammed spread more)" if med > TOL else "FAIL (below)"))
    with open(os.path.join(OUT, "grade_capvol.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
