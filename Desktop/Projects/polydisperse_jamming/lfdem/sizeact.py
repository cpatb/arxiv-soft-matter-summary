"""sizeact.py -- POST HOC, descriptive (not a registered gate): WHICH grains carry the frictional contacts as the
stress rises, in lubricated LF_DEM shear.

    python lfdem/sizeact.py [lfsig|lfsig2]   ->  out/sizeact_<campaign>.txt

For each cell (x_m 0.90 arms), interaction snapshots in strain [2, 4]:
  per grain i: k_i = number of touching contacts (state >= 1, i.e. h < 0), time-averaged over snapshots.
  Grains are split into radius QUARTILES of their own configuration (conf.dat).
  Printed per quartile and stress: <k> (mean touching contacts per grain), P(k >= 1) (activated fraction),
  and the quartile's share of all touching-contact ends, against its share of grains (0.25) and of surface area.
The claim this is built to test (Amendment 1 M1, and Paper C's mechanism): in a polydisperse pack the first contacts to
activate belong to the LARGEST grains, and the activated set spreads down the size distribution as the stress rises.
A size-blind activation would give every quartile the same P(k >= 1) at every stress.
The lubricated rattler analogue of `rattler_closure` (mean degree ~ a^p) is fitted per stress: p(sigma) = slope of
ln<k> against ln(mean quartile radius) over the four quartiles (least squares, numpy.polyfit); p -> 0 is size-blind.
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mech_lfsig import snapshots  # noqa: E402  (the same dump parser the mechanism readout uses)

OUT = os.path.join(os.path.dirname(HERE), "out")
W0, W1 = 2.0, 4.0
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def cell(d: str) -> dict | None:
    ints = glob.glob(os.path.join(d, "int_*.dat"))
    if not ints or not os.path.exists(os.path.join(d, "conf.dat")):
        return None
    kv = dict(re.findall(r"(\w+)=(\S+)", open(os.path.join(d, "META"), encoding="utf-8").read()))
    r = np.loadtxt(os.path.join(d, "conf.dat"), comments="#", ndmin=2)[:, 3]
    N = r.size
    ks = []
    for strain, _rate, s in snapshots(ints[0]):
        if not (W0 <= strain <= W1):
            continue
        t = s[:, 2] >= 1
        k = np.bincount(s[t, 0].astype(int), minlength=N) + np.bincount(s[t, 1].astype(int), minlength=N)
        ks.append(k)
    if not ks:
        return None
    K = np.array(ks, dtype=float)             # snapshots x grains
    kbar = K.mean(axis=0)
    pact = (K >= 1).mean(axis=0)
    q = np.quantile(r, [0.25, 0.5, 0.75])
    qi = np.digitize(r, q)                    # 0..3, smallest to largest
    rows = []
    for j in range(4):
        m = qi == j
        rows.append({"rmean": float(r[m].mean()), "n": int(m.sum()), "k": float(kbar[m].mean()),
                     "pact": float(pact[m].mean()), "share_k": float(kbar[m].sum() / max(kbar.sum(), 1e-12)),
                     "share_area": float((r[m] ** 2).sum() / (r**2).sum())})
    kk = np.array([x["k"] for x in rows]); rr = np.array([x["rmean"] for x in rows])
    p = float(np.polyfit(np.log(rr), np.log(kk), 1)[0]) if np.all(kk > 0) else float("nan")
    return {"label": kv["LABEL"], "cv": float(kv["CV"]), "sig": float(kv["SIG"]), "seed": int(kv.get("SEED", 1)),
            "nsnap": len(ks), "rows": rows, "p": p, "zc": float(kbar.mean())}


def collapse(cells: list) -> None:
    """Does a grain's activation depend on the stress only through sigma * a^m?  The Derjaguin argument: the load a
    stress sigma puts on a contact of grains of size a scales as sigma a^2, the threshold as R* ~ a, so load/threshold
    ~ sigma a (m = 1); m = 0 is size-blind, m = 2 would be load-without-threshold-scaling.  Points: one per (cell,
    quartile), abscissa ln(sigma * r_q^m), ordinate logit P (P clipped to [0.005, 0.995], points with P outside
    (0.02, 0.98) dropped since logit is uninformative there).  Collapse quality: rms residual about ONE least-squares
    line (numpy.polyfit, degree 1) through all points; the m that minimises it is scanned on a grid (0 to 3, step 0.05)
    and also reported with the residual at m = 0, 1, 2.  Pooled over CVs and seeds, so a CV-independent m is part of
    what is being tested; per-CV scans are printed beside it."""
    # activation RELATIVE TO EACH QUARTILE'S OWN PLATEAU (sigma 100 of the same CV and seed): small grains rattle even on
    # the plateau (rattling ~ a^2), so raw P mixes activation with rattling; P/P(100) isolates the activation.
    plat = {(c["cv"], c["seed"]): c for c in cells if c["sig"] == 100.0}
    pts = []
    for c in cells:
        pc = plat.get((c["cv"], c["seed"]))
        if pc is None or c["sig"] == 100.0:
            continue
        for j, x in enumerate(c["rows"]):
            P = x["pact"] / pc["rows"][j]["pact"]
            if 0.02 < P < 0.98:
                pts.append((c["cv"], c["sig"], x["rmean"], P))
    if len(pts) < 8:
        say("\n  collapse: too few points")
        return
    arr = np.array(pts)

    def rms(m, sel):
        a = arr[sel]
        X = np.log(a[:, 1] * a[:, 2] ** m)
        Pc = np.clip(a[:, 3], 0.005, 0.995)
        Y = np.log(Pc / (1 - Pc))
        coef = np.polyfit(X, Y, 1)
        return float(np.sqrt(np.mean((Y - np.polyval(coef, X)) ** 2))), float(coef[0])

    say("\n  COLLAPSE of activation onto sigma * a^m (logit P vs ln sigma a^m, one line; rms residual, lower is better)")
    grid = np.round(np.arange(-1.0, 5.0001, 0.05), 2)
    groups = [("pooled", np.ones(len(arr), bool))] + [(f"CV {cv:.2f}", arr[:, 0] == cv) for cv in sorted(set(arr[:, 0]))]
    for name, sel in groups:
        if sel.sum() < 6:
            continue
        res = [rms(m, sel)[0] for m in grid]
        mb = float(grid[int(np.argmin(res))])
        r0, r1, r2 = rms(0.0, sel)[0], rms(1.0, sel)[0], rms(2.0, sel)[0]
        say(f"    {name:8s} n {int(sel.sum()):3d}: best m {mb:4.2f} (rms {min(res):.3f}, slope {rms(mb, sel)[1]:.2f}); "
            f"rms at m 0 / 1 / 2 = {r0:.3f} / {r1:.3f} / {r2:.3f}")


def main() -> int:
    camp = sys.argv[1] if len(sys.argv) > 1 else "lfsig"
    say(f"sizeact -- POST HOC, descriptive: size-resolved activation in LF_DEM [{camp}], strain [{W0}, {W1}], x_m 0.90 arms")
    say("  quartiles Q1 (smallest) .. Q4 (largest) of each configuration's radii; <k> touching contacts per grain;")
    say("  P = fraction of time a grain has >= 1 touching contact; share = the quartile's share of contact ends (area share in brackets);")
    say("  p = slope of ln<k> vs ln(quartile mean radius): 0 is size-blind, 2 is 'contacts follow surface area'")
    dirs = [x for x in sorted(glob.glob(os.path.join(OUT, camp, "c*"))) if ".killed_" not in x and "_x95_" not in x and "c50x" not in x]
    cells = [c for c in (cell(d) for d in dirs) if c]
    for c in sorted(cells, key=lambda c: (c["cv"], c["seed"], c["sig"])):
        rw = c["rows"]
        say(f"  {c['label']:18s} Zc {c['zc']:5.2f} p {c['p']:5.2f} | " +
            "  ".join(f"Q{j + 1} r{x['rmean']:.2f} <k>{x['k']:.2f} P{x['pact']:.2f} sh{x['share_k']:.2f}({x['share_area']:.2f})"
                      for j, x in enumerate(rw)))
    # the headline contrast: at each CV, the ratio P(Q4)/P(Q1) as the stress rises (size selectivity of activation)
    say("\n  size selectivity P(Q4)/P(Q1) and exponent p vs stress, per CV (mean over seeds present)")
    for cv in sorted({c["cv"] for c in cells}):
        sub = [c for c in cells if c["cv"] == cv]
        for sig in sorted({c["sig"] for c in sub}):
            g = [c for c in sub if c["sig"] == sig]
            sel = [c["rows"][3]["pact"] / c["rows"][0]["pact"] for c in g if c["rows"][0]["pact"] > 0]
            ps = [c["p"] for c in g if np.isfinite(c["p"])]
            say(f"    CV {cv:4.2f} sigma {sig:6.1f}: n {len(g)}  P(Q4)/P(Q1) {np.mean(sel) if sel else float('nan'):7.2f}"
                f"  p {np.mean(ps) if ps else float('nan'):5.2f}  Zc {np.mean([c['zc'] for c in g]):5.2f}")
    collapse(cells)
    outp = os.path.join(OUT, f"sizeact_{camp}.txt")
    with open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
