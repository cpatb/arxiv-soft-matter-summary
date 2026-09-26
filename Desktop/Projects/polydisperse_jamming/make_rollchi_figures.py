"""make_rollchi_figures.py -- figure R1, the geometric window (PREREG_ROLLCHI_2026-09-25).

R1 (a) phi vs Z = z_nr at p_hat 3e-5 for every converged DEM population (seed means, bars = seed sd),
       with phi = Z/(Z + 2 sqrt3 w) drawn at w = 1 (Song-Wang-Makse, monodisperse spheres) and at the
       lognormal CV-0.5 and pure-doublet filling factors measured at mu_p 0 (printed, not chosen).
   (b) window width 1 - chi = 1 - phi(mu_p .3)/phi(mu_p 0) against contact lever arm: rolling spheres at
       their etaRoll, clumps at eta_eff = (pi/4)(alpha - 1)(d_s/d_eq) (derivations/geometric_window C3).
Explanations in FIGURES_ROLLCHI.md; house style from make_figures.py. K4-excluded clump cells are dropped
(the rule and the count are printed).

    python make_rollchi_figures.py [--check]
"""
import argparse
import glob
import json
import math
import os
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

from make_figures import COL, FIG, md5, save   # noqa: E402
from shape import sub_sphere_diameter          # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
RUNG = 3e-5
SQ3 = 2 * math.sqrt(3)
K4_DEPTH, K4_FRAC = 0.05, 1e-3


def k4_ok(r):
    if r.get("doublet_w", 0) > 0:
        return r["mu_p"] < 1.0          # every kmer mu_p 3 multiplet cell is K4-excluded (grade_kmer.txt)
    if r.get("alpha", 1.0) in (1.0,) or r["cv"] != 0.0:
        return True
    p = os.path.join(OUT, "achi_dump", r["key"] + ".npz")
    if not os.path.exists(p):
        return r["mu_p"] < 1.0          # shapechi alpha 2 at mu <= 0.3 is K4-clean (PREREG_KMER Am. 2)
    fn = np.load(p)["fn"]
    rr = 0.5 * sub_sphere_diameter(1.0, r["alpha"])
    return not len(fn) or float(np.mean(fn / rr ** 2 > K4_DEPTH)) <= K4_FRAC


def load():
    recs, n0 = [], 0
    for f in glob.glob(os.path.join(OUT, "rollchi", "phase2_*_lognormal*.jsonl")) + \
            glob.glob(os.path.join(OUT, "shape", "phase2_*_lognormal*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if r.get("tag") not in ("kmer", "shapechi") and not str(r.get("tag", "")).startswith(("rchi", "achi")):
                continue
            n0 += 1
            if r.get("failed") or not r.get("converged") or r.get("triplet_w", 0) != 0 or r["n_requested"] != 1700:
                continue
            recs.append(r)
    by = {}
    for r in recs:
        by[r["key"]] = r
    recs = list(by.values())
    n1 = len(recs)
    recs = [r for r in recs if k4_ok(r)]
    print(f"R1: N records read {n0}, converged/unique/N1700/no-triplet {n1}, after K4 {len(recs)}")
    return recs


def pop(r):
    if r.get("doublet_w", 0) > 0:
        return f"doublet {r['doublet_w']:.2f}"
    if r.get("alpha", 1.0) != 1.0:
        return f"clump {r['alpha']:.2f}"
    if r.get("eta_roll") not in (None, 0.0):
        return f"roll {r['eta_roll']:.2f}"
    return "sphere"


def build():
    recs = load()
    g = defaultdict(list)
    for r in recs:
        if r["p_hat_target"] == RUNG:
            g[(pop(r), r["cv"], r["mu_p"])].append((r["phi"], r["z_nr"]))
    rows = []
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.8))
    # (a)
    style = {"sphere": (COL[0], "o"), "roll": (COL[1], "s"), "clump": (COL[3], "^"), "doublet": (COL[4], "v")}
    w_meas = {}
    for (p, cv, mu), v in sorted(g.items()):
        v = np.array(v)
        phi, z = v[:, 0].mean(), v[:, 1].mean()
        if mu == 0.0 and (p, cv) in (("sphere", 0.5), ("doublet 1.00", 0.0)):
            w_meas[(p, cv)] = float(np.mean(v[:, 1] * (1 - v[:, 0]) / (SQ3 * v[:, 0])))
        c, m = style[p.split()[0]]
        face = c if cv == 0.0 else "white"
        ax[0].errorbar(z, phi, xerr=v[:, 1].std(ddof=1) if len(v) > 1 else 0,
                       yerr=v[:, 0].std(ddof=1) if len(v) > 1 else 0, fmt=m, ms=3.5, mfc=face, mec=c,
                       ecolor=c, elinewidth=0.6, lw=0)
        rows.append({"panel": "a", "population": p, "cv": cv, "mu_p": mu, "Z": z, "phi": phi, "lever": "", "window": "", "n": len(v)})
    zz = np.linspace(0.01, 9.0, 300)
    for w, ls in ((1.0, "-"), (w_meas.get(("sphere", 0.5), np.nan), "--"), (w_meas.get(("doublet 1.00", 0.0), np.nan), ":")):
        if np.isfinite(w):
            ax[0].plot(zz, zz / (zz + SQ3 * w), ls, color="0.45", lw=0.8)
            ax[0].text(8.9, 8.9 / (8.9 + SQ3 * w) + 0.012, f"w {w:.2f}", ha="right", fontsize=7, color="0.35")
    ax[0].set_xlim(0, 9)
    ax[0].set_ylim(0, 0.75)
    ax[0].set_xlabel("non-rattler coordination $Z$")
    ax[0].set_ylabel(r"$\phi$ at $\hat p = 3\times10^{-5}$")
    # (b) phi_m(alpha) at mu_p 0.3, both rungs; rolling-saturated sphere and pure doublet as references
    def seedmean(sel):
        v = [r["phi"] for r in recs if sel(r)]
        return (float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else 0.0, len(v)) if v else None
    for rung, c, m in ((1e-3, COL[2], "o"), (3e-5, COL[3], "^")):
        pts = []
        for a in (1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 1.9, 1.95, 2.0):
            v = seedmean(lambda r, a=a: r["p_hat_target"] == rung and r["cv"] == 0.0 and r["mu_p"] == 0.3
                         and r.get("doublet_w", 0) == 0 and abs(r.get("alpha", 1.0) - a) < 1e-9
                         and r.get("eta_roll") in (None, 0.0))
            if v:
                pts.append((a,) + v)
                rows.append({"panel": "b", "population": f"clump {a}", "cv": 0.0, "mu_p": 0.3, "Z": "", "phi": v[0],
                             "lever": "", "window": "", "n": v[2]})
        a_, mm, sd, _ = zip(*pts)
        ax[1].errorbar(a_, mm, yerr=sd, fmt=m + "-", color=c, ms=3.5, lw=0.8, elinewidth=0.6)
        roll = seedmean(lambda r: r["p_hat_target"] == rung and r["cv"] == 0.0 and r["mu_p"] == 0.3
                        and r.get("alpha", 1.0) == 1.0 and r.get("eta_roll") == 1.0)
        if roll:
            ax[1].axhline(roll[0], color=c, lw=0.7, ls="--")
    ax[1].set_xlim(0.95, 2.05)
    ax[1].set_ylim(0, 0.7)
    ax[1].set_xlabel(r"aspect ratio $\alpha$")
    ax[1].set_ylabel(r"$\phi_m$ ($\mu_p$ 0.3)")
    # (c) window width 1 - chi(alpha) at p_hat 1e-3 (the rung where both branches exist for every alpha)
    wx, wy = [], []
    for a in (1.0, 1.1, 1.2, 1.3, 1.5, 1.75, 2.0):
        f3 = seedmean(lambda r, a=a: r["p_hat_target"] == 1e-3 and r["cv"] == 0.0 and r["mu_p"] == 0.3
                      and r.get("doublet_w", 0) == 0 and abs(r.get("alpha", 1.0) - a) < 1e-9 and r.get("eta_roll") in (None, 0.0))
        f0 = seedmean(lambda r, a=a: r["p_hat_target"] == 1e-3 and r["cv"] == 0.0 and r["mu_p"] == 0.0
                      and r.get("doublet_w", 0) == 0 and abs(r.get("alpha", 1.0) - a) < 1e-9 and r.get("eta_roll") in (None, 0.0))
        if f3 and f0:
            wx.append(a)
            wy.append(1 - f3[0] / f0[0])
            rows.append({"panel": "c", "population": f"clump {a}", "cv": 0.0, "mu_p": 0.3, "Z": "", "phi": "",
                         "lever": "", "window": wy[-1], "n": min(f3[2], f0[2])})
    ax[2].plot(wx, wy, "o-", color=COL[0], ms=3.5, lw=0.8)
    ax[2].set_xlim(0.95, 2.05)
    ax[2].set_ylim(0, None)
    ax[2].set_xlabel(r"aspect ratio $\alpha$")
    ax[2].set_ylabel(r"window $1-\chi$ ($\hat p = 10^{-3}$)")
    for a_, lab in zip(ax, "abc"):
        a_.text(0.02, 0.97, f"({lab})", transform=a_.transAxes, va="top")
    fig.tight_layout()
    save(fig, "R1_geometric_window", rows,
         "make_rollchi_figures.py; inputs out/rollchi/*.jsonl, out/shape/phase2_{kmer,shapechi}*.jsonl; "
         "p_hat 3e-5, N 1700, seed means (bars = seed sd); K4 applied to clumps")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    build()
    if a.check:
        h1 = md5(os.path.join(FIG, "R1_geometric_window.png"))
        build()
        h2 = md5(os.path.join(FIG, "R1_geometric_window.png"))
        print(f"R1 png md5 {h1} / {h2}: {'IDENTICAL' if h1 == h2 else 'DIFFERS'}")


if __name__ == "__main__":
    main()
