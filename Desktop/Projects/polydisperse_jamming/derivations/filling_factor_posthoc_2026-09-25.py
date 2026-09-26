"""filling_factor_posthoc_2026-09-25.py -- POST HOC, existing records only (no new data).

Question: does the size distribution act on the jamming packing fraction as ONE Voronoi filling
factor w, independent of the contact constraints? In the map phi = Z/(Z + 2 sqrt3 w)
(Song, Wang & Makse 2008, doi:10.1038/nature06981, w = 1 for monodisperse spheres), invert per cell

    w = Z (1 - phi) / (2 sqrt3 phi),     Z = z_nr (non-rattler mechanical coordination).

If size and constraints separate, w depends on (family, CV) and NOT on friction mu_p.
Registered for new data as PREREG_ROLLCHI R6; this script is the exploratory look at the
records already on disk and is labelled POST HOC wherever quoted.

Loader: analyse_phase2.load_all (deduplicated), converged sphere records only, N printed.
"""
from __future__ import annotations

import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from analyse_phase2 import load_all  # noqa: E402

SQ3 = 2 * math.sqrt(3)
L = []


def say(s=""):
    print(s)
    L.append(s)


recs = load_all()
say(f"N records loaded (dedup): {len(recs)}")
recs = [r for r in recs if not r.get("failed") and r.get("converged")]
say(f"N converged: {len(recs)}")
recs = [r for r in recs if r.get("alpha", 1.0) == 1.0 and r.get("eta_roll") is None
        and r.get("doublet_w", 0) == 0 and r.get("triplet_w", 0) == 0
        and r.get("z_nr") is not None and np.isfinite(r["z_nr"])]
say(f"N plain FrictMat spheres with finite z_nr: {len(recs)}")

g = defaultdict(list)
for r in recs:
    k = (r.get("family", "lognormal"), round(r["cv"], 3), r["mu_p"], r["p_hat_target"], r["n_requested"])
    g[k].append((r["phi"], r["z_nr"]))

say("\nper cell group: w = Z(1-phi)/(2 sqrt3 phi)  [mean over seeds]; SWM residual at w = 1")
say("family      CV    rung    N     mu 0: phi  Z   w     | mu .3: phi  Z   w     | mu 3: phi  Z   w     | w spread")
spread = []
rows = defaultdict(dict)
for (fam, cv, mu, rung, n), v in g.items():
    phi = np.mean([x[0] for x in v])
    z = np.mean([x[1] for x in v])
    w = np.mean([x[1] * (1 - x[0]) / (SQ3 * x[0]) for x in v])
    rows[(fam, cv, rung, n)][mu] = (phi, z, w, len(v))
for k in sorted(rows):
    d = rows[k]
    if 0.0 not in d or len(d) < 2:
        continue
    cells = []
    for mu in (0.0, 0.3, 3.0):
        if mu in d:
            phi, z, w, n = d[mu]
            cells.append(f"{phi:.4f} {z:.2f} {w:.3f}")
        else:
            cells.append("   --            ")
    ws = [d[mu][2] for mu in d]
    sp_ = max(ws) - min(ws)
    spread.append((k, sp_, d))
    say(f"{k[0]:11s} {k[1]:.2f} {k[2]:.0e} {k[3]:5d}  " + " | ".join(cells) + f" | {sp_:.3f}")

say("\nsummary: w spread across friction at fixed (family, CV, rung, N)")
s = np.array([x[1] for x in spread])
say(f"   N groups {len(s)}: median {np.median(s):.3f}, max {s.max():.3f}")
# the size effect on w itself, mu 0, lognormal, for scale
base = {k: d for k, _, d in spread if k[0] == "lognormal"}
say("   scale: lognormal mu-0 w by CV (should fall with CV if CV = better filling):")
for k in sorted(base):
    if 0.0 in base[k]:
        say(f"      CV {k[1]:.2f} rung {k[2]:.0e} N {k[3]}: w {base[k][0.0][2]:.3f}")

say("\nRUNG DRIFT (referee pass 1): w against p_hat at fixed (family, CV, mu)")
for fam in ("lognormal", "uniform"):
    for cv in (0.0, 0.5):
        for mu in (0.0, 0.3):
            pts = [(k[2], d[mu][2]) for k, d in rows.items() if k[0] == fam and k[1] == cv and mu in d]
            if pts:
                pts.sort()
                say(f"   {fam:9s} CV {cv} mu {mu}: " + "  ".join(f"{p:.0e}: {w:.3f}" for p, w in pts))
say("   -> (G) is a jamming-point relation; w drifts upward with p_hat, so every comparison is at FIXED rung")

say("\nSEPARABILITY AS A PREDICTION: w from the mu-0 cell of each group, measured Z at mu 0.3 -> phi_m")
res = []
for k, _, d in spread:
    if 0.0 in d and 0.3 in d:
        w0 = d[0.0][2]
        z3 = d[0.3][1]
        pred = z3 / (z3 + SQ3 * w0)
        res.append((k, d[0.3][0] - pred))
rr = np.array([x[1] for x in res])
say(f"   N groups {len(rr)}: phi_m(measured) - phi_m(predicted): mean {rr.mean():+.4f}, sd {rr.std(ddof=1):.4f}, "
    f"max |.| {np.abs(rr).max():.4f}; within 2 PHI_SCALE (0.0042): {np.mean(np.abs(rr) < 0.0042):.2f}")
for fam in ("lognormal", "uniform", "bidisperse"):
    v = [x[1] for x in res if x[0][0] == fam]
    say(f"      {fam}: n {len(v)} mean {np.mean(v):+.4f} max |.| {np.max(np.abs(v)):.4f}")

# ---------------------------------------------------------------- SHAPE
import glob  # noqa: E402
import json  # noqa: E402

say("\nSHAPE (existing kmer / shapechi records, N 1700, CV 0): the same inversion")
sg = defaultdict(list)
for f in glob.glob(os.path.join(os.path.dirname(HERE), "out", "shape", "phase2_*_lognormal*.jsonl")):
    for line in open(f, encoding="utf-8"):
        r = json.loads(line)
        if r.get("failed") or not r.get("converged") or r.get("tag") not in ("kmer", "shapechi"):
            continue
        if r.get("triplet_w", 0) > 0 or r["n_requested"] != 1700 or r["cv"] != 0.0:
            continue
        lab = f"doublet w2 {r['doublet_w']:.2f}" if r["tag"] == "kmer" else f"clump alpha {r['alpha']:.1f}"
        sg[(lab, r["p_hat_target"], r["mu_p"])].append((r["phi"], r["z_nr"]))
say("   grain                rung    mu 0: phi   Z     w     | mu .3: phi   Z     w     | Z_m/Z_0  chi   chi from (Z_0, Z_m, w_0)")
for lab in sorted({k[0] for k in sg}):
    for rung in (1e-3, 3e-5):
        d = {}
        for mu in (0.0, 0.3):
            v = sg.get((lab, rung, mu))
            if v:
                phi_ = np.mean([x[0] for x in v])
                z_ = np.mean([x[1] for x in v])
                d[mu] = (phi_, z_, z_ * (1 - phi_) / (SQ3 * phi_))
        if len(d) < 2:
            continue
        chi_ = d[0.3][0] / d[0.0][0]
        w0, z0, zm = d[0.0][2], d[0.0][1], d[0.3][1]
        chi_pred = zm * (z0 + SQ3 * w0) / (z0 * (zm + SQ3 * w0))
        say(f"   {lab:20s} {rung:.0e}  {d[0.0][0]:.4f} {z0:.2f} {w0:.3f} | {d[0.3][0]:.4f} {zm:.2f} {d[0.3][2]:.3f} | "
            f"{zm / z0:.3f}   {chi_:.4f} {chi_pred:.4f}")

with open(os.path.join(HERE, "filling_factor_posthoc_2026-09-25.txt"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
