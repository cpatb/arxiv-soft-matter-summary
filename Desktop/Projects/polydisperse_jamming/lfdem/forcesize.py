"""forcesize.py -- POST HOC, descriptive (not registered): how the normal LOAD on a pair scales with its size, in
lubricated LF_DEM shear, and whether that explains the CV-dependent activation exponent m of `sizeact.py`.

    python lfdem/forcesize.py [lfsig2]   ->  out/forcesize_<campaign>.txt

Reading. A pair activates (touches, becomes frictional) when its normal load exceeds its repulsive threshold. With a
Derjaguin threshold F_th ~ R* = a_i a_j/(a_i + a_j) and a load that scales as F ~ sigma R*^beta, activation happens when
sigma R*^(beta - 1) exceeds a constant, so the activation exponent is m = beta - 1 (a constant-force threshold gives
m = beta). Mean-field (load = stress x area) gives beta = 2, m = 1.

Measured: per interaction snapshot in strain [2, 4], every listed pair (all pairs inside lub_max_gap), total normal load
F = lubrication normal (col 8, positive compressive) + contact spring normal (col 12, touching only) + repulsion (col 16).
Only COMPRESSIVE pairs (F > 0) enter. beta = slope of ln<F | R* bin> against ln R*, 8 equal-count bins of R*, numpy
least squares; also reported for the touching subset alone and for the non-touching subset alone (the latter is the
population that has NOT yet activated, i.e. the one whose load decides the next activation).
Compared against sizeact's best m per CV: the prediction m = beta - 1.
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mech_lfsig import snapshots  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "out")
W0, W1 = 2.0, 4.0
NBIN = 8
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def beta_fit(rs: np.ndarray, f: np.ndarray) -> tuple[float, int]:
    m = f > 0
    rs, f = rs[m], f[m]
    if rs.size < 4 * NBIN:
        return float("nan"), int(rs.size)
    edges = np.quantile(rs, np.linspace(0, 1, NBIN + 1))
    idx = np.clip(np.digitize(rs, edges[1:-1]), 0, NBIN - 1)
    xb = np.array([rs[idx == k].mean() for k in range(NBIN)])
    yb = np.array([f[idx == k].mean() for k in range(NBIN)])
    return float(np.polyfit(np.log(xb), np.log(yb), 1)[0]), int(rs.size)


def cell(d: str) -> dict | None:
    ints = glob.glob(os.path.join(d, "int_*.dat"))
    if not ints or not os.path.exists(os.path.join(d, "conf.dat")):
        return None
    kv = dict(re.findall(r"(\w+)=(\S+)", open(os.path.join(d, "META"), encoding="utf-8").read()))
    r = np.loadtxt(os.path.join(d, "conf.dat"), comments="#", ndmin=2)[:, 3]
    RS, F, T = [], [], []
    for strain, _rate, s in snapshots(ints[0]):
        if not (W0 <= strain <= W1):
            continue
        i, j = s[:, 0].astype(int), s[:, 1].astype(int)
        RS.append(r[i] * r[j] / (r[i] + r[j]))
        F.append(s[:, 7] + s[:, 11] + s[:, 15])
        T.append(s[:, 2] >= 1)
    if not RS:
        return None
    rs, f, t = np.concatenate(RS), np.concatenate(F), np.concatenate(T)
    # per-PAIR activation: fraction of listed pairs (inside lub_max_gap) that touch, in 8 equal-count R* bins
    edges = np.quantile(rs, np.linspace(0, 1, NBIN + 1))
    idx = np.clip(np.digitize(rs, edges[1:-1]), 0, NBIN - 1)
    pair_r = np.array([rs[idx == k].mean() for k in range(NBIN)])
    pair_f = np.array([t[idx == k].mean() for k in range(NBIN)])
    b_all, n_all = beta_fit(rs, f)
    b_t, n_t = beta_fit(rs[t], f[t])
    b_n, n_n = beta_fit(rs[~t], f[~t])
    return {"label": kv["LABEL"], "cv": float(kv["CV"]), "sig": float(kv["SIG"]), "seed": int(kv.get("SEED", 1)),
            "b_all": b_all, "b_t": b_t, "b_n": b_n, "n_all": n_all, "n_t": n_t, "n_n": n_n,
            "pair_r": pair_r, "pair_f": pair_f}


def pair_collapse(cells: list) -> None:
    """Per-PAIR activation collapse: logit(fraction of neighbouring pairs touching) against ln(sigma R*^m), one line,
    rms residual scanned over m.  Derjaguin predicts m = 1 (load ~ sigma R*^2, threshold ~ R*), a size-independent
    critical load m = 2.  Unlike the per-grain version in sizeact.py, no neighbour-count multiplicity enters.
    Pairs with fraction outside (0.02, 0.98) dropped (logit uninformative)."""
    pts = [(c["cv"], c["sig"], r, f) for c in cells for r, f in zip(c["pair_r"], c["pair_f"]) if 0.02 < f < 0.98]
    if len(pts) < 8:
        say("\n  pair collapse: too few points")
        return
    arr = np.array(pts)

    def rms(m, sel):
        a = arr[sel]
        X = np.log(a[:, 1] * a[:, 2] ** m)
        Y = np.log(a[:, 3] / (1 - a[:, 3]))
        c = np.polyfit(X, Y, 1)
        return float(np.sqrt(np.mean((Y - np.polyval(c, X)) ** 2)))

    grid = np.round(np.arange(-1.0, 5.0001, 0.05), 2)
    say("\n  PER-PAIR activation collapse onto sigma R*^m (Derjaguin m 1, constant critical load m 2)")
    groups = [("pooled", np.ones(len(arr), bool))] + [(f"CV {cv:.2f}", arr[:, 0] == cv) for cv in sorted(set(arr[:, 0]))]
    for name, sel in groups:
        if sel.sum() < 6:
            continue
        res = [rms(m, sel) for m in grid]
        say(f"    {name:8s} n {int(sel.sum()):3d}: best m {float(grid[int(np.argmin(res))]):4.2f} (rms {min(res):.3f}); "
            f"rms at m 0 / 1 / 2 = {rms(0.0, sel):.3f} / {rms(1.0, sel):.3f} / {rms(2.0, sel):.3f}")


def main() -> int:
    camp = sys.argv[1] if len(sys.argv) > 1 else "lfsig2"
    say(f"forcesize -- POST HOC: load-size exponent beta (F ~ R*^beta, compressive pairs) in LF_DEM [{camp}], strain [{W0}, {W1}]")
    say("  prediction under a Derjaguin threshold: activation exponent m = beta - 1 (mean field beta 2 -> m 1)")
    dirs = [x for x in sorted(glob.glob(os.path.join(OUT, camp, "c*"))) if ".killed_" not in x and "c50x" not in x]
    cells = [c for c in (cell(d) for d in dirs) if c]
    say(f"  {'label':18s} {'b_all':>6s} {'b_touch':>7s} {'b_notyet':>8s}   n_all n_touch n_notyet")
    for c in sorted(cells, key=lambda c: (c["cv"], c["seed"], c["sig"])):
        say(f"  {c['label']:18s} {c['b_all']:6.2f} {c['b_t']:7.2f} {c['b_n']:8.2f}   {c['n_all']:6d} {c['n_t']:7d} {c['n_n']:8d}")
    say("\n  per CV, cells below the plateau (sigma <= 10, where activation happens): mean beta and the implied m = beta - 1")
    for cv in sorted({c["cv"] for c in cells}):
        g = [c for c in cells if c["cv"] == cv and c["sig"] <= 10 and np.isfinite(c["b_n"])]
        if not g:
            continue
        bn = np.mean([c["b_n"] for c in g]); ba = np.mean([c["b_all"] for c in g])
        say(f"    CV {cv:4.2f}: n {len(g):2d} cells  beta(not-yet-touching) {bn:.2f} -> m {bn - 1:.2f};  beta(all) {ba:.2f} -> m {ba - 1:.2f}")
    pair_collapse(cells)
    outp = os.path.join(OUT, f"forcesize_{camp}.txt")
    with open(outp, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {outp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
