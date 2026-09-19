"""geom_common_2026-09-19.py -- shared loaders and closed-form PSD volume-share functions.

Imported by every geom_*_2026-09-19.py in this directory. READ ONLY on out/.

Two things live here:

1. `load_fits()` / `load_shards()`. Deliberately NOT `sweep_descriptive.load_shards`,
   for one reason: that function globs `out/phase2*.jsonl`, which on 2026-09-19 also
   matches `out/phase2_wcbidi*.jsonl`, a campaign that is MID-FLIGHT. A half-written
   shard would silently join every mean. This loader excludes it by name and prints the
   exclusion count. The dedup rule is otherwise identical to the shipped one (drop
   failed, key-dedup keeping the highest MAX_ITER, then collapse one physical cell per
   (family, cv, seed, p_hat, mu_p) -- the 2026-09-17 `_it{cap}` defect).

2. `Psi(family, cv, q)` -- the VOLUME share of the smallest number-fraction `q` of a
   PSD, in closed form for each of the three families. This is the object the rattler
   decomposition needs and it is not a standard moment; see
   geom_rattler_decomposition_2026-09-19.md.
"""

from __future__ import annotations

import glob
import json
import math
import os

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "out")
N_MATCH = 1500
T_CRIT_2DOF = 4.302652729911275   # scipy.stats.t.ppf(0.975, 2); 3 seeds -> 2 dof


# ----------------------------------------------------------------- loaders
def load_fits(verbose=True) -> pd.DataFrame:
    """Per-cell extrapolated phi_J_lin. Same rule as sweep_descriptive.load_fits."""
    d = pd.read_csv(os.path.join(OUT, "phase2_fits.csv"))
    n0 = len(d)
    d = d[(d.n_requested == N_MATCH) & np.isfinite(d.phi_J_lin)]
    n1 = len(d)
    d = d.loc[sorted(d.groupby(["family", "cv_nom", "mu_p", "seed"])["max_iter"].idxmax())]
    if verbose:
        print(f"[load_fits] N rows in phase2_fits.csv: {n0}")
        print(f"[load_fits] N after n_requested=={N_MATCH} and finite phi_J_lin: {n1}")
        print(f"[load_fits] N after highest-MAX_ITER per (family,cv,mu,seed): {len(d)}")
    return d


def load_shards(verbose=True) -> pd.DataFrame:
    """Per-rung shard records at N=1500, converged, one record per physical cell."""
    allf = sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl")))
    files = [f for f in allf if "wcbidi" not in os.path.basename(f)]
    rows = []
    for f in files:
        for ln in open(f, encoding="utf-8"):
            ln = ln.strip()
            if not ln:
                continue
            try:
                rows.append(json.loads(ln))
            except json.JSONDecodeError:
                continue
    by = {}
    for r in rows:
        if r.get("failed"):
            continue
        k = r.get("key")
        if k in by and r.get("max_iter", 0) <= by[k].get("max_iter", 0):
            continue
        by[k] = r
    sel = [r for r in by.values()
           if r.get("n_requested") == N_MATCH and r.get("converged")]
    cap = {}
    for r in sorted(sel, key=lambda x: x.get("key", "")):
        pk = ((r.get("family") or "lognormal"), r["cv"], r["seed"],
              r["p_hat_target"], r.get("mu_p"))
        if pk not in cap or r.get("max_iter", 0) > cap[pk].get("max_iter", 0):
            cap[pk] = r
    d = pd.DataFrame([{
        "family": (r.get("family") or "lognormal"), "cv": r["cv"], "seed": r["seed"],
        "mu_p": r.get("mu_p"), "p": r["p_hat_target"], "phi": r["phi"],
        "rf": r["rattler_frac"], "z_nr": r["z_nr"], "z_all": r["z_all"],
        "n_contacts": r["n_contacts"], "n_placed": r["n_placed"],
        "cv_real": r.get("cv_realised"), "d10": r.get("d10_realised"),
        "d32": r.get("d32_realised"), "d43": r.get("d43_realised"),
        "max_iter": r.get("max_iter"),
    } for r in cap.values()])
    if verbose:
        print(f"[load_shards] files matching out/phase2*.jsonl: {len(allf)}")
        print(f"[load_shards] EXCLUDED as mid-flight (name contains 'wcbidi'): "
              f"{len(allf) - len(files)}")
        print(f"[load_shards] N records read: {len(rows)}")
        print(f"[load_shards] N after dropping failed + key-dedup (highest cap): {len(by)}")
        print(f"[load_shards] N after n_requested=={N_MATCH} and converged: {len(sel)}")
        print(f"[load_shards] N after collapse to one physical cell: {len(d)}")
    return d


# ------------------------------------------------- closed-form PSD volume share
def _lognormal_sigma(cv: float) -> float:
    return math.sqrt(math.log(1.0 + cv * cv))


def Psi(family: str, cv: float, q: float) -> float:
    """Volume share of the SMALLEST number-fraction `q` of the PSD.

        Psi(q) = int_0^{F^{-1}(q)} d^3 f(d) dd / <d^3>

    Closed form per family; `q` in [0, 1]. Psi(0)=0, Psi(1)=1, and Psi(q) <= q with
    equality iff the PSD is a point mass (CV = 0).
    """
    q = float(np.clip(q, 0.0, 1.0))
    if q <= 0.0:
        return 0.0
    if q >= 1.0:
        return 1.0
    if cv == 0.0:
        return q                      # point mass: volume share == number share
    if family == "lognormal":
        # d = exp(mu + sigma Z); E[d^3 1(Z<z)] = E[d^3] Phi(z - 3 sigma).
        s = _lognormal_sigma(cv)
        return float(norm.cdf(norm.ppf(q) - 3.0 * s))
    if family == "uniform":
        # d ~ U[1-w, 1+w], w = cv*sqrt(3). F^{-1}(q) = 1-w+2wq.
        w = cv * math.sqrt(3.0)
        lo, hi = 1.0 - w, 1.0 + w
        t = lo + 2.0 * w * q
        return float((t**4 - lo**4) / (hi**4 - lo**4))
    if family == "bidisperse":
        # atoms at 1-cv and 1+cv, 50:50 BY NUMBER.
        a, b = 1.0 - cv, 1.0 + cv
        ws = a**3 / (a**3 + b**3)     # small grains' share of the solid volume
        if q <= 0.5:
            return float(2.0 * q * ws)
        return float(ws + (2.0 * q - 1.0) * (1.0 - ws))
    raise ValueError(f"unknown family {family!r}")


def psi_bracket(family: str, cv: float, q: float) -> tuple[float, float]:
    """Exact bounds on the rattler VOLUME fraction given the rattler NUMBER fraction q.

    Rearrangement: the volume of any subset of number-fraction q is minimised by taking
    the q smallest grains and maximised by taking the q largest.
    """
    return Psi(family, cv, q), 1.0 - Psi(family, cv, 1.0 - q)


def tstat(x) -> tuple[float, float, float, int]:
    """mean, se, t, n for a 1-D sample. t is mean/se; compare against T_CRIT_2DOF at n=3."""
    x = np.asarray(x, float)
    n = len(x)
    m = float(x.mean())
    se = float(x.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")
    return m, se, (m / se if se > 0 else float("nan")), n


if __name__ == "__main__":
    print("self-test of Psi: Psi(q) <= q, Psi(0)=0, Psi(1)=1, monotone")
    for fam, cv in [("lognormal", 0.5), ("uniform", 0.5), ("bidisperse", 0.5)]:
        qs = np.linspace(0, 1, 11)
        v = [Psi(fam, cv, q) for q in qs]
        print(f"  {fam:11s} cv={cv}: " + " ".join(f"{x:.4f}" for x in v))
        assert all(v[i] <= v[i + 1] + 1e-12 for i in range(len(v) - 1))
        assert all(v[i] <= qs[i] + 1e-12 for i in range(len(v)))
    for fam in ("lognormal", "uniform", "bidisperse"):
        assert abs(Psi(fam, 0.0, 0.37) - 0.37) < 1e-12
    print("  OK")
