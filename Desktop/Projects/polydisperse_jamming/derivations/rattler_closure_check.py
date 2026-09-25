"""rattler_closure_check.py -- checks for derivations/rattler_closure_2026-09-24.md.

    python derivations/rattler_closure_check.py  ->  out/rattler_closure.txt

Model (no fitted parameter):
  closure C (sigma_star_geometry eq. 12): pair weight w(a_i, a_j) ~ a_i a_j (a_i + a_j), so the mean
  contact number of a grain of radius a is  lam(a) = c * a (a m1 + m2)            (eq. 3)
  (normalised: <lam> = c * 2 m1 m2).
  per-grain contact counts k ~ Poisson(lam(a))                                     (assumption P)
  rattler (measure.py primary definition, NOT iterated here): k < 2                (eq. 4)
Two closures:
  (i)  c from the measured all-grain Z:     <lam> = Z                               -> f_r
  (ii) c from the measured BACKBONE Z_nr:   Z_nr = E[k; k>=2] / P(k>=2)             -> f_r, Z
Null: size-blind lam(a) = const.
Then sigma*(CV)/sigma*(0) from eq. 13 of sigma_star_geometry with Z predicted by (ii):
  R_pred = [Z phi / <a>]_CV / [Z phi / <a>]_0 x (C/q ratio measured)  vs  R_meas (actshape dumps).
Radii: psd.LogNormalPSD sampler with the project's truncation, at the record's NOMINAL cv (the
realised cv of each cell is printed; it differs by finite N).
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from psd import LogNormalPSD  # noqa: E402

K_MIN = 2           # measure.RATTLER_MIN_CONTACTS
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


_cache: dict[float, np.ndarray] = {}


def radii(cv: float) -> np.ndarray:
    if cv not in _cache:
        rng = np.random.default_rng(20260924)
        _cache[cv] = LogNormalPSD(d10=2.0, cv=cv).radii_for_yade(400_000, rng)
    return _cache[cv]


def shape(a: np.ndarray, sizeblind: bool) -> np.ndarray:
    if sizeblind:
        return np.ones_like(a)
    m1, m2 = a.mean(), (a * a).mean()
    s = a * (a * m1 + m2)
    return s / s.mean()


def fr_from_Z(a, Z, sizeblind=False):
    lam = Z * shape(a, sizeblind)
    return float(np.mean(poisson.cdf(K_MIN - 1, lam)))


def from_Znr(a, Znr, sizeblind=False):
    sh = shape(a, sizeblind)

    def backbone_mean(Z):
        lam = Z * sh
        p_ok = poisson.sf(K_MIN - 1, lam)                   # P(k >= 2)
        e_ok = lam - lam * poisson.pmf(0, lam)             # E[k; k>=2] = lam - 1*P(k=1) = lam - lam e^-lam
        return float(np.sum(e_ok) / np.sum(p_ok))

    Z = brentq(lambda z: backbone_mean(z) - Znr, 1e-3, 50.0)
    return Z, fr_from_Z(a, Z, sizeblind)


def main() -> int:
    rows = []
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase9_actshape_*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if "error" in r or not (r["steady"] and r["servo_held"] and r["quasistatic"]):
                continue
            rows.append(r)
    rows.sort(key=lambda r: (r["mu_p"], r["p_hat_target"], r["cv"], r["n_placed"], r["seed"]))
    say("rattler closure -- Poisson x closure C, rattler = k < 2 (non-iterated); records phase9_actshape_*")
    say(f"healthy records: {len(rows)}")
    say(" mu   p      cv  cvreal N     Z     Znr   | fr_meas | (i) C    null  | (ii) Z_pred  fr_C   fr_null")
    res = {"i": [], "i0": [], "ii": [], "ii0": [], "Z": []}
    for r in rows:
        a = radii(r["cv"])
        Z, Znr, fr = r["z_ss"], r["z_nr_ss"], r["rattler_frac_ss"]
        f1, f10 = fr_from_Z(a, Z), fr_from_Z(a, Z, True)
        Z2, f2 = from_Znr(a, Znr)
        _, f20 = from_Znr(a, Znr, True)
        res["i"].append(f1 - fr); res["i0"].append(f10 - fr)
        res["ii"].append(f2 - fr); res["ii0"].append(f20 - fr); res["Z"].append(Z2 / Z - 1)
        say(f" {r['mu_p']:.1f} {r['p_hat_target']:.0e} {r['cv']:.1f} {r['cv_realised']:.3f} {r['n_placed']:5d} "
            f"{Z:.3f} {Znr:.3f} | {fr:.3f}   | {f1:.3f}   {f10:.3f} | {Z2:.3f}       {f2:.3f}  {f20:.3f}")
    say()
    for k, name in (("i", "(i) closure C, from Z"), ("i0", "(i) size-blind null, from Z"),
                    ("ii", "(ii) closure C, from Z_nr"), ("ii0", "(ii) size-blind null, from Z_nr")):
        e = np.array(res[k])
        say(f"  {name:34s}: fr_pred - fr_meas  mean {e.mean():+.4f}  rms {np.sqrt(np.mean(e**2)):.4f}  "
            f"max|.| {np.abs(e).max():.4f}")
    e = np.array(res["Z"])
    say(f"  (ii) Z_pred/Z_meas - 1: mean {e.mean():+.4f}  rms {np.sqrt(np.mean(e**2)):.4f}  max|.| {np.abs(e).max():.4f}")
    # by CV
    for cv in (0.0, 0.2, 0.5):
        idx = [i for i, r in enumerate(rows) if r["cv"] == cv]
        for k in ("ii", "ii0"):
            e = np.array(res[k])[idx]
            say(f"    CV {cv}: {k:4s} mean {e.mean():+.4f}  rms {np.sqrt(np.mean(e**2)):.4f}  (n {len(idx)})")
    with open(os.path.join(ROOT, "out", "rattler_closure.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
