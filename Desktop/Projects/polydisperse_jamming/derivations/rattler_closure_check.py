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
from psd import BidispersePSD, LogNormalPSD, UniformPSD  # noqa: E402

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


_fam_cache: dict = {}


def radii_family(family: str, cv: float) -> np.ndarray:
    """Radii (number mean 1/2 of d10 = 1) for the three phase-2 families, via psd.py samplers."""
    key = (family, cv)
    if key not in _fam_cache:
        rng = np.random.default_rng(20260925)
        if family == "lognormal":
            d = LogNormalPSD(d10=1.0, cv=cv).sample(200_000, rng)
        elif family == "uniform":
            d = UniformPSD(d10=1.0, cv=cv).sample(200_000, rng)
        elif family == "bidisperse":
            d = BidispersePSD(d10=1.0, cv=cv).sample(200_000, rng)
        else:
            raise ValueError(family)
        _fam_cache[key] = d / 2.0
    return _fam_cache[key]


def isotropic_block() -> None:
    """OUT OF SAMPLE: isotropically compressed phase-2 packs (a different protocol from the sheared
    actshape packs the closure was found on), three PSD families, rattler = iterated k < 2 (measure.py).
    Closure (i) only (from the record's own z_all). Converged, alpha 1, no multiplets, CV > 0."""
    say("\n=== OUT OF SAMPLE: isotropic phase-2 records, closure (i) from z_all ===")
    recs = []
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase2*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if ("error" in r or not r.get("converged") or r.get("family") not in ("lognormal", "uniform", "bidisperse")
                    or r.get("alpha", 1.0) != 1.0 or r.get("doublet_w", 0) or r.get("triplet_w", 0)
                    or "z_all" not in r or "rattler_frac" not in r):
                continue
            recs.append(r)
    # dedupe by key (a key can repeat across shards; keep the last, as the loaders do)
    recs = list({r["key"]: r for r in recs}.values())
    say(f"  records: {len(recs)}")
    say("  family      mu   cv    n    fr_meas(mean)  (i) C: mean err  rms     null: mean err  rms")
    groups: dict = {}
    for r in recs:
        groups.setdefault((r["family"], r["mu_p"], round(r["cv"], 3)), []).append(r)
    allC, all0 = [], []
    for (fam, mu, cv), rs in sorted(groups.items()):
        a = radii_family(fam, cv)
        eC = [fr_from_Z(a, r["z_all"]) - r["rattler_frac"] for r in rs]
        e0 = [fr_from_Z(a, r["z_all"], True) - r["rattler_frac"] for r in rs]
        allC += eC
        all0 += e0
        eC, e0 = np.array(eC), np.array(e0)
        say(f"  {fam:10s} {mu:.1f} {cv:.2f} {len(rs):4d}  {np.mean([r['rattler_frac'] for r in rs]):.4f}"
            f"         {eC.mean():+.4f}  {np.sqrt(np.mean(eC**2)):.4f}   {e0.mean():+.4f}  {np.sqrt(np.mean(e0**2)):.4f}")
    allC, all0 = np.array(allC), np.array(all0)
    say(f"  ALL: closure C rms {np.sqrt(np.mean(allC**2)):.4f} (mean {allC.mean():+.4f});  size-blind null rms "
        f"{np.sqrt(np.mean(all0**2)):.4f} (mean {all0.mean():+.4f})")
    backbone_block(groups)


def psi_from_Z(a: np.ndarray, Z: float) -> float:
    """Model rattler VOLUME share: <a^3 P(k<2 | lam(a))>/<a^3>."""
    pr = poisson.cdf(K_MIN - 1, Z * shape(a, False))
    return float(np.sum(a**3 * pr) / np.sum(a**3))


def backbone_block(groups: dict) -> None:
    """geom_rattler_decomposition_2026-09-19 left backbone invariance UNDETERMINED because no shard records
    the rattler VOLUME share psi. The count-validated model supplies it: phi_B = phi (1 - psi_model).
    Printed per (family, mu, CV) at the p_hat 3e-5 rung (jam_map.RUNG, the window convention), seed means.
    This is a MODEL-derived psi, not a measurement; its validation is the count test above."""
    say("\n=== backbone solid fraction phi_B = phi (1 - psi_model), p_hat 3e-5 rung, seed means ===")
    say("  family      mu   cv    n   phi      f_r(meas)  psi_model  phi_B     d(phi_B) vs CV 0")
    base = {}
    for (fam, mu, cv), rs in sorted(groups.items()):
        rs = [r for r in rs if abs(r["p_hat_target"] - 3e-5) < 1e-9]
        if not rs:
            continue
        a = radii_family(fam, cv)
        phi = np.mean([r["phi"] for r in rs])
        psi = np.mean([psi_from_Z(a, r["z_all"]) for r in rs])
        pb = np.mean([r["phi"] * (1 - psi_from_Z(a, r["z_all"])) for r in rs])
        if cv == 0.0:
            base[(fam, mu)] = pb
        d = pb - base.get((fam, mu), np.nan)
        say(f"  {fam:10s} {mu:.1f} {cv:.2f} {len(rs):4d}  {phi:.4f}   {np.mean([r['rattler_frac'] for r in rs]):.4f}     "
            f"{psi:.4f}     {pb:.4f}    {d:+.4f}")


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
    isotropic_block()
    with open(os.path.join(ROOT, "out", "rattler_closure.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
