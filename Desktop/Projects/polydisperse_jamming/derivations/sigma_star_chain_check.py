"""sigma_star_chain_check.py -- the closed chain for sigma*(CV)/sigma*(0) at fixed number-mean size.

    python derivations/sigma_star_chain_check.py  ->  out/sigma_star_chain.txt

Inputs per sheared actshape record: the backbone count Z_nr (z_nr_ss), phi (phi_ss) and the PSD (nominal
CV, project sampler). Steps, each from a refereed derivation:
  1. degree law k(a) ~ Poisson(Z a^2/<a^2>)           (rattler_closure_2026-09-24, best OOS form)
     Z solved so that the mean degree over NON-rattlers (k >= 2) equals Z_nr.
  2. configuration model: <a_i a_j>_c = a_k^2, a_k = <k a>/<k> = <a^3>/<a^2> for k ~ a^2   (sigma_star §5b)
  3. sigma*/kappa = Z phi a_k^2 / (8 pi <a^3>) * C/q    (eq. 15; C/q ratio taken MEASURED from the dumps)
Printed: the predicted geometric ratio G_pred = [Z phi a_k^2/<a^3>]_CV / [...]_0 (at <a> = 1/2), the same
with C/q, and the measured R (seed-paired means from out/sigma_star_geometry.txt ratio table).
Z_nr is a single-snapshot record field (sigma_star §7), so this inherits its noise.
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys

import numpy as np
from scipy.optimize import brentq
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
from rattler_closure_check import radii  # noqa: E402

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def z_from_znr(a: np.ndarray, znr: float) -> float:
    sh = a * a / np.mean(a * a)

    def f(Z):
        lam = Z * sh
        return float(np.sum(lam - lam * np.exp(-lam)) / np.sum(poisson.sf(1, lam))) - znr

    return brentq(f, 1e-3, 50.0)


# mean-degree exponent measured on the isotropic per-grain dumps (out/degree_exponent_probe.txt: 1.55, 1.70)
P_MEAN = (1.55 + 1.70) / 2


def geo(a: np.ndarray, Z: float, phi: float, p: float = 2.0) -> float:
    m1, m3 = a.mean(), np.mean(a**3)
    ak = np.mean(a ** (p + 1)) / np.mean(a**p)
    s = 0.5 / m1                                   # rescale to <a> = 1/2 (d10 = 1)
    return Z * phi * (ak * s) ** 2 / (m3 * s**3)


def main() -> int:
    recs = []
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase9_actshape_*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if "error" not in r and r["steady"] and r["servo_held"] and r["quasistatic"]:
                recs.append(r)
    txt = open(os.path.join(ROOT, "out", "sigma_star_geometry.txt"), encoding="utf-8").read()
    meas = {}
    for m in re.finditer(r"^\s+(0\.[03]) (\de-0\d)\s+(0\.[25])\s+\d+\s+(\S+)\s+=\s+\S+\s+x\s+(\S+)\s+x\s+(\S+)", txt, re.M):
        mu, p, cv, R, C, iq = m.groups()
        meas[(float(mu), float(p), float(cv))] = (float(R), float(C) * float(iq))
    groups = {}
    for r in recs:
        groups.setdefault((r["mu_p"], r["p_hat_target"]), {}).setdefault(r["cv"], []).append(r)
    say("closed chain: Z_nr, phi, PSD -> Z (a^2 Poisson) -> configuration model -> sigma* ratio vs CV 0")
    say("  mu   p      CV   G_pred   x C/q(meas)   R_meas   pred/meas  | Z_pred/Z_meas(CV)")
    for (mu, p), byc in sorted(groups.items()):
        if 0.0 not in byc:
            continue
        g0 = np.mean([geo(radii(0.0), z_from_znr(radii(0.0), r["z_nr_ss"]), r["phi_ss"]) for r in byc[0.0]])
        g0p = np.mean([geo(radii(0.0), z_from_znr(radii(0.0), r["z_nr_ss"]), r["phi_ss"], P_MEAN) for r in byc[0.0]])
        for cv in (0.2, 0.5):
            if cv not in byc or (mu, p, cv) not in meas:
                continue
            a = radii(cv)
            zs = [(z_from_znr(a, r["z_nr_ss"]), r["z_ss"], r["phi_ss"]) for r in byc[cv]]
            g = np.mean([geo(a, z, ph) for z, _, ph in zs])
            R, cq = meas[(mu, p, cv)]
            G = g / g0
            Gp = np.mean([geo(a, z, ph, P_MEAN) for z, _, ph in zs]) / g0p
            say(f"  {mu:.1f} {p:.0e}  {cv:.1f}  {G:.4f}   {G * cq:.4f}        {R:.4f}   {G * cq / R:.4f}     | "
                f"{np.mean([z / zm for z, zm, _ in zs]):.3f}   drop captured (1-pred)/(1-R) {(1 - G * cq) / (1 - R):.2f}"
                f"   | a_k with p={P_MEAN:.3f}: pred {Gp * cq:.4f} pred/meas {Gp * cq / R:.4f}")
    with open(os.path.join(ROOT, "out", "sigma_star_chain.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
