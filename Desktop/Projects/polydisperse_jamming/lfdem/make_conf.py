"""Write an LF_DEM text configuration for a lognormal PSD at a target volume fraction.

LF_DEM (Seto/Mari/Morris, github.com/ryseto/LF_DEM, commit 00fa7cc0) reads a text
configuration of the form

    # np1 np2 vf lx ly lz vf1 vf2 disp
    # <values>
    x y z radius
    ...

(`Configuration.h::readTxtBaseConfiguration`; only lx, ly, lz, vf are read from the
header, every particle carries its own radius). LF_DEM's own generator is bidisperse
only, so polydisperse packs are made here.

Radii come from THE project's sampler, `psd.LogNormalPSD(d10=2, cv).radii_for_yade`,
with the project's default truncation (0.2, 5.0) x d10 and number-mean radius 1 --
i.e. the same draw `packing.py` makes, in LF_DEM length units (radius 1 = the
number-mean radius). Overlaps from random insertion are removed by minimising the
harmonic overlap energy with scipy's L-BFGS (a library optimiser) on a periodic
cKDTree neighbour list; the configuration is accepted only if the largest remaining
overlap is below `tol` (relative to the pair's mean radius), which is checked, not
assumed.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from psd import LogNormalPSD  # noqa: E402

TRUNCATE = (0.2, 5.0)  # packing.py default, in units of d10


def radii(n: int, cv: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    r = LogNormalPSD(d10=2.0, cv=cv).radii_for_yade(n, rng, truncate=TRUNCATE)
    assert r.shape == (n,) and np.all(r > 0)
    return r


def _energy(x: np.ndarray, r: np.ndarray, L: float, pairs: np.ndarray):
    p = x.reshape(-1, 3)
    i, j = pairs[:, 0], pairs[:, 1]
    d = p[j] - p[i]
    d -= L * np.round(d / L)
    dist = np.linalg.norm(d, axis=1)
    ov = r[i] + r[j] - dist
    m = ov > 0
    e = 0.5 * np.sum(ov[m] ** 2)
    g = np.zeros_like(p)
    if np.any(m):
        u = d[m] / dist[m, None]
        f = ov[m, None] * u  # dE/dp_j = -ov*u ; dE/dp_i = +ov*u
        np.add.at(g, i[m], f)
        np.add.at(g, j[m], -f)
    return e, g.ravel()


def relax(r: np.ndarray, phi: float, seed: int, tol: float = 1e-6):
    n = r.size
    L = float((np.sum(4.0 / 3.0 * np.pi * r**3) / phi) ** (1.0 / 3.0))
    rng = np.random.default_rng(seed + 7919)
    p = rng.uniform(0, L, (n, 3))
    cut = 2.0 * r.max()
    for outer in range(200):
        p %= L
        tree = cKDTree(p, boxsize=L)
        pairs = tree.query_pairs(cut + 0.5, output_type="ndarray")
        res = minimize(_energy, p.ravel(), args=(r, L, pairs), jac=True,
                       method="L-BFGS-B", options={"maxiter": 5000, "gtol": 1e-12, "ftol": 1e-16})
        p = res.x.reshape(-1, 3) % L
        tree = cKDTree(p, boxsize=L)
        pr = tree.query_pairs(cut, output_type="ndarray")
        d = p[pr[:, 1]] - p[pr[:, 0]]
        d -= L * np.round(d / L)
        ov = (r[pr[:, 0]] + r[pr[:, 1]] - np.linalg.norm(d, axis=1)) / (0.5 * (r[pr[:, 0]] + r[pr[:, 1]]))
        worst = float(ov.max()) if ov.size else 0.0
        if worst < tol:
            return p, L, worst, outer + 1
    raise RuntimeError(f"relax failed: worst relative overlap {worst:.3e} after {outer+1} rounds")


def write(path: Path, p: np.ndarray, r: np.ndarray, L: float, phi_real: float) -> None:
    with open(path, "w", newline="\n") as fh:
        fh.write("# np1 np2 vf lx ly lz vf1 vf2 disp\n")
        fh.write(f"# {r.size} 0 {phi_real:.10g} {L:.10g} {L:.10g} {L:.10g} {phi_real:.10g} 0 0\n")
        for (x, y, z), a in zip(p, r):
            fh.write(f"{x:.10g} {y:.10g} {z:.10g} {a:.10g}\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--cv", type=float, required=True)
    ap.add_argument("--phi", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    r = radii(a.n, a.cv, a.seed)
    p, L, worst, rounds = relax(r, a.phi, a.seed)
    phi_real = float(np.sum(4.0 / 3.0 * np.pi * r**3) / L**3)
    assert abs(phi_real - a.phi) < 1e-9, (phi_real, a.phi)
    write(a.out, p, r, L, phi_real)
    d = 2 * r
    meta = {"n": a.n, "cv_req": a.cv, "phi": phi_real, "seed": a.seed, "L": L,
            "worst_rel_overlap": worst, "relax_rounds": rounds,
            "cv_real": float(np.std(d) / d.mean()), "a_mean": float(r.mean()),
            "a32": float(np.mean(r**3) / np.mean(r**2)), "a_max": float(r.max()),
            "L_over_dmax": float(L / d.max())}
    a.out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
