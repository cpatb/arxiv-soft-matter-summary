"""Write an LF_DEM configuration for a BINARY suspension (PREREG_LFSIG2 Amendment 7).

Radii are deterministic: N_s grains of radius 1 (the length unit, so sigma_0 = F_R0/a_s^2 here) and N_l of radius
`ratio`, with the small species carrying the fraction `x_small` of the solid volume
(n_s = x_s r^3/(1 - x_s + x_s r^3), `psd.BinaryPSD`'s relation), rounded to integers. Positions are seeded and relaxed to
non-overlap by `make_conf.relax` (the same routine and tolerance the lognormal campaigns use), written by
`make_conf.write`. No random radius draw, so seeds differ only in positions.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from make_conf import relax, write  # noqa: E402


def radii(n: int, ratio: float, x_small: float) -> np.ndarray:
    ns_frac = x_small * ratio**3 / (1 - x_small + x_small * ratio**3)
    ns = int(round(n * ns_frac))
    r = np.concatenate([np.ones(ns), np.full(n - ns, ratio)])
    assert r.size == n and 0 < ns < n
    return r


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--ratio", type=float, required=True)
    ap.add_argument("--xsmall", type=float, required=True)
    ap.add_argument("--phi", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()
    r = radii(a.n, a.ratio, a.xsmall)
    p, L, worst, rounds = relax(r, a.phi, a.seed)
    phi_real = float(np.sum(4.0 / 3.0 * np.pi * r**3) / L**3)
    assert abs(phi_real - a.phi) < 1e-9, (phi_real, a.phi)
    write(a.out, p, r, L, phi_real)
    ns = int((r == 1.0).sum())
    xs_real = float(ns * 1.0 / np.sum(r**3))
    meta = {"n": a.n, "ratio": a.ratio, "x_small_req": a.xsmall, "x_small_real": xs_real, "n_small": ns,
            "phi": phi_real, "seed": a.seed, "L": L, "worst_rel_overlap": worst, "relax_rounds": rounds,
            "cv_real": float(np.std(2 * r) / np.mean(2 * r)), "a_mean": float(r.mean()),
            "a32": float(np.mean(r**3) / np.mean(r**2)), "a_max": float(r.max()), "L_over_dmax": float(L / (2 * r.max()))}
    a.out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
