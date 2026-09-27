"""facet_census.py -- the registered F5 null: face-alignment of a contact normal under UNIFORM-AREA surface sampling.

For the equal-axis superellipsoid |x|^(2/eps) + |y|^(2/eps) + |z|^(2/eps) = 1 (eps 1 sphere, eps -> 0 cube), the outward
normal at a surface point is grad F / |grad F| with F the left-hand side. m = max_i |n_i| is 1 at a face centre and
1/sqrt(3) on a corner diagonal. Area-uniform sampling: draw points uniformly in a thin shell |G - 1| < DELTA of the
degree-1 homogeneous G = F^(eps/2); the shell's normal thickness is DELTA/|grad G|, so a volume-uniform draw has
density proportional to 1/|grad G| per unit area, and weighting each point by |grad G| makes it area-uniform.

Anchors (checked in main): eps = 1 must give the uniform-sphere law P(max|n_i| > c) (computed independently by drawing
Gaussian unit vectors); the area-weighted mean of n must vanish; and P(m > c) must rise as eps falls (flatter faces).
"""
from __future__ import annotations

import numpy as np

DELTA = 0.01
C_FACE = 0.95  # registered F5 threshold (PREREG_FACET)


def area_null(eps: float, c: float = C_FACE, n_draw: int = 4_000_000, seed: int = 0) -> dict:
    rng = np.random.default_rng(seed)
    # the box must contain the whole shell: G < 1 + DELTA reaches |x_i| = 1 + DELTA at the face centres, and a box of
    # half-width 1 truncates exactly the face-centred region the statistic measures (first version: 0.1443 vs 0.1501)
    x = rng.uniform(-(1.0 + 2 * DELTA), 1.0 + 2 * DELTA, size=(n_draw, 3))
    u = np.abs(x)
    F = np.sum(u ** (2.0 / eps), axis=1)
    G = F ** (eps / 2.0)
    keep = np.abs(G - 1.0) < DELTA
    x, u, F = x[keep], u[keep], F[keep]
    gF = (2.0 / eps) * u ** (2.0 / eps - 1.0) * np.sign(x)
    gG = (eps / 2.0) * F[:, None] ** (eps / 2.0 - 1.0) * gF
    w = np.linalg.norm(gG, axis=1)
    n = gF / np.linalg.norm(gF, axis=1)[:, None]
    m = np.max(np.abs(n), axis=1)
    p = float(np.sum(w * (m > c)) / np.sum(w))
    mean_n = np.sum(w[:, None] * n, axis=0) / np.sum(w)
    return dict(eps=eps, p_face=p, f_ff0=p * p, n_shell=int(keep.sum()), mean_n_abs=float(np.max(np.abs(mean_n))))


def sphere_law(c: float = C_FACE, n_draw: int = 2_000_000, seed: int = 1) -> float:
    v = np.random.default_rng(seed).normal(size=(n_draw, 3))
    v /= np.linalg.norm(v, axis=1)[:, None]
    return float(np.mean(np.max(np.abs(v), axis=1) > c))


def main() -> int:
    ref = sphere_law()
    rows = [area_null(e) for e in (1.0, 0.9, 0.8, 0.7, 0.5)]
    print(f"uniform-sphere law P(m > {C_FACE}) = {ref:.4f} (Gaussian draw, independent anchor)")
    for r in rows:
        print(f"eps {r['eps']:.1f}: P_area(m > {C_FACE}) = {r['p_face']:.4f}  f_ff0 = {r['f_ff0']:.5f}  "
              f"(shell points {r['n_shell']}, |<n>| {r['mean_n_abs']:.1e})")
    assert abs(rows[0]["p_face"] - ref) < 0.01, "eps 1 area null must equal the uniform-sphere law"
    assert all(r["mean_n_abs"] < 5e-3 for r in rows), "area-weighted mean normal must vanish"
    ps = [r["p_face"] for r in rows]
    assert all(b > a for a, b in zip(ps, ps[1:])), "P(m > c) must rise as the faces flatten"
    print("anchors: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
