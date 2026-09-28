"""q6_order.py -- Amendment 4 ordering check: global bond-orientational order Q6 (Steinhardt, Nelson & Ronchetti 1983)
from LF_DEM interaction dumps.

    python lfdem/q6_order.py lfsig2m lfsig2 ...   ->  out/q6_order.txt

Neighbours: every pair in the dump with dimensionless gap s - 2 < 0.4, i.e. centre distance < 1.2 x contact distance
(s = 2 r/(a_i + a_j); the dump lists every pair inside lub_max_gap 0.5, so this shell is complete). Bond directions are
the dump's normal vectors (cols 4-6). Q6 = sqrt(4 pi/13 sum_m |<Y_6m(bond)>|^2), averaged over bonds (bond-order is
invariant to the bond's sign because l = 6 is even). Read from the LAST snapshot of each cell. Registered bar: Q6 > 0.3
= ORDERED (random isotropic: ~ 1/sqrt(N_bonds); fcc 0.575). Library spherical harmonics (scipy.special).

Anchors, printed every run (outside the code under test): isotropic random unit vectors (Q6 ~ 1/sqrt(n)), and the 12 fcc
nearest-neighbour directions (Q6 = 0.5745, the textbook value).
"""
from __future__ import annotations

import glob
import os
import re
import sys

import numpy as np
from scipy import special

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mech_lfsig import snapshots  # noqa: E402

OUT = os.path.join(os.path.dirname(HERE), "out")
GAP_MAX = 0.4
BAR = 0.3


def ylm(l: int, m: int, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Y_lm(polar theta, azimuth phi) through whichever scipy API is present."""
    if hasattr(special, "sph_harm_y"):
        return special.sph_harm_y(l, m, theta, phi)
    return special.sph_harm(m, l, phi, theta)   # legacy argument order


def q6(n: np.ndarray) -> float:
    n = n / np.linalg.norm(n, axis=1, keepdims=True)
    theta = np.arccos(np.clip(n[:, 2], -1, 1))
    phi = np.arctan2(n[:, 1], n[:, 0])
    s = sum(abs(ylm(6, m, theta, phi).mean()) ** 2 for m in range(-6, 7))
    return float(np.sqrt(4 * np.pi / 13 * s))


def anchors() -> list[str]:
    rng = np.random.default_rng(20260927)
    v = rng.normal(size=(2000, 3))
    fcc = np.array([(a, b, 0) for a in (1, -1) for b in (1, -1)] + [(a, 0, b) for a in (1, -1) for b in (1, -1)]
                   + [(0, a, b) for a in (1, -1) for b in (1, -1)], dtype=float)
    q_r, q_f = q6(v), q6(fcc)
    assert abs(q_f - 0.5745) < 1e-3, q_f
    return [f"  anchors: random 2000 bonds Q6 {q_r:.4f} (expect ~{1 / np.sqrt(2000):.3f}); fcc 12 bonds Q6 {q_f:.4f} (textbook 0.5745)"]


def main() -> int:
    sets = sys.argv[1:] or ["lfsig2m"]
    lines = ["q6_order -- Amendment 4 ordering check (global Q6 over bonds with gap < 0.4 at the last snapshot; ORDERED if > 0.3)"]
    lines += anchors()
    for st in sets:
        for d in sorted(glob.glob(os.path.join(OUT, st, "c*"))):
            if ".killed_" in d:
                continue
            ints = glob.glob(os.path.join(d, "int_*.dat"))
            if not ints:
                continue
            last = None
            for strain, _r, s in snapshots(ints[0]):
                last = (strain, s)
            if last is None:
                continue
            strain, s = last
            sel = s[:, 6] < GAP_MAX
            q = q6(s[sel, 3:6])
            nb = int(sel.sum())
            lab = dict(re.findall(r"(\w+)=(\S+)", open(os.path.join(d, "META"), encoding="utf-8").read()))["LABEL"]
            lines.append(f"  {st:8s} {lab:20s} strain {strain:4.1f} bonds {nb:5d} Q6 {q:.4f} (random ~{1 / np.sqrt(nb):.3f}) "
                         f"{'ORDERED' if q > BAR else 'random'}")
    txt = "\n".join(lines) + "\n"
    print(txt)
    with open(os.path.join(OUT, "q6_order.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
