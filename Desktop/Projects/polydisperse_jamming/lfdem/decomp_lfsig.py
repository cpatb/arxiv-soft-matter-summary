"""decomp_lfsig.py -- POST HOC (not in PREREG_LFSIG): the contact-pair closure of
derivations/sigma_star_geometry_2026-09-24.md tested on the LUBRICATED LF_DEM cells.

    python lfdem/decomp_lfsig.py  ->  out/decomp_lfsig.txt

Per cell (out/lfsig/<label>/: conf.dat radii, int_*.dat snapshots), for snapshots with strain in [1.5, 4]:
  contact pairs = contact state >= 1 (int_ col 3: 1 frictionless, 2 stuck, 3 sliding)
  Z_c = 2 N_contact / N;  <a_i a_j>_c;  closure C predicts <a_i a_j>_c = <a^3>/<a> (any PSD);
  closure A predicts (m3 m1 + m2^2)/(m2 + m1^2).
Reports <aa>_c / (m3/m1) and <aa>_c / A, per cell.
"""
from __future__ import annotations

import glob
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "out")
W0, W1 = 1.5, 4.0
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def read_radii(path: str) -> np.ndarray:
    a = np.loadtxt(path, comments="#", ndmin=2)
    return a[:, 3]


def snapshots(path: str):
    """Yield (strain, array of rows) for each block of an int_ file."""
    strain, rows = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("# cumulated strain"):
                if strain is not None and rows:
                    yield strain, np.array(rows)
                strain, rows = float(line.split(":")[1]), []
            elif line.startswith("#") or not line.strip():
                continue
            elif strain is not None:
                rows.append([float(x) for x in line.split()[:3]])
    if strain is not None and rows:
        yield strain, np.array(rows)


def main() -> int:
    say("decomp_lfsig -- POST HOC closure-C test on LF_DEM contacts (not registered)")
    say("  label           strain_snaps  Z_c     <aa>_c   <aa>/(m3/m1)  <aa>/A")
    for d in sorted(glob.glob(os.path.join(OUT, "lfsig", "c*"))):
        ints = glob.glob(os.path.join(d, "int_*.dat"))
        conf = os.path.join(d, "conf.dat")
        if not ints or not os.path.exists(conf):
            continue
        r = read_radii(conf)
        m1, m2, m3 = r.mean(), (r**2).mean(), (r**3).mean()
        A = (m3 * m1 + m2**2) / (m2 + m1**2)
        zs, aas, n = [], [], 0
        for s, arr in snapshots(ints[0]):
            if not (W0 <= s <= W1):
                continue
            c = arr[arr[:, 2] >= 1]
            i, j = c[:, 0].astype(int), c[:, 1].astype(int)
            zs.append(2 * len(c) / len(r))
            aas.append(float(np.mean(r[i] * r[j])) if len(c) else np.nan)
            n += 1
        if n == 0:
            continue
        aa = float(np.nanmean(aas))
        lab = os.path.basename(d)
        say(f"  {lab:15s} {n:3d}           {np.mean(zs):.3f}   {aa:.4f}   {aa / (m3 / m1):.4f}        {aa / A:.4f}")
    with open(os.path.join(OUT, "decomp_lfsig.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
