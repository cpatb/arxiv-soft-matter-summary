"""decomp_lfsig.py -- POST HOC (not in PREREG_LFSIG): the geometry of derivations/sigma_star_geometry and
rattler_closure applied to the LUBRICATED LF_DEM cells.

    python lfdem/decomp_lfsig.py  ->  out/decomp_lfsig.txt

Per cell (out/lfsig/<label>/: conf.dat radii, int_*.dat snapshots), snapshots with strain in [1.5, 4]:
  contacts = int_ col 3 state >= 1 (1 frictionless, 2 stuck, 3 sliding), i.e. surfaces touching (h < 0).
  Z_c = 2 N_contact / N;  <a_i a_j>_c against closure C (m3/m1), uncorrelated a^2b^2 ((m3/m2)^2) and the
  configuration model with the measured degrees ((sum k a / sum k)^2).
  Force-quantile sigma* in the actshape convention: at a touching contact the normal load is the contact
  spring force (col 12) plus the repulsion at contact (col 16, = F_R R* in LF_DEM units), so
  g = F_n/(R* sigma_int) with sigma_int = sigma_input/(6 pi) (SimulationInit.cpp:289) and
  sigma*/kappa = 1/Q_{1-1/e}(g) -- the stress at which 1/e of the TOUCHING contacts exceed F_R R*... note
  that in LF_DEM every touching contact already exceeds the repulsion by construction, so this is a
  formal transfer of the Yade estimator, not the rheological onset; it is printed to compare RATIOS.
"""
from __future__ import annotations

import glob
import math
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "out")
W0, W1 = 1.5, 4.0
INV_E = math.exp(-1.0)
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def snapshots(path: str):
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
                rows.append([float(x) for x in line.split()[:17]])
    if strain is not None and rows:
        yield strain, np.array(rows)


def main() -> int:
    say("decomp_lfsig -- POST HOC geometry on LF_DEM contacts (not registered; strain [1.5, 4])")
    say("  label           snaps  Z_c    <aa>_c  /C      /a2b2   /config  r_assort  s*/k(quantile)  n_contacts")
    res = {}
    for d in [x for x in sorted(glob.glob(os.path.join(OUT, "lfsig", "c*"))) if ".killed_" not in x]:
        ints = glob.glob(os.path.join(d, "int_*.dat"))
        conf = os.path.join(d, "conf.dat")
        if not ints or not os.path.exists(conf):
            continue
        r = np.loadtxt(conf, comments="#", ndmin=2)[:, 3]
        lab = os.path.basename(d)
        sig_in = float(re.search(r"_s([\d.]+)$", lab).group(1))
        sig = sig_in / (6 * math.pi)
        m1, m2, m3 = r.mean(), (r**2).mean(), (r**3).mean()
        Zs, aas, confs, rs, gs, n = [], [], [], [], [], 0
        for s, arr in snapshots(ints[0]):
            if not (W0 <= s <= W1):
                continue
            c = arr[arr[:, 2] >= 1]
            if len(c) == 0:
                continue
            i, j = c[:, 0].astype(int), c[:, 1].astype(int)
            k = np.bincount(np.concatenate([i, j]), minlength=len(r)).astype(float)
            Zs.append(k.mean())
            aas.append(float(np.mean(r[i] * r[j])))
            confs.append(float((np.sum(k * r) / np.sum(k)) ** 2))
            x, y = np.concatenate([r[i], r[j]]), np.concatenate([r[j], r[i]])
            rs.append(float(np.corrcoef(x, y)[0, 1]))
            Rs = r[i] * r[j] / (r[i] + r[j])
            gs.append((c[:, 11] + c[:, 15]) / (Rs * sig))
            n += 1
        if n == 0:
            continue
        g = np.concatenate(gs)
        sstar = 1.0 / float(np.quantile(g, 1 - INV_E))
        aa = float(np.mean(aas))
        res[lab] = dict(Z=np.mean(Zs), sstar=sstar)
        say(f"  {lab:15s} {n:4d}  {np.mean(Zs):.3f}  {aa:.4f}  {aa / (m3 / m1):.4f}  {aa / (m3 / m2) ** 2:.4f}  "
            f"{aa / np.mean(confs):.4f}   {np.mean(rs):+.4f}   {sstar:.5g}        {g.size}")
    say("\n  CV 0.5 / CV 0.2 at matched (x_m, sigma):  Z_c ratio   sigma*/kappa(quantile) ratio")
    for lab in sorted(res):
        if lab.startswith("c50"):
            other = "c20" + lab[3:]
            if other in res:
                say(f"    {lab[4:]:12s}  {res[lab]['Z'] / res[other]['Z']:.4f}     {res[lab]['sstar'] / res[other]['sstar']:.4f}")
    with open(os.path.join(OUT, "decomp_lfsig.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
