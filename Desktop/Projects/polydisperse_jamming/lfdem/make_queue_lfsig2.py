"""Write lfdem/queue_lfsig2.txt for PREREG_LFSIG2_2026-09-27.md.

Arms: CV {0.2, 0.35, 0.5} at x_m 0.90 (Yade mu_p 0.3, p_hat 3e-5 phi_m convention), seeds {1, 2, 3},
stress ladder {1, 2, 3, 4.5, 7, 10, 15, 30, 100} sigma_0, strain 4, mu 0.3. 81 cells, ordered so each
seed's crossing-region cells run first. phi_m at CV 0.35 is INTERPOLATED (quadratic through the three
measured phi_m; registered as sign-only). Prints every phi it writes and the interpolation check.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
XM = 0.90
# phi_m(mu_p 0.3, p_hat 3e-5), converged records, from out/predict_cvmap.txt (sd 0.00125 / 0.00044 / 0.00323)
PHIM = {0.0: 0.57535, 0.2: 0.58755, 0.5: 0.63667}
SIGMAS = (1, 2, 3, 4.5, 7, 10, 15, 30, 100)
ORDER = (4.5, 7, 3, 10, 100, 2, 15, 30, 1)   # decisive first; plateau (slowest) before the outer rungs
SEEDS = (1, 2, 3)
STRAIN = 4
MU = 0.3


def phim_interp(cv: float) -> float:
    cvs = np.array(sorted(PHIM)); ph = np.array([PHIM[c] for c in cvs])
    coef = np.polyfit(cvs, ph, 2)          # exact through three points
    assert np.allclose(np.polyval(coef, cvs), ph, atol=1e-12)
    return float(np.polyval(coef, cv))


def main() -> int:
    assert set(ORDER) == set(SIGMAS) and len(ORDER) == len(SIGMAS)
    phi = {0.2: XM * PHIM[0.2], 0.5: XM * PHIM[0.5]}
    # lfsig used 0.52880 / 0.57300 (predict_cvmap rounding); reproduce those literals so the seed-1 cells are
    # the SAME cells (P4 bit-identity needs the same phi string in make_conf)
    assert abs(phi[0.2] - 0.52880) < 6e-6 and abs(phi[0.5] - 0.57300) < 6e-6, phi
    phi[0.2], phi[0.5] = 0.52880, 0.57300
    pm35 = phim_interp(0.35)
    lin35 = PHIM[0.2] + (PHIM[0.5] - PHIM[0.2]) * (0.35 - 0.2) / 0.3
    phi[0.35] = round(XM * pm35, 5)
    print(f"phi_m(0.35): quadratic {pm35:.5f}, linear {lin35:.5f}, difference {pm35 - lin35:+.5f} "
          f"(= {XM * (pm35 - lin35):+.5f} in phi, {(pm35 - lin35) / PHIM[0.2] * XM:+.4f} in x_m)")
    lines = []
    for seed in SEEDS:
        for sig in ORDER:
            for cv in (0.2, 0.35, 0.5):
                lab = f"c{int(round(cv * 100)):02d}_x90_s{sig:g}_r{seed}"
                cmd = (f"bash ~/lfdem/run_cell2.sh {lab} {cv} {phi[cv]:.5f} {sig:g} {MU} {seed} {STRAIN} "
                       f"> ~/lfdem/runs/lfsig2_{lab}.out 2>&1")
                lines.append(cmd)
    out = HERE / "queue_lfsig2.txt"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {out} ({len(lines)} cells); phi: " + ", ".join(f"CV {c}: {p:.5f}" for c, p in sorted(phi.items())))
    print("first 6:"); print("\n".join(lines[:6]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
