"""Write lfdem/queue_lfsig2x.txt: the x-matched CV 0.5 arm of PREREG_LFSIG2 Amendment 2.

The registered CV 0.2 and CV 0.5 arms share the Yade x_m = 0.90 but sit at different LF-DEM-relative distances
from their own jamming points (lfsig P2 inversion at x_m 0.90: phi_m,LF = 0.6258 for CV 0.2, 0.6585 for CV 0.5;
out/grade_lfsig.txt). This arm re-runs CV 0.5 at the phi that gives it the CV 0.2 arm's LF-relative x:
    x_LF(0.2) = 0.52880 / 0.6258;   phi_x = x_LF(0.2) * 0.6585
on the five rungs {3, 4.5, 7, 10, 100} (the crossing region plus the normalising plateau), seeds 1-3: 15 cells.
Every number is printed. Labels c50x_x85_s<sigma>_r<seed>; run set lfsig2x (separate directory, separate loader key).
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PHI_02 = 0.52880
PHIM_LF = {0.2: 0.6258, 0.5: 0.6585}      # lfsig P2 inversion, x_m 0.90 arms (out/grade_lfsig.txt)
SIGMAS = (4.5, 7, 3, 10, 100)             # decisive first
SEEDS = (1, 2, 3)
STRAIN, MU = 4, 0.3


def main() -> int:
    x_lf = PHI_02 / PHIM_LF[0.2]
    phi_x = round(x_lf * PHIM_LF[0.5], 5)
    print(f"x_LF(CV 0.2 arm) = {x_lf:.5f}; x-matched CV 0.5 phi = {phi_x:.5f} (registered arm: 0.57300; "
          f"difference {phi_x - 0.573:+.5f}, i.e. Yade x_m {phi_x / 0.63667:.4f})")
    lines = []
    for seed in SEEDS:
        for sig in SIGMAS:
            lab = f"c50x_x85_s{sig:g}_r{seed}"
            lines.append(f"bash ~/lfdem/run_cell2x.sh {lab} 0.5 {phi_x:.5f} {sig:g} {MU} {seed} {STRAIN} "
                         f"> ~/lfdem/runs/lfsig2x_{lab}.out 2>&1")
    out = HERE / "queue_lfsig2x.txt"
    out.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {out} ({len(lines)} cells)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
