"""Figure 3 of the strain-clock note: constant-stress start-up on silica (a second protocol).

Richards et al. (2019) deposit, 4 um silica in 87 wt% glycerol-water, phi = 0.574.
(a) relative viscosity against strain after each stress step, for a subset of stresses
    below (dashed) and above (solid) the thickening onset of the steady curve.
(b) the fit-free jam strain gamma_j (steepest fall of the rate, filled) and the recoil R
    during the preceding rest (open), against stress; the grey band is the half-recovery
    strain gamma_50 after REVERSAL on Lin's 2 um silica (same material class, different
    sample), for scale.

Numbers come from scripts/probe_startup_clock.py and probe_startup_p4.py (imported).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
for p in (HERE, HERE.parent / "scripts", HERE.parent / "src"):
    sys.path.insert(0, str(p))

import _common as C  # noqa: E402
import fig_clock as FC  # noqa: E402
import probe_startup_clock as S  # noqa: E402
import probe_startup_p4 as S4  # noqa: E402

SLUG = "startup"
SHOW = (1.0, 3.0, 4.5, 6.0, 8.0, 11.0)


def main() -> None:
    _, _, _, eta_s_rows = S.steady()
    eta_s = float(np.median(eta_s_rows))
    rows, onset = S4.jam_table()
    rev = [x for x in FC.clocks()["Fig_2b"] if x["included"]]
    g50 = np.array([x["g50"] for x in rev])

    fig, ax = plt.subplots(1, 2, figsize=(5.8, 2.4), constrained_layout=True)
    a = ax[0]
    out = []
    j = 0
    for k, sig, p in S.step_files():
        if sig not in SHOW:
            continue
        g, gd, eta, _ = S.load(p)
        m = g > 0
        a.plot(g[m], eta[m] / eta_s, color=C.CB[j % len(C.CB)], lw=0.8,
               ls="-" if sig >= onset else "--", label=f"{sig:g}")
        out += [("curve", sig, float(x), float(y)) for x, y in zip(g[m][::20], eta[m][::20] / eta_s)]
        j += 1
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlim(1e-2, 3)
    a.set_xlabel(r"strain $\gamma$ after the step")
    a.set_ylabel(r"$\eta_r$")
    a.legend(title=r"$\sigma$ (Pa)", fontsize=6, title_fontsize=6, loc="upper left", ncol=2)
    a.text(0.95, 0.06, "a", transform=a.transAxes, ha="right", fontweight="bold")

    b = ax[1]
    sig = np.array([r["sigma"] for r in rows])
    gj = np.array([r["gamma_j"] for r in rows])
    R = np.array([r["R"] for r in rows])
    above = np.array([r["above"] for r in rows])
    b.axhspan(g50.min(), g50.max(), color="0.85", lw=0)
    b.plot(sig[above], gj[above], "o", color=C.CB[0], ms=3.5, label=r"$\gamma_j$")
    b.plot(sig[~above], gj[~above], "o", color=C.CB[0], ms=3.5, mfc="none")
    b.plot(sig, R, "s", color=C.CB[1], ms=3, mfc="none", label=r"recoil $R$")
    b.axvline(onset, color="k", ls=":", lw=0.7)
    b.set_yscale("log")
    b.set_ylim(0.03, 3)
    b.set_xlim(0, None)
    b.set_xlabel(r"$\sigma$ (Pa)")
    b.set_ylabel("strain")
    b.legend(loc="center right")
    b.text(0.95, 0.06, "b", transform=b.transAxes, ha="right", fontweight="bold")
    out += [("jam", float(s_), float(g_), float(r_)) for s_, g_, r_ in zip(sig, gj, R)]

    prov = C.provenance(Path(__file__), "Richards2019 DataShare 10283/3357 Figure3/Data/silica87glyc; Lin2015 Fig2b",
                        f"eta_s={eta_s:.2f}mPas onset={onset:g}Pa")
    C.DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(C.DATA_DIR / f"fig_{SLUG}.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# {prov}\n# curve,sigma,strain,eta_r (every 20th) | jam,sigma,gamma_j,R\n")
        csv.writer(fh).writerows(out)
    print(*C.save(fig, SLUG, prov))


if __name__ == "__main__":
    main()
