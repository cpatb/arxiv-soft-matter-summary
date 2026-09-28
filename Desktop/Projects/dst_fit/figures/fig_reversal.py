"""Figure 1 of the strain-clock note: reversal transients and which clock governs.

(a,b) Lin et al. (2015) post-reversal relative viscosity against accumulated strain,
      2 um silica (Fig_2b) and 1.4 um PMMA (Fig_2d), one line per shear rate. Dashed:
      Wyart-Cates's own frictionless floor (1 - phi/phi_0)^-n from the steady fit the
      gate uses (G10/G11). Dotted: the lowest viscosity on the same sample's STEADY
      rate sweep -- a model-free reference.
(c)   the fit-free half-recovery landmark read as a strain (gamma_50, filled) and as a
      time (t_50 = gamma_50 / gamma_dot, open), each normalised by its value at the
      lowest admitted rate. A pure strain clock is flat; a pure time clock falls as
      1/gamma_dot (grey line, slope -1).

All numbers come from scripts/gate_g11.py's own functions (imported, not copied).
Linear y axes start at zero (house rule).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "scripts"))
sys.path.insert(0, str(HERE.parent / "src"))

import _common as C  # noqa: E402
import gate_g11 as G  # noqa: E402
from dst_fit.io_literature import read_lin_transients  # noqa: E402
from dst_fit.transient import frictionless_relative_viscosity  # noqa: E402

SLUG = "reversal"
LABELS = {"Fig_2b": "silica", "Fig_2d": "PMMA"}


def compute():
    out = {}
    for panel, steady_panel in G.PAIRS:
        split, phi, phi_m, sigma_star = G.steady_fit(steady_panel, G.PHI0_DEFAULT)
        trans = read_lin_transients(C.LIN_ROOT, panel)
        floor = frictionless_relative_viscosity(phi, G.PHI0_DEFAULT)
        marks = G.landmarks(trans, phi, phi_m, G.PHI0_DEFAULT)
        records = [(rate, np.asarray(g), np.asarray(e)) for rate, g, e in trans]
        out[panel] = dict(phi=phi, phi_m=phi_m, sigma_star=sigma_star, floor=floor,
                          sweep_min=float(np.min(split.sweep_eta_r)),
                          records=records, marks=marks)
        print(f"{panel}: phi={phi:.4f} phi_m={phi_m:.4f} floor={floor:.3f} "
              f"sweep min eta_r={out[panel]['sweep_min']:.3f} "
              f"N records={len(records)} N landmarks={len(marks)}")
        for rate, g, e in records:
            print(f"   rate {rate:g}: N={g.size} strain {g.min():.4g}-{g.max():.4g} "
                  f"eta_r {e.min():.3f}-{e.max():.3f} min/floor={e.min() / floor:.3f}")
    return out


def main() -> None:
    d = compute()
    fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.3), constrained_layout=True)
    rows = []
    for k, panel in enumerate(("Fig_2b", "Fig_2d")):
        a = ax[k]
        dd = d[panel]
        for j, (rate, g, e) in enumerate(dd["records"]):
            a.plot(g, e, color=C.CB[j], lw=0.8, label=f"{rate:g}")
            step = max(1, g.size // 400)
            rows += [(panel, rate, float(x), float(y)) for x, y in zip(g[::step], e[::step])]
        a.axhline(dd["floor"], color="k", ls="--", lw=0.8)
        a.axhline(dd["sweep_min"], color="k", ls=":", lw=0.8)
        a.set_xscale("log")
        a.set_xlim(1e-2, 10)
        a.set_ylim(0, None)
        a.set_xlabel(r"strain $\gamma$")
        a.set_ylabel(r"$\eta_r$")
        a.legend(title=r"$\dot\gamma$ (1/s)", loc="upper left", title_fontsize=6.5)
        a.text(0.97, 0.05, "abc"[k], transform=a.transAxes, ha="right", fontweight="bold")

    a = ax[2]
    for k, panel in enumerate(("Fig_2b", "Fig_2d")):
        m = d[panel]["marks"]
        r = np.array(sorted(m))
        g50 = np.array([m[x][0] for x in r])
        t50 = g50 / r
        a.plot(r, g50 / g50[0], "o-", color=C.CB[k], ms=3.5, label=f"{LABELS[panel]} $\\gamma_{{50}}$")
        a.plot(r, t50 / t50[0], "s--", color=C.CB[k], ms=3.5, mfc="none",
               label=f"{LABELS[panel]} $t_{{50}}$")
        for x, y1, y2 in zip(r, g50, t50):
            rows.append((panel + "_landmark", float(x), float(y1), float(y2)))
    rr = np.array([0.03, 3.0])
    a.plot(rr, (rr / 0.05) ** -1.0 * 1.0, color="0.6", lw=0.7)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel(r"$\dot\gamma$ (1/s)")
    a.set_ylabel("normalised landmark")
    a.legend(loc="lower left")
    a.text(0.97, 0.9, "c", transform=a.transAxes, ha="right", fontweight="bold")

    prov = C.provenance(Path(__file__), "Lin2015 DataShare 10283/886 Fig2a-d",
                        f"phi_0={G.PHI0_DEFAULT} n=WC_EXPONENT")
    C.DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(C.DATA_DIR / f"fig_{SLUG}.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# {prov}\n# rows: panel,rate,strain,eta_r (decimated) | "
                 f"panel_landmark,rate,gamma_50,t_50\n")
        csv.writer(fh).writerows(rows)
    for panel in ("Fig_2b", "Fig_2d"):
        print(f"{panel} floor={d[panel]['floor']:.4f} sweep_min={d[panel]['sweep_min']:.4f}")
    pdf, png = C.save(fig, SLUG, prov)
    print(pdf, png)


if __name__ == "__main__":
    main()
