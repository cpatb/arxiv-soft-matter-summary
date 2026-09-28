"""Figure 2 of the strain-clock note: the measured clock, and the definition control.

(a) Three estimators of the strain clock per shear rate, both materials (Lin 2015):
    the exponential constant fitted in f-space over the whole record (filled, with its
    profile interval), the same fit truncated at gamma = 3 (open), and the fit-free
    half-recovery strain gamma_50 (crosses). Excluded, drawn grey: PMMA at 0.05 1/s,
    whose inferred plateau in f is outside [0, 1] (gate G11 rule). Horizontal line:
    gamma_0 = 0.2, the value the drop-impact solver ran. NOT "Richards' corn starch":
    Richards et al. 2019 (PRL 123 038004) used gamma_0 = O(0.1) (~0.1 by back-solving
    their t_c and epsilon); 0.197 is Han et al. 2018 (PRF 3 073301)'s START-UP scale in
    g = 1 - exp(-gamma/gamma*); simulations use 0.023 (Mari et al. 2015). See
    PRIOR_ART_STRAIN_CLOCK_2026-09-27.md.
(b) Lin's SIMULATION record (Fig_3b, gamma_dot = 0.74 gamma_dot_0) read two ways, each
    normalised to its own fitted start and plateau: the directly evaluated contact
    stress fraction s = sigma_c / sigma, and the Wyart-Cates inversion of the total
    viscosity. Normalisation makes the inverted curve gauge-free (f is affine in u), so
    the curve does not depend on the packing triple. Lines: the fitted exponentials.

Every number is produced by the package's own fit_recovery / half_recovery_strain.
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
import probe_g11_direct_control as P  # noqa: E402
from dst_fit.transient import fit_recovery  # noqa: E402
from dst_fit.transient import frictional_fraction_from_viscosity as inv  # noqa: E402

SLUG = "clock"
LABELS = {"Fig_2b": "silica", "Fig_2d": "PMMA"}


def clocks():
    out = {}
    for panel, steady_panel in G.PAIRS:
        _, trans, phi, phi_m, _, full, _ = G.measure(panel, steady_panel, G.PHI0_DEFAULT)
        _, _, _, _, _, trunc, _ = G.measure(panel, steady_panel, G.PHI0_DEFAULT,
                                            gamma_fit_max=G.TRUNCATE_AT)
        marks = G.landmarks(trans, phi, phi_m, G.PHI0_DEFAULT)
        tr = {r.rate: r for r in trunc}
        rows = []
        for r in full:
            ok = 0.0 <= r.f_plateau <= 1.0
            rows.append(dict(rate=r.rate, full=r.gamma_eff, lo=r.gamma_eff_lo,
                             hi=r.gamma_eff_hi,
                             trunc=tr[r.rate].gamma_eff if r.rate in tr else np.nan,
                             g50=marks.get(r.rate, (np.nan,))[0], included=ok))
        out[panel] = rows
        for x in rows:
            print(f"{panel} {x['rate']:g}: full {x['full']:.4f} [{x['lo']:.4f},{x['hi']:.4f}] "
                  f"trunc {x['trunc']:.4f} g50 {x['g50']:.4f} included={x['included']}")
    inc = [v for rows in out.values() for x in rows if x["included"]
           for v in (x["full"], x["trunc"], x["g50"]) if np.isfinite(v)]
    print(f"included estimates N={len(inc)} min={min(inc):.4f} max={max(inc):.4f} "
          f"min/0.2={min(inc) / G.SOLVER_GAMMA_0:.2f}")
    return out


def control():
    contact = P.load("Fig_3b_contact.csv")
    total = P.load("Fig_3b_total.csv")
    g = contact[:, 0]
    s = contact[:, 1] / total[:, 1]
    eta = total[:, 1]
    phi, phi_m, phi_0 = P.TRIPLES[0]
    r_inv = fit_recovery(g, eta, P.RATE, phi, phi_m, phi_0)
    r_dir = fit_recovery(g, P.wc_forward(s, phi, phi_m, phi_0), P.RATE, phi, phi_m, phi_0)
    f = inv(eta, phi, phi_m, phi_0)
    curves = {}
    for name, y, r in (("direct", s, r_dir), ("inverted", f, r_inv)):
        norm = (y - r.f_start) / (r.f_plateau - r.f_start)
        gg = np.geomspace(r.gamma_fit_min, g[-1], 200)
        model = 1.0 - np.exp(-(gg - r.gamma_fit_min) / r.gamma_eff)
        curves[name] = (norm, gg, model, r.gamma_eff)
        print(f"control {name}: gamma_eff={r.gamma_eff:.4f} start at gamma={r.gamma_fit_min:.4f}")
    return g, curves


def main() -> None:
    cl = clocks()
    g, curves = control()
    fig, ax = plt.subplots(1, 2, figsize=(5.2, 2.3), constrained_layout=True)
    a = ax[0]
    rows_out = []
    for k, panel in enumerate(("Fig_2b", "Fig_2d")):
        for x in cl[panel]:
            col = C.CB[k] if x["included"] else "0.6"
            a.errorbar(x["rate"], x["full"], yerr=[[x["full"] - x["lo"]], [x["hi"] - x["full"]]],
                       fmt="o", color=col, ms=3.5, capsize=1.5, lw=0.7)
            a.plot(x["rate"] * 1.06, x["trunc"], "o", mfc="none", color=col, ms=3.5)
            a.plot(x["rate"] * 0.94, x["g50"], "x", color=col, ms=3.5)
            rows_out.append((panel, x["rate"], x["full"], x["lo"], x["hi"], x["trunc"],
                             x["g50"], int(x["included"])))
        a.plot([], [], "o", color=C.CB[k], ms=3.5, label=LABELS[panel])
    a.axhline(G.SOLVER_GAMMA_0, color="k", ls="--", lw=0.8)
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_ylim(0.1, 3)
    a.set_xlabel(r"$\dot\gamma$ (1/s)")
    a.set_ylabel(r"strain clock $\gamma_0$")
    a.legend(loc="lower left")
    a.text(0.95, 0.9, "a", transform=a.transAxes, ha="right", fontweight="bold")

    b = ax[1]
    step = max(1, g.size // 600)
    for k, name in enumerate(("direct", "inverted")):
        norm, gg, model, ge = curves[name]
        b.plot(g[::step], norm[::step], ".", ms=1.2, color=C.CB[k + 2], alpha=0.5)
        b.plot(gg, model, color=C.CB[k + 2], lw=1.0,
               label=rf"{name}, $\gamma_0$ = {ge:.2f}")
        rows_out += [("control_" + name, float(x), float(y)) for x, y in zip(g[::step], norm[::step])]
    b.set_xscale("log")
    b.set_xlim(1e-2, 10)
    b.set_ylim(-0.5, 1.3)
    b.axhline(0, color="0.8", lw=0.5)
    b.set_xlabel(r"strain $\gamma$")
    b.set_ylabel("normalised recovery")
    b.legend(loc="upper left")
    b.text(0.95, 0.06, "b", transform=b.transAxes, ha="right", fontweight="bold")

    prov = C.provenance(Path(__file__), "Lin2015 DataShare 10283/886 Fig2a-d, Fig3b",
                        f"phi_0={G.PHI0_DEFAULT} truncate={G.TRUNCATE_AT} control_triple={P.TRIPLES[0]}")
    C.DATA_DIR.mkdir(parents=True, exist_ok=True)
    with open(C.DATA_DIR / f"fig_{SLUG}.csv", "w", newline="", encoding="utf-8") as fh:
        fh.write(f"# {prov}\n# panel,rate,full,lo,hi,trunc,g50,included | "
                 f"control_<channel>,strain,normalised (decimated)\n")
        csv.writer(fh).writerows(rows_out)
    print(*C.save(fig, SLUG, prov))


if __name__ == "__main__":
    main()
