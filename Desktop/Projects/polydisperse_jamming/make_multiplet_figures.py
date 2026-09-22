"""make_multiplet_figures.py — the aggregation figures.

House style is `make_figures.py`'s and is imported from it rather than restated: minimal
text on the axes (the explanation lives in `derivations/multiplet_psd_2026-09-22.md`),
linear axes start at zero, log only above a decade, and every figure saves .pdf + .png +
the PLOTTED DATA as .csv.

    python make_multiplet_figures.py            # build
    python make_multiplet_figures.py --check    # build twice, compare md5

Figures
-------
A1  THE AGGREGATION PLANE. `Lambda/G` over `(d_f, <k>)`, with the `Lambda = G` crossing
    drawn as a contour and the region where the model is INADMISSIBLE (`phi_eff > 1`)
    hatched out. One figure carries the whole result: the sign of the aggregation effect
    is set by the aggregate's OPENNESS, not by how much it has aggregated, and the model
    destroys itself before the open-aggregate corner is reached.
A2  THE TWO LEVERS SEPARATELY, against `d_f` at three aggregation numbers, with the
    MEASURED ceiling on `G` drawn as a horizontal line. This is the closure-free form:
    `G` never leaves a narrow band that measurement already bounds, while `Lambda` runs
    away, and where the two curves cross is `d_f*`.
"""

from __future__ import annotations

import argparse
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

from make_figures import COL, FIG, md5, save          # noqa: E402
from multiplet_levers import BETA, critical_df, levers  # noqa: E402

# The solid loading the admissibility bound is drawn at. It is the Wyart-Cates solver's
# own default PHI, so the hatched region is the one a DNS run would actually hit.
PHI_SOLID = 0.52

# The MEASURED ceiling on G: the widest pack this project has ever made, lognormal
# CV = 0.5 at this rung, divided by the measured monodisperse beta. Source:
# out/mixture_closure.txt. Not a fit.
G_MAX_MEASURED = 0.6489 / BETA


def fig_A1() -> None:
    """The aggregation plane, with the crossing and the inadmissible corner.

    The colour is `log2(Lambda/G)` so that the crossing sits at 0 and the two signs are
    symmetric about it; plotting the raw ratio would put the interesting boundary at 1
    on a scale running to 25 and hide it.
    """
    dfs = np.linspace(1.2, 3.0, 73)
    kbs = np.linspace(1.05, 5.0, 61)
    ratio = np.zeros((kbs.size, dfs.size))
    lam = np.zeros_like(ratio)
    for i, kb in enumerate(kbs):
        for j, df in enumerate(dfs):
            lv = levers(float(kb), float(df))
            ratio[i, j] = lv["ratio"]
            lam[i, j] = lv["lam"]

    fig, ax = plt.subplots(figsize=(3.6, 2.9), constrained_layout=True)
    # A diverging map centred on the crossing is the right form, but the data is very
    # asymmetric: log2(Lambda/G) runs about -0.05 to +4.7, so a symmetric scale at the
    # full range leaves the ENTIRE suppressing side in one shade and hides the structure
    # around the boundary the figure exists to show. Clipped to +-CLIP with the colorbar
    # marked `extend`, so the clipping is declared rather than silent; the .csv carries
    # the unclipped values.
    CLIP = 2.0
    lo, hi = float(np.log2(ratio).min()), float(np.log2(ratio).max())
    pc = ax.pcolormesh(dfs, kbs, np.log2(ratio), cmap="RdBu_r",
                       vmin=-CLIP, vmax=CLIP, shading="gouraud", rasterized=True)
    cb = fig.colorbar(pc, ax=ax, pad=0.02,
                      extend="max" if hi > CLIP else "neither")
    cb.set_label(r"$\log_2(\Lambda/G)$", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    print(f"      A1: log2(Lambda/G) spans {lo:+.3f} .. {hi:+.3f}, clipped at "
          f"+-{CLIP} for display")

    # the crossing itself
    cs = ax.contour(dfs, kbs, ratio, levels=[1.0], colors=COL[0], linewidths=1.4)
    ax.clabel(cs, fmt={1.0: r"$\Lambda=G$"}, fontsize=7, inline=True)

    # where the effective-sphere model destroys itself: phi_eff = Lambda*phi_solid > 1
    ax.contourf(dfs, kbs, lam * PHI_SOLID, levels=[1.0, 1e9],
                colors="none", hatches=["///"], alpha=0.0)
    ax.contour(dfs, kbs, lam * PHI_SOLID, levels=[1.0], colors=COL[2],
               linewidths=1.0, linestyles="--")

    ax.set_xlabel(r"$d_f$")
    ax.set_ylabel(r"$\langle k\rangle$")
    ax.set_xlim(dfs[0], dfs[-1])
    ax.set_ylim(kbs[0], kbs[-1])
    ax.tick_params(labelsize=8)
    save(fig, "A1_aggregation_plane",
         [dict(d_f=float(dfs[j]), k_bar=float(kbs[i]),
               lambda_over_G=float(ratio[i, j]), Lambda=float(lam[i, j]),
               phi_eff=float(lam[i, j] * PHI_SOLID))
          for i in range(0, kbs.size, 4) for j in range(0, dfs.size, 4)],
         f"Lambda/G over the aggregation plane. Solid line: the Lambda = G crossing. "
         f"Dashed line: phi_eff = Lambda*phi_solid = 1 at phi_solid = {PHI_SOLID}, "
         f"beyond which the no-interpenetration model is INADMISSIBLE. "
         f"Closure: mixture_model refit, beta = {BETA}.")


def fig_A2() -> None:
    """The two levers separately, with the MEASURED ceiling on G.

    This is the closure-free form of the argument in one picture: G is penned into a
    narrow band that measurement already bounds, Lambda is not, and d_f* is where they
    meet. The y axis is log because Lambda spans more than a decade -- the house rule is
    linear-from-zero unless the data spans a decade, and it does.
    """
    dfs = np.linspace(1.2, 3.0, 145)
    fig, ax = plt.subplots(figsize=(3.6, 2.9), constrained_layout=True)
    rows = []
    for c, kb in enumerate((1.2, 2.0, 5.0)):
        lam = np.array([levers(float(kb), float(d))["lam"] for d in dfs])
        g = np.array([levers(float(kb), float(d))["G"] for d in dfs])
        ax.plot(dfs, lam, color=COL[c + 1], lw=1.4,
                label=rf"$\langle k\rangle={kb:g}$")
        ax.plot(dfs, g, color=COL[c + 1], lw=1.0, ls=":")
        d = critical_df(float(kb))
        if np.isfinite(d):
            ax.plot([d], [levers(float(kb), d)["lam"]], "o", color=COL[c + 1], ms=4)
        rows += [dict(k_bar=kb, d_f=float(x), Lambda=float(a), G=float(b))
                 for x, a, b in zip(dfs[::4], lam[::4], g[::4])]

    ax.axhline(G_MAX_MEASURED, color=COL[0], lw=1.0, ls="--")
    ax.set_yscale("log")
    ax.set_xlabel(r"$d_f$")
    ax.set_ylabel(r"$\Lambda$  (solid),   $G$  (dotted)")
    ax.set_xlim(dfs[0], dfs[-1])
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    ax.tick_params(labelsize=8)
    save(fig, "A2_two_levers", rows,
         f"Lambda (solid) and G (dotted) against d_f at three aggregation numbers. "
         f"Markers: the Lambda = G crossing d_f*. Dashed horizontal line: the MEASURED "
         f"ceiling on G, {G_MAX_MEASURED:.4f}, from the widest pack this project has "
         f"made (lognormal CV = 0.5). Log y: Lambda spans more than a decade.")


def build() -> None:
    for f in (fig_A1, fig_A2):
        try:
            f()
            print(f"  {f.__name__} ok")
        except Exception as exc:      # noqa: BLE001
            print(f"  {f.__name__} FAILED: {exc!r}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="build twice and compare md5 -- the cheapest refactor test")
    args = ap.parse_args()
    build()
    if not args.check:
        return 0
    names = ["A1_aggregation_plane", "A2_two_levers"]
    before = {n: md5(os.path.join(FIG, f"{n}.pdf")) for n in names}
    build()
    bad = [n for n in names if md5(os.path.join(FIG, f"{n}.pdf")) != before[n]]
    for n in names:
        print(f"   {n}: {'IDENTICAL' if n not in bad else '🔴 CHANGED'}")
    if bad:
        print(f"🔴 {len(bad)} figure(s) are not byte-reproducible: {bad}")
        return 1
    print("all figures byte-identical on rebuild")
    return 0


if __name__ == "__main__":
    sys.exit(main())
