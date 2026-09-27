"""make_facet_figures.py -- figure FC1, convex facets and the window (PREREG_FACET_2026-09-27, rotations live).

FC1 (a) phi at p_hat 1e-3 against the superball exponent p = 1/eps: frictionless (phi_0) and mu_p 0.3 (phi_m),
        N 100 seed means (filled, bars = seed sd) and N 200 (open). Grey crosses: the WITHDRAWN frozen-rotation pilot
        (sphericalBodies=True, N 100), drawn to show the size of the artefact.
    (b) window width W = 1 - phi_m/phi_0 per seed-mean, same symbols; linear axis from zero.
    (c) face-face contact enrichment E_ff = f_ff / f_ff0 (area null, facet_census.py), pooled over seeds.
Explanations in FIGURES_FACET.md. Health and primary values come from grade_facet (stop-rule phi).

    python make_facet_figures.py [--check]
"""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

import facet_census               # noqa: E402
import grade_facet as gf          # noqa: E402
from make_figures import COL, FIG, md5, save   # noqa: E402

PILOT_FILES = ("ls_ladder.jsonl", "ls_repl.jsonl", "ls_extra.jsonl")


def seedmean(recs, n, e, m):
    v = [gf.phi_primary(recs[k]) for k in recs if k[:3] == (n, e, m) and gf.healthy(recs[k])]
    return (float(np.mean(v)), float(np.std(v, ddof=1)) if len(v) > 1 else 0.0, len(v)) if v else None


def pilot_phi0():
    out = {}
    for f in PILOT_FILES:
        p = os.path.join(gf.LSDIR, f)
        if not os.path.exists(p):
            continue
        for line in open(p, encoding="utf-8"):
            r = json.loads(line)
            if r["mu_p"] == 0.0 and r["n"] == 100 and r["nodes"] == 1000 and r["converged"]:
                out.setdefault(r["eps"], []).append(r["phi"])
    return {e: float(np.mean(v)) for e, v in out.items()}


def build() -> None:
    recs = gf.load()
    rows = []
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
    for n, filled in ((100, True), (200, False)):
        for m, c, lab in ((0.0, COL[1], r"$\mu_p$ 0"), (0.3, COL[2], r"$\mu_p$ 0.3")):
            pts = [(1 / e, seedmean(recs, n, e, m)) for e in gf.EPS_LADDER]
            pts = [(p, s) for p, s in pts if s]
            if not pts:
                continue
            ax[0].errorbar([p for p, _ in pts], [s[0] for _, s in pts], yerr=[s[1] for _, s in pts], fmt="o-" if filled else "s--",
                           color=c, mfc=c if filled else "white", ms=5, capsize=2,
                           label=f"{lab}, N {n}")
            for p, s in pts:
                rows.append(dict(panel="a", n=n, mu_p=m, p=p, phi=s[0], sd=s[1], n_seeds=s[2]))
        ws = []
        for e in gf.EPS_LADDER:
            s0, sm = seedmean(recs, n, e, 0.0), seedmean(recs, n, e, 0.3)
            if s0 and sm:
                ws.append((1 / e, 1 - sm[0] / s0[0]))
                rows.append(dict(panel="b", n=n, mu_p=-1, p=1 / e, phi=1 - sm[0] / s0[0], sd=0.0, n_seeds=min(s0[2], sm[2])))
        if ws:
            ax[1].plot(*zip(*ws), "o-" if filled else "s--", color=COL[0], mfc=COL[0] if filled else "white", ms=5,
                       label=f"N {n}")
    pil = pilot_phi0()
    if pil:
        ax[0].plot([1 / e for e in sorted(pil)], [pil[e] for e in sorted(pil)], "x:", color="0.6", ms=6,
                   label="frozen rotations (withdrawn)")
        for e, v in pil.items():
            rows.append(dict(panel="a_pilot", n=100, mu_p=0.0, p=1 / e, phi=v, sd=0.0, n_seeds=0))
    nulls = {e: facet_census.area_null(e)["f_ff0"] for e in gf.EPS_LADDER}
    for m, c in ((0.0, COL[1]), (0.3, COL[2])):
        es = []
        for e in gf.EPS_LADDER:
            keys = [k for k in recs if k[1] == e and k[2] == m and gf.healthy(recs[k]) and recs[k].get("census")]
            if not keys:
                continue
            ff = []
            for k in keys:
                z = np.load(os.path.join(gf.LSDIR, recs[k]["census"]))
                ff.append((z["m1"] > gf.C_FACE) & (z["m2"] > gf.C_FACE))
            E = float(np.mean(np.concatenate(ff))) / nulls[e]
            es.append((1 / e, E))
            rows.append(dict(panel="c", n=0, mu_p=m, p=1 / e, phi=E, sd=0.0, n_seeds=len(keys)))
        if es:
            ax[2].plot(*zip(*es), "o-", color=c, ms=5)
    ax[2].axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax[0].set_xlabel("$p = 1/\\varepsilon$"); ax[0].set_ylabel(r"$\phi$ at $\hat p = 10^{-3}$")
    ax[1].set_xlabel("$p$"); ax[1].set_ylabel(r"$W = 1-\phi_m/\phi_0$"); ax[1].set_ylim(0, None)
    ax[2].set_xlabel("$p$"); ax[2].set_ylabel(r"$E_{ff}$"); ax[2].set_ylim(0, None)
    ax[0].legend(fontsize=6, frameon=False)
    for i, t in enumerate("abc"):
        ax[i].text(0.03, 0.95, f"({t})", transform=ax[i].transAxes, va="top")
    fig.tight_layout()
    save(fig, "FC1_facet_window", rows, "PREREG_FACET_2026-09-27; stop-rule phi; see FIGURES_FACET.md")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true", help="regenerate twice and compare PNG md5")
    a = ap.parse_args()
    build()
    if a.check:
        h1 = md5(os.path.join(FIG, "FC1_facet_window.png"))
        build()
        h2 = md5(os.path.join(FIG, "FC1_facet_window.png"))
        print(f"FC1 png md5 {h1} / {h2}: {'IDENTICAL' if h1 == h2 else 'DIFFERS'}")
        return 0 if h1 == h2 else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
