"""make_facet_figures.py -- figure FC1, convex facets and the window (PREREG_FACET_2026-09-27, rotations live).

FC1 (a) phi at p_hat 1e-3 against the superball exponent p = 1/eps, N 100, seeds 1-3, one point per seed:
        mu_p 0.3 (orange) and mu_p 0.1 (blue); filled = met the stop rule, open = unconverged, read from the progress
        log's last quarter (Am. 6, descriptive). Grey diamonds: frictionless seed-1 trajectories at their last log
        (never converged, Am. 4). Grey crosses: the WITHDRAWN frozen-rotation pilot (frictionless), for scale.
    (b) W_0.1 = 1 - phi(mu 0.3)/phi(mu 0.1) per seed where both edges exist; open if either edge is unconverged.
    (c) face-face contact enrichment E_ff = f_ff / f_ff0 (area null, facet_census.py), mu_p 0.3, pooled over seeds.
Explanations in FIGURES_FACET.md. Health and primary values come from grade_facet.

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


def edge(recs, e, mu, s):
    """(phi, converged) at N 100 for one cell: stop-rule value, else the Am. 6 log readout, else None."""
    r = recs.get((100, e, mu, s))
    if r is not None and gf.healthy(r):
        return gf.phi_primary(r), True
    if mu == 0.1 and e != 1.0:
        lp = gf.log_phi(os.path.join(gf.LSDIR, f"facet2mu01_e{e}_s{s}.log"), e)
        if lp:
            return lp[0], False
    return None


def pilot_phi0():
    out = {}
    for f in PILOT_FILES:
        p = os.path.join(gf.LSDIR, f)
        if not os.path.exists(p):
            continue
        for line in open(p, encoding="utf-8"):
            r = json.loads(line)
            if r["mu_p"] == 0.0 and r["n"] == 100 and r["nodes"] == 1000 and r["converged"] and r["eps"] < 1.0:
                out.setdefault(r["eps"], []).append(r["phi"])
    return {e: float(np.mean(v)) for e, v in out.items()}


def build() -> None:
    recs = gf.load()
    rows = []
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
    for mu, c in ((0.3, COL[2]), (0.1, COL[1])):
        for e in gf.EPS_LADDER:
            for s in (1, 2, 3):
                v = edge(recs, e, mu, s)
                if v is None:
                    continue
                ax[0].plot(1 / e, v[0], "o", color=c, mfc=c if v[1] else "white", ms=5)
                rows.append(dict(panel="a", mu_p=mu, p=1 / e, seed=s, value=v[0], converged=int(v[1])))
    for e in (0.9, 0.8, 0.7, 0.5):
        lp = gf.log_phi(os.path.join(gf.LSDIR, f"facet2vis_n100_e{e}_m0.0_s1.log"), e)
        if lp:
            ax[0].plot(1 / e, lp[0], "D", color="0.35", mfc="white", ms=5)
            rows.append(dict(panel="a_frictionless_traj", mu_p=0.0, p=1 / e, seed=1, value=lp[0], converged=0))
    s0 = [recs[k] for k in recs if k[:3] == (100, 1.0, 0.0) and gf.healthy(recs[k])]
    if s0:
        v = float(np.mean([gf.phi_primary(r) for r in s0]))
        ax[0].plot(1.0, v, "D", color="0.35", ms=5)
        rows.append(dict(panel="a_frictionless_sphere", mu_p=0.0, p=1.0, seed=0, value=v, converged=1))
    pil = pilot_phi0()
    if pil:
        ax[0].plot([1 / e for e in sorted(pil)], [pil[e] for e in sorted(pil)], "x", color="0.6", ms=6)
        for e, v in pil.items():
            rows.append(dict(panel="a_pilot_withdrawn", mu_p=0.0, p=1 / e, seed=0, value=v, converged=1))
    for e in gf.EPS_LADDER:
        for s in (1, 2, 3):
            a, b = edge(recs, e, 0.1, s), edge(recs, e, 0.3, s)
            if a and b:
                w = 1 - b[0] / a[0]
                ok = a[1] and b[1]
                ax[1].plot(1 / e, w, "o", color=COL[0], mfc=COL[0] if ok else "white", ms=5)
                rows.append(dict(panel="b", mu_p=-1, p=1 / e, seed=s, value=w, converged=int(ok)))
    nulls = {e: facet_census.area_null(e)["f_ff0"] for e in gf.EPS_LADDER}
    es = []
    for e in gf.EPS_LADDER:
        keys = [k for k in recs if k[1] == e and k[2] == 0.3 and gf.healthy(recs[k]) and recs[k].get("census")]
        if not keys:
            continue
        ff = []
        for k in keys:
            z = np.load(os.path.join(gf.LSDIR, recs[k]["census"]))
            ff.append((z["m1"] > gf.C_FACE) & (z["m2"] > gf.C_FACE))
        es.append((1 / e, float(np.mean(np.concatenate(ff))) / nulls[e]))
        rows.append(dict(panel="c", mu_p=0.3, p=1 / e, seed=0, value=es[-1][1], converged=1))
    if es:
        ax[2].plot(*zip(*es), "o-", color=COL[2], ms=5)
    ax[2].axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax[0].set_xlabel("$p = 1/\\varepsilon$"); ax[0].set_ylabel(r"$\phi$ at $\hat p = 10^{-3}$"); ax[0].set_ylim(0, 0.8)
    ax[1].set_xlabel("$p$"); ax[1].set_ylabel(r"$W_{0.1} = 1-\phi_{0.3}/\phi_{0.1}$"); ax[1].set_ylim(0, None)
    ax[2].set_xlabel("$p$"); ax[2].set_ylabel(r"$E_{ff}$"); ax[2].set_ylim(0, None)
    for i, t in enumerate("abc"):
        ax[i].text(0.03, 0.95, f"({t})", transform=ax[i].transAxes, va="top")
    fig.tight_layout()
    save(fig, "FC1_facet_window", rows, "PREREG_FACET_2026-09-27; see FIGURES_FACET.md")


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
