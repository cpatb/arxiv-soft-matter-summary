"""probe_backbone_bound.py -- EXPLORATORY. Is a backbone-volume mechanism the right SIZE?

STATUS 2026-09-13: NOT PURSUED, kept as a recorded negative. The candidate mechanism was
that polydispersity moves phi by burying small grains as rattlers, so the LOAD-BEARING
solid fraction phi*(1 - f_rat_vol) is the invariant. This probe brackets f_rat_vol two
ways and finds the bracket far too wide to decide anything -- 0.06 to 0.45 for the CV 0.5
lognormal -- so the question needs the rattler IDENTITIES, which are not recorded, i.e. a
re-run with a new observable. The composition campaign (PREREG_MIXTURE_2026-09-13.md) was
run instead because it settles a sharper question with the measurements already defined.
The one thing this probe does establish: under the SIZE-ORDERED assignment the family
spread shrinks only 0.0447 -> 0.0242 at CV 0.5, so the rattler-volume account cannot
collapse the families at that end of its bracket.


NOT A MEASUREMENT. The rattler IDENTITIES are not on disk -- only the rattler NUMBER
fraction is -- so this script asks what the volume fraction of the rattler set WOULD be
under two bracketing assignments, and whether the resulting backbone solid fraction
phi_bb = phi * (1 - f_rat_vol) is closer to collapsing across families than phi is.

    SIZE-ORDERED   the rattlers are exactly the smallest f_n of the grains by number.
                   This is the LOWER bound on f_rat_vol (any other set of the same
                   number carries more volume) and therefore the UPPER bound on phi_bb.
    UNBIASED       rattler status is independent of size, f_rat_vol = f_n.
                   UPPER bound on f_rat_vol, LOWER bound on phi_bb.

The truth is strictly between them, because small grains are preferentially rattlers but
not exclusively. If the family spread in phi_bb fails to shrink under BOTH bounds, the
mechanism is dead before any DEM runs. If it shrinks under one and not the other, the
measurement decides. Run before pre-registering anything.
"""
from __future__ import annotations

import collections

import numpy as np

import analyse_phase2 as a
from psd import make_psd

TRUNC = (0.2, 5.0)
NDRAW = 400_000
RNG_SEED = 20260913


def sampled_diameters(family: str, cv: float) -> np.ndarray:
    """The truncated sample the packings were drawn from, at large N."""
    rng = np.random.default_rng(RNG_SEED)
    p = make_psd(family, d10=1.0, cv=cv)
    return np.asarray(p.sample(NDRAW, rng, truncate=TRUNC), dtype=float)


def vol_frac_of_smallest(d: np.ndarray, f_n: float) -> float:
    """Volume fraction carried by the smallest `f_n` of the grains BY NUMBER."""
    if f_n <= 0:
        return 0.0
    ds = np.sort(d)
    k = int(round(f_n * ds.size))
    v = ds ** 3
    return float(v[:k].sum() / v.sum())


def main() -> None:
    rows = a.load_all()
    sel = [r for r in rows
           if r.get("n_requested") == 1500
           and abs(r.get("mu_p", 9) - 0.3) < 1e-9
           and abs(r.get("p_hat_target", 0) - 1e-3) < 1e-12
           and r.get("converged") and r.get("eta_roll") is None]
    g: dict[tuple[str, float], list[dict]] = collections.defaultdict(list)
    for r in sel:
        g[(r.get("family") or "lognormal", round(float(r["cv"]), 3))].append(r)

    print("=== probe_backbone_bound.py -- EXPLORATORY, no DEM run ===")
    print(f"   isotropic cells, N=1500, mu_p=0.3, p_hat=1e-3, converged; {len(sel)} cells")
    print("   f_n is MEASURED (iterated k=2 rattlers). f_vol is ASSUMED, two ways.")
    print()
    hdr = (f"{'family':<11}{'CV':>5}{'phi':>8}{'f_n':>8}{'f_vol_ord':>11}"
           f"{'f_vol_unb':>11}{'phi_bb_ord':>12}{'phi_bb_unb':>12}")
    print(hdr)
    out: dict[tuple[str, float], tuple[float, float, float]] = {}
    for k in sorted(g, key=lambda t: (t[0], t[1])):
        fam, cv = k
        rs = g[k]
        phi = float(np.mean([r["phi"] for r in rs]))
        f_n = float(np.mean([r["rattler_frac"] for r in rs]))
        d = sampled_diameters(fam, cv)
        f_ord = vol_frac_of_smallest(d, f_n)
        bb_ord = phi * (1.0 - f_ord)
        bb_unb = phi * (1.0 - f_n)
        out[k] = (phi, bb_ord, bb_unb)
        print(f"{fam:<11}{cv:5.2f}{phi:8.4f}{f_n:8.4f}{f_ord:11.4f}"
              f"{f_n:11.4f}{bb_ord:12.4f}{bb_unb:12.4f}")

    print()
    print("=== does the FAMILY SPREAD shrink? (max-min over the three families) ===")
    print(f"{'CV':>5}{'spread phi':>13}{'spread bb_ord':>16}{'spread bb_unb':>16}")
    for cv in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5):
        fams = [f for f in ("lognormal", "uniform", "bidisperse") if (f, cv) in out]
        if len(fams) < 3:
            continue
        p = [out[(f, cv)][0] for f in fams]
        bo = [out[(f, cv)][1] for f in fams]
        bu = [out[(f, cv)][2] for f in fams]
        print(f"{cv:5.2f}{max(p)-min(p):13.4f}{max(bo)-min(bo):16.4f}"
              f"{max(bu)-min(bu):16.4f}")

    print()
    print("=== the headline contrast, CV 0 -> 0.5, per family ===")
    print(f"{'family':<11}{'dphi':>9}{'d bb_ord':>11}{'d bb_unb':>11}")
    for fam in ("lognormal", "uniform", "bidisperse"):
        if (fam, 0.0) not in out or (fam, 0.5) not in out:
            continue
        lo, hi = out[(fam, 0.0)], out[(fam, 0.5)]
        print(f"{fam:<11}{hi[0]-lo[0]:9.4f}{hi[1]-lo[1]:11.4f}{hi[2]-lo[2]:11.4f}")
    print()
    print("   COLLAPSE would mean the family spread near 0 and the three dphi equal.")
    print("   Read this as a MAGNITUDE check only: the rattler identities are assumed.")


if __name__ == "__main__":
    main()
