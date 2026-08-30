"""Grid-convergence study for the disc-altitude field: observed order and GCI.

WHY THIS REPLACES THE OLD TEST
------------------------------
`validate_cfd.test_grid_independence` compared the production mesh with one coarser
mesh at **five hand-picked stations** and reported the worst of the five as a fraction
of the LOCAL speed. That has three problems, and together they make its FAIL neither
trustworthy nor interpretable:

1. **Two grids give a difference, not an error.** |f_fine - f_coarse| says the two
   meshes disagree; it says nothing about how far the fine mesh is from the continuum
   solution. Only three geometrically similar grids give an observed order of
   convergence and a Richardson extrapolation.

2. **Normalising by the local speed.** The field is heavily sheltered, with
   disc-altitude speeds of 0.02-1.8 m/s. Dividing by a local speed that approaches zero
   in a recirculation makes the metric diverge exactly where the flow is least
   resolved, so the verdict is set by whichever station happened to land in a wake. The
   old test's 57.1% is one station; normalised by u_ref the same station is 10.3%.

3. **A max over five samples** is a noisy statistic, and it under-sampled: over the
   whole pitch at 2 m the fine-vs-coarse difference has p95 = 14.3% and max = 21.2% of
   u_ref, both worse than the five stations suggested.

This script uses the four-member family that moves ONLY the fieldBox refinement level
(so the meshes stay geometrically similar and the refinement ratio is exactly 2):

    gridL1  level 1   8.0  x 2.00 m
    coarse  level 2   4.0  x 1.00 m
    production level 3   2.0  x 0.50 m      <- the mesh the paper reports
    gridL4  level 4   1.0  x 0.25 m

METHOD (Roache; ASME V&V 20-2009)
---------------------------------
For three solutions f1 (finest), f2, f3 (coarsest) with constant refinement ratio r:

    eps21 = f2 - f1,  eps32 = f3 - f2
    p     = ln|eps32 / eps21| / ln r                    observed order
    f_ext = f1 + (f1 - f2) / (r^p - 1)                  Richardson extrapolation
    e_a21 = |(f1 - f2) / f1|                            approximate relative error
    GCI21 = Fs * e_a21 / (r^p - 1),  Fs = 1.25          error band on the FINE grid

The sign of `eps32/eps21` matters and is checked, not assumed. A negative ratio means
the solutions are oscillating rather than converging monotonically; `p` is then not
defined by the formula above and the GCI must not be quoted as if it were. This is a
separated, limit-cycle-averaged flow, so oscillatory behaviour is a real possibility
and is reported rather than hidden.

With four grids there are two independent triplets. Comparing their `p` values tests
whether the fine end of the family is in the asymptotic range at all -- a check a
three-grid study cannot make.

Usage:
    python grid_convergence.py --heading 90 --out grid_convergence.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "P5_Synthesis"))
from conventions import (disc_altitude_index, horizontal_speed,  # noqa: E402
                         pitch_mask)

# Family members, coarsest first. `h` is the fieldBox horizontal cell size in metres;
# the vertical size scales with it by the same factor, so the cells stay geometrically
# similar and a single h characterises the refinement.
FAMILY = [
    ("gridL1", 8.0, "flow_fields_dir_090_gridL1.npz"),
    ("coarse", 4.0, "flow_fields_coarse090.npz"),
    ("production", 2.0, "flow_fields.npz"),
    ("gridL4", 1.0, "flow_fields_dir_090_gridL4.npz"),
]
FS = 1.25          # Roache's safety factor for a three-grid study
R_REFINE = 2.0     # exact, by construction of the family


# COST 732 / AIJ practice: at least this many cells across a building's narrowest
# horizontal dimension before its wake can be considered resolved. Below it, a mesh is
# not a coarser member of a convergence family -- it is a different, unresolved flow.
COST732_MIN_CELLS_PER_BUILDING = 10


def building_resolution(F, heading, h_list):
    """Cells across the narrowest horizontal side of each large building, per grid.

    This is the context that decides whether a family member belongs in a Richardson
    fit at all. Campus Towers -- the 45 m block nine metres off the eastern sideline
    that the paper's headline claim is about -- is only 24 m across its narrow
    dimension, so at h = 8 m and h = 4 m it spans 3 and 6 cells: below the accepted
    guideline, and therefore not a resolved solution of the same problem.
    """
    from scipy import ndimage
    x, y = F["x"], F["y"]
    sol = F[f"dir_{int(heading):03d}__solid"].astype(bool)
    lab, n = ndimage.label(sol)
    sizes = ndimage.sum(sol, lab, range(1, n + 1))
    out = []
    for k in np.argsort(sizes)[::-1][:6]:
        m = lab == k + 1
        xs, ys = x[m.any(axis=1)], y[m.any(axis=0)]
        ex = float(xs.max() - xs.min() + 1)
        ey = float(ys.max() - ys.min() + 1)
        mn = min(ex, ey)
        out.append(dict(area_cells=int(sizes[k]), extent_x=ex, extent_y=ey,
                        min_side_m=mn, centroid=[float(xs.mean()), float(ys.mean())],
                        cells_across={f"h={h:g}": mn / h for h in h_list}))
    return out


def load_member(fname, heading):
    p = HERE / fname
    if not p.exists():
        return None
    F = np.load(p, allow_pickle=False)
    key = f"dir_{int(heading):03d}__U"
    if key not in F.files:
        return None
    return F


def quantities(F, heading):
    """The three aggregates the paper's claims actually rest on.

    `U_mean`   pitch-mean horizontal speed at 2 m -- the sheltering headline
               (the UfieldMeanH macros), and what wind_factor normalises by u_ref.
    `k_mean`   pitch-mean turbulent kinetic energy at 2 m.
    `I_median` MEDIAN over pitch cells of sqrt(2k/3)/|U_h| -- the turbulence
               intensity that drives turbulence_factor, and hence the 270-is-worst
               verdict.

    The MEDIAN, not the mean, of the per-cell intensity: the mean-of-ratio has exactly
    the local-normalisation fragility this script exists to criticise. Over this pitch
    it swings 0.33 -> 0.56 across the family while the median moves only 0.23 -> 0.22,
    because a handful of near-stagnant cells dominate the mean of a ratio whose
    denominator approaches zero.
    """
    zk, _, ok = disc_altitude_index(F["z"], 2.0)
    if not ok:
        raise ValueError("field is not sampled at 2 m")
    u = horizontal_speed(F[f"dir_{int(heading):03d}__U"][:, :, zk, :])
    k = F[f"dir_{int(heading):03d}__k"][:, :, zk]
    solid = F[f"dir_{int(heading):03d}__solid"].astype(bool)
    fx, fy = pitch_mask(F["x"], F["y"])
    pm = np.zeros_like(solid)
    pm[np.ix_(fx, fy)] = True
    sel = pm & ~solid
    per_cell_I = np.sqrt(2 * np.maximum(k[sel], 0.0) / 3.0) / np.maximum(u[sel], 1e-3)
    return {
        "U_mean": float(u[sel].mean()),
        "k_mean": float(k[sel].mean()),
        "I_median": float(np.median(per_cell_I)),
    }, u, sel


def scalar_and_field(F, heading):
    """The quantity the paper reports, plus the full plane behind it.

    The reported scalar is the mean horizontal speed at 2 m over the playing surface,
    excluding cells inside buildings -- the same definition `make_results_tex.py` uses
    for \\UfieldMeanH*, so the GCI applies to the published number itself rather than
    to a proxy.
    """
    zk, _, ok = disc_altitude_index(F["z"], 2.0)
    if not ok:
        raise ValueError("field is not sampled at 2 m")
    u = horizontal_speed(F[f"dir_{int(heading):03d}__U"][:, :, zk, :])
    solid = F[f"dir_{int(heading):03d}__solid"].astype(bool)
    fx, fy = pitch_mask(F["x"], F["y"])
    pm = np.zeros_like(solid)
    pm[np.ix_(fx, fy)] = True
    sel = pm & ~solid
    return float(u[sel].mean()), u, sel


def triplet(f1, f2, f3, r=R_REFINE, fs=FS):
    """GCI on the FINEST of three solutions. f1 finest, f3 coarsest.

    Returns a dict; `valid` is False when the solutions are not monotonically
    converging, in which case p, f_ext and gci are not meaningful and are returned as
    NaN rather than as a number a reader might quote.
    """
    eps21 = f2 - f1
    eps32 = f3 - f2
    out = dict(f1=f1, f2=f2, f3=f3, eps21=eps21, eps32=eps32)

    if eps21 == 0.0:
        out.update(valid=False, reason="f1 == f2 exactly; p undefined",
                   p=float("nan"), f_ext=float("nan"), gci=float("nan"),
                   e_a21=0.0, ratio=float("nan"))
        return out

    ratio = eps32 / eps21
    out["ratio"] = ratio
    if ratio <= 0:
        out.update(
            valid=False,
            reason=("oscillatory: eps32/eps21 < 0, the solutions do not converge "
                    "monotonically, so the observed order and GCI are undefined"),
            p=float("nan"), f_ext=float("nan"), gci=float("nan"),
            e_a21=abs(eps21 / f1))
        return out

    p = np.log(abs(ratio)) / np.log(r)
    e_a21 = abs(eps21 / f1)
    denom = r ** p - 1.0
    if denom <= 0:
        out.update(valid=False, reason=f"r^p - 1 <= 0 (p = {p:.3f})",
                   p=p, f_ext=float("nan"), gci=float("nan"), e_a21=e_a21)
        return out

    f_ext = f1 + (f1 - f2) / denom
    gci = fs * e_a21 / denom
    out.update(valid=True, reason="monotone", p=p, f_ext=f_ext, gci=gci,
               e_a21=e_a21)
    return out


def pointwise_gci(u_list, sel, r=R_REFINE, fs=FS):
    """Cellwise GCI over the pitch, using the three finest available solutions.

    Reported as a distribution rather than a single number: the field is not uniformly
    resolved, and a max over a separated wake is not representative of the pitch.
    """
    u1, u2, u3 = u_list[-1], u_list[-2], u_list[-3]
    a, b, c = u1[sel], u2[sel], u3[sel]
    e21, e32 = b - a, c - b
    ok = (np.abs(e21) > 1e-9) & (e32 / np.where(np.abs(e21) < 1e-30, np.nan, e21) > 0)
    p = np.full(a.shape, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        p[ok] = np.log(np.abs(e32[ok] / e21[ok])) / np.log(r)
    # Clamp p into a physically sensible band before forming GCI. Cellwise p in a
    # separated flow is extremely noisy; values outside [0.5, 3] indicate the cell is
    # not in the asymptotic range, and are reported as such rather than propagated.
    good = ok & np.isfinite(p) & (p >= 0.5) & (p <= 3.0)
    gci = np.full(a.shape, np.nan)
    denom = r ** p[good] - 1.0
    gci[good] = fs * np.abs(e21[good] / np.maximum(np.abs(a[good]), 1e-9)) / denom
    return dict(
        n_cells=int(a.size),
        frac_monotone=float(ok.mean()),
        frac_asymptotic=float(good.mean()),
        p_median=float(np.nanmedian(p[good])) if good.any() else float("nan"),
        gci_median=float(np.nanmedian(gci[good])) if good.any() else float("nan"),
        gci_p95=float(np.nanpercentile(gci[good], 95)) if good.any() else float("nan"),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--heading", type=int, default=90)
    ap.add_argument("--out", default="grid_convergence.json")
    args = ap.parse_args()

    print("=" * 76)
    print(f"GRID CONVERGENCE at disc altitude, heading {args.heading} deg")
    print("=" * 76)

    members, scalars, fields, sel_ref = [], [], [], None
    member_fields = []
    for name, h, fname in FAMILY:
        F = load_member(fname, args.heading)
        if F is None:
            print(f"  {name:11s} h = {h:4.1f} m   MISSING ({fname})")
            continue
        s, u, sel = scalar_and_field(F, args.heading)
        if sel_ref is None:
            sel_ref = sel
        else:
            assert sel.shape == sel_ref.shape, "grids differ between members"
        members.append((name, h))
        member_fields.append(F)
        scalars.append(s)
        fields.append(u)
        print(f"  {name:11s} h = {h:4.1f} m   pitch-mean |U_h| at 2 m = {s:.4f} m/s")

    if len(members) < 3:
        print("\n  need at least three family members for a convergence study; "
              f"have {len(members)}. Run the missing cases first.")
        return

    # --- is each member even resolving the buildings? ------------------------------
    # This decides whether a member belongs in a Richardson fit at all, and it is the
    # context the original two-mesh test lacked.
    Fref = load_member(FAMILY[2][2], args.heading)
    br = None
    if Fref is not None:
        h_all = [h for _, h, _ in FAMILY]
        br = building_resolution(Fref, args.heading, h_all)
        print("\n  --- cells across each large building's narrow side ---")
        print("      (COST 732 / AIJ practice: at least "
              f"{COST732_MIN_CELLS_PER_BUILDING} before a wake is resolved)")
        print("      min side   " + "".join(f"h={h:<5g}" for h in h_all))
        for b in br:
            cx, cy = b["centroid"]
            tag = "  <- Campus Towers" if (50 < cx < 90 and -30 < cy < 30) else ""
            cells = "".join(f"{b['min_side_m'] / h:<7.1f}" for h in h_all)
            print(f"      {b['min_side_m']:6.0f} m   {cells}{tag}")
        below = [h for h in h_all
                 if any(b["min_side_m"] / h < COST732_MIN_CELLS_PER_BUILDING
                        for b in br[:5])]
        if below:
            print(f"      -> h = {below} m leave at least one major building below "
                  "the guideline.")
            print("         Those are not coarser solutions of the same problem, they "
                  "are UNRESOLVED")
            print("         ones, and a Richardson fit spanning them crosses a regime "
                  "change.")

    # --- all three claim-bearing quantities, not just the mean speed ----------------
    # The sheltering headline rests on U_mean, but the 270-is-worst verdict rests on
    # the turbulence factor (1 - I), so a study that convergence-tests only the speed
    # says nothing about the claim the paper actually resolves.
    qnames = ["U_mean", "k_mean", "I_median"]
    qvals = {q: [] for q in qnames}
    for (name, _h), F in zip(members, member_fields):
        qq, _, _ = quantities(F, args.heading)
        for q in qnames:
            qvals[q].append(qq[q])

    print("\n  --- the three quantities the claims rest on ---")
    print(f"      {'grid':11s} " + "".join(f"{q:>12s}" for q in qnames))
    for i, (name, _h) in enumerate(members):
        print(f"      {name:11s} "
              + "".join(f"{qvals[q][i]:12.4f}" for q in qnames))
    if len(members) >= 3:
        print(f"\n      {'quantity':10s} {'ratio':>8s} {'p':>7s} {'f_ext':>10s} "
              f"{'err(production)':>16s} {'err(finest)':>12s}")
        quantity_fits = {}
        for q in qnames:
            v = qvals[q]
            t = triplet(v[-1], v[-2], v[-3])
            quantity_fits[q] = t
            if not t["valid"]:
                print(f"      {q:10s} {t['ratio']:8.3f} {'--':>7s} {'--':>10s} "
                      f"{'oscillatory: not converging':>29s}")
                continue
            fe = t["f_ext"]
            prod = v[members.index(("production", 2.0))] \
                if ("production", 2.0) in members else float("nan")
            ep = abs(prod - fe) / fe * 100 if fe else float("nan")
            ef = abs(v[-1] - fe) / fe * 100 if fe else float("nan")
            # An observed order far below the scheme's formal order means the family is
            # NOT in the asymptotic range for this quantity. Richardson then projects a
            # large change from small grid-to-grid differences, and the extrapolated
            # error is an artefact of the tiny exponent rather than a measurement. Flag
            # it instead of printing a number a reader would quote.
            t["asymptotic"] = bool(0.5 <= t["p"] <= 3.0)
            flag = "" if t["asymptotic"] else "  <- p far from formal order:"
            print(f"      {q:10s} {t['ratio']:8.3f} {t['p']:7.2f} {fe:10.4f} "
                  f"{ep:15.1f}% {ef:11.1f}%{flag}")
            if not t["asymptotic"]:
                adj = abs(v[-1] - prod) / abs(prod) * 100
                print(f"      {'':10s} extrapolation NOT trustworthy. What is "
                      f"measured: the production->finest")
                print(f"      {'':10s} change is {adj:.1f}%.")

    print(f"\n  refinement ratio r = {R_REFINE:g} (exact, by construction)")
    print(f"  safety factor  Fs = {FS:g}")

    # Triplets, finest-last ordering in `scalars`.
    results = []
    for i in range(len(scalars) - 2):
        f3, f2, f1 = scalars[i], scalars[i + 1], scalars[i + 2]
        n3, n2, n1 = members[i][0], members[i + 1][0], members[i + 2][0]
        t = triplet(f1, f2, f3)
        t["grids"] = [n1, n2, n3]
        t["h"] = [members[i + 2][1], members[i + 1][1], members[i][1]]
        results.append(t)
        print(f"\n  --- triplet ({n1} / {n2} / {n3}) -- GCI on {n1} ---")
        print(f"      f = {f1:.4f}, {f2:.4f}, {f3:.4f} m/s")
        print(f"      eps21 = {t['eps21']:+.5f}   eps32 = {t['eps32']:+.5f}   "
              f"ratio = {t['ratio']:+.3f}")
        if not t["valid"]:
            print(f"      -> NOT VALID: {t['reason']}")
            print("      -> no observed order, no GCI; the family is not converging "
                  "monotonically here")
            continue
        print(f"      observed order p   = {t['p']:.3f}")
        print(f"      Richardson f_ext   = {t['f_ext']:.4f} m/s")
        print(f"      approx rel error   = {100 * t['e_a21']:.2f}%")
        print(f"      GCI on {n1:11s} = {100 * t['gci']:.2f}%")

    valid = [t for t in results if t["valid"]]
    if len(valid) >= 2:
        ps = [t["p"] for t in valid]
        print(f"\n  asymptotic-range check: observed orders {', '.join(f'{v:.2f}' for v in ps)}")
        spread = max(ps) - min(ps)
        print(f"      spread {spread:.2f} -> "
              f"{'consistent' if spread < 0.5 else 'INCONSISTENT: the family is not '
                 'in the asymptotic range, and the GCI is an estimate at best'}")

    pw = None
    if len(fields) >= 3:
        pw = pointwise_gci(fields, sel_ref)
        print("\n  --- cellwise over the pitch (three finest grids) ---")
        print(f"      cells                         {pw['n_cells']}")
        print(f"      monotonically converging      {100 * pw['frac_monotone']:.1f}%")
        print(f"      in a usable asymptotic range  {100 * pw['frac_asymptotic']:.1f}%")
        print(f"      median observed order p       {pw['p_median']:.2f}")
        print(f"      median GCI                    {100 * pw['gci_median']:.2f}%")
        print(f"      p95    GCI                    {100 * pw['gci_p95']:.2f}%")

    payload = dict(
        heading=args.heading,
        members=[dict(name=n, h=h, pitch_mean_U=s)
                 for (n, h), s in zip(members, scalars)],
        refinement_ratio=R_REFINE, safety_factor=FS,
        triplets=results, pointwise=pw,
        building_resolution=br,
        cost732_min_cells_per_building=COST732_MIN_CELLS_PER_BUILDING,
        quantity="mean horizontal |U| at z = 2 m over the playing surface, "
                 "open cells only (the definition behind \\UfieldMeanH*)",
    )
    Path(args.out).write_text(json.dumps(payload, indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
