"""Turn the sweep into the theory's answers, and grade the theory against them.

FOUR QUESTIONS, IN ORDER OF HOW MUCH THEY COULD EMBARRASS THE THEORY
--------------------------------------------------------------------
1. Does the gust-driven landing dispersion grow LINEARLY with wind speed, as
   `wind_scaling.py` + `dispersion_law.py` jointly predict, and where does it
   stop?
2. Does the measured dispersion match the ballistic law's MAGNITUDE?
3. Do the two laws' discriminating EXPONENTS -- sigma ~ V_d^-3 vs V_d^-1, and
   sigma ~ L_f^1.7 vs L_f^0.5 -- come out right? These are prefactor-free, so
   they survive any disagreement found in (2) and they are what actually
   separates the two laws. This is the brief's WP2 P2.
4. What fraction of a thrower's total error is the wind at all?

Question 4 is deliberately last and it is the one with the most practical
content: at this site the release scatter alone is ~7 m, so a result about
wind that never compares itself to that floor is describing a correction, not
an effect.

ZERO-BASED AXES. Every linear axis in the figures starts at zero. A wind-speed
response plotted from its own data range makes any curve look dramatic; the
whole point here is which effects are large compared with the floor. Log axes
are used only where the quantity spans more than a decade, and the figure's
sibling .md says so.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import dispersion_law as dl                                # noqa: E402

BG = "#ffffff"


def load(path):
    d = json.loads(Path(path).read_text())
    return d["meta"], d["rows"]


def group(rows, *keys):
    g = defaultdict(list)
    for r in rows:
        # `.get`, not `[]`: the spin axis was added to the sweep after the
        # first speed sweep had already been written, so older result files
        # have no `spin` key and must still analyse.
        g[tuple(r.get(k) for k in keys)].append(r)
    for k in g:
        g[k].sort(key=lambda r: r["u_ref"])
    return g


# ---------------------------------------------------------------------------
# Q1  the wind-speed law
# ---------------------------------------------------------------------------

def wind_speed_law(rows, u_min=1.0):
    """Fit sigma_gust ~ Uref^p and bias ~ Uref^q per (heading, throw, spin).

    The still-air cell is excluded because a power law cannot be fitted through
    Uref = 0; it is kept in the data as the floor and reported separately.
    `u_min` also drops the very lowest speeds where the paired gust
    displacement is comparable to its own standard error.
    """
    out = []
    for (h, th, sp), rs in group(rows, "heading", "throw", "spin").items():
        u = np.array([r["u_ref"] for r in rs])
        m = u >= u_min
        if m.sum() < 3 or len(np.unique(u[m])) < 3:
            continue
        rec = dict(heading=h, throw=th, spin=sp, u_ref=list(u[m]))
        for name, key in (("sigma_gust", "gust__sigma"),
                          ("bias_wind", "wind__bias"),
                          ("sigma_wind", "wind__sigma"),
                          ("sigma_total_extra", "total__sigma")):
            y = np.array([r[key] for r in rs])[m]
            p, lo, hi, r2 = dl.fit_exponent(u[m], y)
            rec[name] = dict(exponent=p, ci=[lo, hi], r2=r2,
                             values=[float(v) for v in y])
            rec[f"{name}_local"] = [float(v) for v in
                                    dl.local_exponent(u[m], np.maximum(y, 1e-9))]
        still = [r for r in rs if r["u_ref"] == 0.0]
        rec["sigma_still"] = float(still[0]["still__sigma"]) if still else float("nan")
        out.append(rec)
    return out


# ---------------------------------------------------------------------------
# Q2  magnitude
# ---------------------------------------------------------------------------

def law_comparison(rows, kappa=dl.KAPPA_ANISOTROPIC):
    """Measured gust dispersion against both laws, per cell.

    Uses the PATH-averaged sigma_u and L (measured on the same flights) and the
    flight's own T_f, V_d and L_f. Nothing is fitted.
    """
    out = []
    for r in rows:
        if r["u_ref"] <= 0:
            continue
        su = r.get("mean__path_sigma_u")
        L = r.get("mean__path_L_int")
        if su is None or not np.isfinite(su) or not np.isfinite(L):
            continue
        Tf, Vd, Lf = r["full__T_f"], r["full__V_d"], r["full__L_f"]
        b = float(dl.sigma_ballistic(su, Tf, Vd, Lf, L, kappa=kappa))
        t = float(dl.sigma_tracer(su, Tf, Lf, L))
        m = r["gust__sigma"]
        out.append(dict(heading=r["heading"], throw=r["throw"], spin=r.get("spin"),
                        u_ref=r["u_ref"], N=Lf / L, sigma_u=su, L_int=L,
                        T_f=Tf, V_d=Vd, L_f=Lf,
                        measured=m, ballistic=b, tracer=t,
                        measured_se=r["gust__sigma_se"],
                        ratio_ballistic=b / max(m, 1e-12),
                        ratio_tracer=t / max(m, 1e-12)))
    return out


# ---------------------------------------------------------------------------
# Q3  the discriminating exponents
# ---------------------------------------------------------------------------

def shape_test(comp):
    """Regress log(measured) on log(predicted), for each law.

    THIS IS THE DISCRIMINATOR, and it is better than fitting sigma against L_f
    and V_d separately.

    Q2 already showed that the ballistic law's MAGNITUDE is wrong by a roughly
    constant factor. A constant factor is the least interesting kind of wrong:
    it is one number, it could be kappa, and it cannot distinguish a law that
    has the right physics from one that does not. What distinguishes them is
    whether the prediction tracks the measurement as the flight and the wind
    change -- that is, whether

        log(sigma_measured) = a + b * log(sigma_predicted)

    has b = 1. A law with the right shape and a wrong prefactor gives b = 1 with
    small scatter and a != 0. A law with the wrong shape gives b != 1 no matter
    what prefactor is chosen, because no constant can fix it.

    Fitting L_f and V_d exponents separately, as the brief's WP2 P2 proposes,
    is the same test done with fewer of the variables and it needs the two to
    be varied independently. This uses every cell in the sweep and every axis
    at once. Both are reported; they should agree.
    """
    out = {}
    if not comp:
        return out
    m = np.array([c["measured"] for c in comp])

    # THE TEST IS VACUOUS UNLESS THE FLIGHT ITSELF VARIES, and it will pass
    # anyway if it is not. Both laws are linear in sigma_u, so on a sweep where
    # only the wind speed moves, both predictions are proportional to sigma_u
    # and both slopes come out at 1 no matter which law is right. That is a
    # test of the wind-speed axis wearing a discriminator's clothes.
    #
    # What separates them lives in T_f, V_d and L_f. So the span of those is
    # reported with the verdict, and the verdict is downgraded to
    # "uninformative" unless the flight quantities move by a decent factor.
    span = {}
    for nm in ("L_f", "V_d", "T_f", "sigma_u"):
        v = np.array([c[nm] for c in comp])
        span[nm] = float(np.max(v) / max(np.min(v), 1e-12))
    # The criterion is not "does L_f move" but "do the two laws DISAGREE across
    # this set". Their ratio
    #
    #     ballistic/tracer = (kappa/2) (g L_f / V_d^2) sqrt(Phi(N) N)
    #
    # contains no sigma_u at all, so it is exactly the part of the prediction
    # that the wind-speed axis cannot exercise. If that ratio is near-constant
    # over the cells, no amount of data separates the laws and both will report
    # b = 1 -- which is what happens on a wind-speed-only sweep, where sigma_u
    # spans 42x and the ratio spans almost nothing. Requiring the ratio to move
    # by 1.6x is requiring the two predictions to be distinguishable by more
    # than the ~1.4x residual scatter these fits carry.
    bal = np.array([c["ballistic"] for c in comp])
    tra = np.array([c["tracer"] for c in comp])
    good = (bal > 0) & (tra > 0)
    lr = bal[good] / tra[good]
    ratio_span = float(np.max(lr) / max(np.min(lr), 1e-12)) if good.any() else 1.0
    informative = ratio_span >= 1.6
    out["_span"] = dict(span, law_ratio_span=ratio_span,
                        informative=bool(informative),
                        note=("the two laws' predictions stay within "
                              f"{ratio_span:.2f}x of each other across this "
                              "set, so both will fit with slope 1 whichever is "
                              "right: this sweep CANNOT discriminate them"
                              if not informative else
                              f"the two laws' predictions separate by "
                              f"{ratio_span:.2f}x across this set, enough to "
                              "discriminate"))

    for law in ("ballistic", "tracer"):
        p = np.array([c[law] for c in comp])
        ok = (m > 0) & (p > 0) & np.isfinite(m) & np.isfinite(p)
        if ok.sum() < 4:
            continue
        b, lo, hi, r2 = dl.fit_exponent(p[ok], m[ok])
        # The prefactor the fit would need, at b forced to 1: the median ratio.
        ratio = float(np.median(m[ok] / p[ok]))
        resid = np.log(m[ok]) - (np.log(p[ok]) + np.log(ratio))
        out[law] = dict(slope=b, ci=[lo, hi], r2=r2, n=int(ok.sum()),
                        prefactor_needed=ratio,
                        scatter_at_unit_slope=float(np.exp(np.std(resid))),
                        shape_ok=bool(lo <= 1.0 <= hi))
    return out


def exponent_tests(rows, u_ref=None):
    """sigma_gust against L_f and against V_d, across the throw catalogue.

    Held at ONE wind speed so sigma_u is common to every point; otherwise the
    fit would mix the wind-speed axis into the flight axis. The measured L_f
    and V_d are used, not the catalogue's nominal targets, because the wind
    changes both and the law is stated in terms of the flight that happened.

    Note that neither law's exponent is a pure power: the ballistic 1.7 is the
    LOCAL slope of L_f^2 sqrt(Phi(L_f/L)) at N ~ 4, and it drifts as N moves.
    The predicted value is therefore recomputed at each cell's own N rather
    than quoted as a constant.
    """
    res = {}
    for (h,), _ in group(rows, "heading").items():
        rs = [r for r in rows if r["heading"] == h]
        if u_ref is not None:
            rs = [r for r in rs if abs(r["u_ref"] - u_ref) < 1e-9]
        rs = [r for r in rs if r["u_ref"] > 0 and np.isfinite(r["gust__sigma"])]
        if len(rs) < 3:
            continue
        Lf = np.array([r["full__L_f"] for r in rs])
        Vd = np.array([r["full__V_d"] for r in rs])
        sg = np.array([r["gust__sigma"] for r in rs])
        Ns = np.array([r["full__L_f"] / r["mean__path_L_int"] for r in rs])
        pL, loL, hiL, r2L = dl.fit_exponent(Lf, sg)
        pV, loV, hiV, r2V = dl.fit_exponent(Vd, sg)
        # Predicted local ballistic L_f exponent at this set's own mean N.
        pred_L = dl.check_exponents_symbolically(float(np.mean(Ns)))["exponent_L_f_local"]
        res[int(h)] = dict(
            n_cells=len(rs), mean_N=float(np.mean(Ns)),
            L_f=dict(fitted=pL, ci=[loL, hiL], r2=r2L,
                     predicted_ballistic=pred_L, predicted_tracer=0.5),
            V_d=dict(fitted=pV, ci=[loV, hiV], r2=r2V,
                     predicted_ballistic=-3.0, predicted_tracer=-1.0),
            throws=sorted({r["throw"] for r in rs}))
    return res


# ---------------------------------------------------------------------------
# Q4  how much of the error is the wind
# ---------------------------------------------------------------------------

def error_budget(rows, r_target=3.0):
    """Split each cell's landing error into thrower, mean wind and gusts.

    `wind_fraction` is the share of the total mean-square landing error that
    would vanish in still air. It is the number a coach or a groundskeeper
    actually wants, and it is the one the literature's turbulence-intensity
    criteria cannot produce.

    `Pi` is the dimensionless playability index: total rms landing error in
    units of the catch radius. Pi < 1 means a typical throw is caught. The
    BAND BOUNDARIES ARE PROPOSALS, not measurements -- there is no established
    projectile-sport equivalent of the Lawson pedestrian-comfort bands, which
    is precisely the gap this is meant to fill, and someone has to write the
    first one down. They are stated here so they can be argued with.
    """
    out = []
    for r in rows:
        still = r["still__sigma"]
        b = r["wind__bias"]
        sw = r["wind__sigma"]
        sg = r["gust__sigma"]
        wind_ms = b ** 2 + sw ** 2 + sg ** 2
        total_ms = still ** 2 + wind_ms
        out.append(dict(
            heading=r["heading"], throw=r["throw"], u_ref=r["u_ref"],
            spin=r.get("spin"),
            sigma_still=still, bias_wind=b, sigma_windshear=sw, sigma_gust=sg,
            sigma_total=float(np.sqrt(total_ms)),
            wind_fraction=float(wind_ms / total_ms) if total_ms > 0 else 0.0,
            # The correctable variant: a thrower who knows the field aims off
            # the mean bias, so only the two random terms remain.
            wind_fraction_aimed=float((sw ** 2 + sg ** 2) / (still ** 2 + sw ** 2
                                                             + sg ** 2)),
            Pi=float(np.sqrt(total_ms) / r_target),
            U_over_Vd=float(r["field__U_mean"] / max(r["full__V_d"], 1e-9))))
    return out


def crossover_speed(rows, sigma_release=None):
    """Wind speed at which the wind terms equal the thrower's own scatter.

    WHY THIS AND NOT `wind_fraction` ALONE. Every absolute verdict in the error
    budget -- the wind fraction, Pi, the band -- is set by `sigma_release`, and
    `sigma_release` here is whatever `trajectories.RELEASE_DISTRIBUTIONS` says
    (speed +/- 1.5 m/s, elevation +/- 3 deg, nose +/- 4 deg, roll +/- 5 deg).
    Those numbers are labelled "practitioner-typical" in their own docstring and
    are NOT calibrated against a measured thrower. They produce a 7.4-7.8 m
    landing scatter on a 27 m throw, which is poor. A better thrower moves every
    absolute number in the budget.

    What does NOT depend on that calibration is the SHAPE. The wind terms all
    grow as Uref^1, so

        sigma_wind(Uref) = c * Uref     =>     U* = sigma_release / c

    and U* scales linearly with whatever release scatter you assume. Reporting
    `c` and `U*` together lets a reader put their own thrower in: halve the
    release scatter and U* halves.

    `sigma_release=None` uses each cell's own measured still-air scatter.
    """
    out = []
    for (h, th, sp), rs in group(rows, "heading", "throw", "spin").items():
        rs = [r for r in rs if r["u_ref"] > 0]
        if len(rs) < 3:
            continue
        u = np.array([r["u_ref"] for r in rs])
        w = np.array([np.sqrt(r["wind__bias"] ** 2 + r["wind__sigma"] ** 2
                              + r["gust__sigma"] ** 2) for r in rs])
        # Slope through the origin, which is what "all terms grow as Uref^1"
        # means; fitting an intercept as well would let a spurious offset
        # absorb the very linearity being claimed.
        c = float((u @ w) / (u @ u))
        resid = w - c * u
        rel_resid = float(np.sqrt(np.mean(resid ** 2)) / np.mean(w))
        sr = (sigma_release if sigma_release is not None
              else float(rs[0]["still__sigma"]))
        out.append(dict(heading=h, throw=th, spin=sp,
                        slope_m_per_ms=c, sigma_release=sr,
                        u_star=sr / c if c > 0 else float("inf"),
                        linear_fit_rel_residual=rel_resid))
    return out


PI_BANDS = [(0.0, 1.0, "routine"), (1.0, 2.0, "demanding"),
            (2.0, 4.0, "unreliable"), (4.0, np.inf, "not a throw")]


def pi_band(pi):
    for lo, hi, name in PI_BANDS:
        if lo <= pi < hi:
            return name
    return "not a throw"


# ---------------------------------------------------------------------------
# figures
# ---------------------------------------------------------------------------

def figures(rows, budget, comp, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    hcol = {90: "#2e8b57", 210: "#1f77b4", 270: "#d62728"}

    # --- fig 1: the three error contributions vs wind speed ----------------
    throws = sorted({r["throw"] for r in rows})
    # sharey, not just sharex. The point of the panel row is that the thrower's
    # own scatter is larger than either wind term until the wind is very
    # strong; three independently scaled panels would hide exactly that, by
    # drawing a 0.2 m gust effect at the same height as a 7.8 m release scatter.
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.4), sharex=True, sharey=True)
    keys = [("still__sigma", "thrower alone"),
            ("wind__bias", "mean-wind bias"),
            ("gust__sigma", "gust dispersion")]
    for ax, (key, title) in zip(axes, keys):
        for h in sorted({r["heading"] for r in rows}):
            for th in throws:
                rs = sorted([r for r in rows
                             if r["heading"] == h and r["throw"] == th],
                            key=lambda r: r["u_ref"])
                if not rs:
                    continue
                ax.plot([r["u_ref"] for r in rs], [r[key] for r in rs],
                        marker="o", ms=3, lw=1.4, color=hcol.get(h, "#555"),
                        ls="-" if th == throws[0] else "--",
                        label=f"{h}" if th == throws[0] else None)
        ax.set_title(title, fontsize=9)
        ax.set_xlabel(r"$U_{ref}$  (m/s)")
        ax.set_xlim(0, None)
        ax.grid(alpha=0.25, lw=0.5)
    # One shared, zero-based limit set from the largest value in ANY panel.
    # Letting each panel autoscale under sharey silently clips the others --
    # the gust panel reaches 16 m and the thrower panel 7.8, so whichever is
    # set last wins and the other loses its top.
    ymax = max(r[k] for r in rows for k, _ in keys)
    axes[0].set_ylim(0, 1.05 * ymax)
    axes[0].set_ylabel("landing error  (m)")
    axes[0].legend(frameon=False, fontsize=8, title="wind from", title_fontsize=8,
                   loc="lower left")
    # The line style carries the throw type and would otherwise be unlabelled.
    from matplotlib.lines import Line2D
    axes[2].legend(handles=[Line2D([], [], color="#555", ls=st, label=th)
                            for th, st in zip(throws, ("-", "--", ":", "-."))],
                   frameon=False, fontsize=8, loc="upper left")
    fig.tight_layout()
    fig.savefig(outdir / "fig_wind_speed.png", dpi=200, facecolor=BG)
    plt.close(fig)

    # --- fig 2: measured vs the two laws -----------------------------------
    if comp:
        fig, ax = plt.subplots(figsize=(4.4, 4.2))
        m = np.array([c["measured"] for c in comp])
        for key, lab, col in (("ballistic", "ballistic  " + r"$\Phi(N)$", "#d62728"),
                              ("tracer", "tracer  " + r"$1/\sqrt{N}$", "#1f77b4")):
            ax.scatter(m, [c[key] for c in comp], s=16, alpha=0.7,
                       color=col, label=lab)
        # Limits from the data on BOTH axes. The first version anchored the
        # lower end at 1e-3 m, which is three decades below anything plotted
        # and squashed every point into the top-right corner.
        allv = np.concatenate([m] + [np.array([c[k] for c in comp])
                                     for k in ("ballistic", "tracer")])
        allv = allv[np.isfinite(allv) & (allv > 0)]
        lim = [0.6 * float(allv.min()), 1.7 * float(allv.max())]
        ax.plot(lim, lim, color="#333", lw=1.0, ls="--")
        # Log axes: the quantity spans well over a decade across the sweep.
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(*lim)
        ax.set_ylim(*lim)
        ax.set_xlabel("measured gust dispersion  (m)")
        ax.set_ylabel("predicted  (m)")
        ax.legend(frameon=False, fontsize=8)
        ax.grid(alpha=0.25, lw=0.5, which="both")
        fig.tight_layout()
        fig.savefig(outdir / "fig_law_test.png", dpi=200, facecolor=BG)
        plt.close(fig)

    # --- fig 3: how much of the error is the wind --------------------------
    fig, ax = plt.subplots(figsize=(4.8, 3.6))
    for h in sorted({b["heading"] for b in budget}):
        for th in throws:
            bs = sorted([b for b in budget
                         if b["heading"] == h and b["throw"] == th],
                        key=lambda b: b["u_ref"])
            if not bs:
                continue
            ax.plot([b["u_ref"] for b in bs],
                    [100 * b["wind_fraction"] for b in bs],
                    marker="o", ms=3, lw=1.4, color=hcol.get(h, "#555"),
                    ls="-" if th == throws[0] else "--",
                    label=f"{h}" if th == throws[0] else None)
    ax.set_xlabel(r"$U_{ref}$  (m/s)")
    ax.set_ylabel("share of landing error due to wind  (%)")
    ax.set_xlim(0, None)
    ax.set_ylim(0, 100)
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(frameon=False, fontsize=8, title="wind dir", title_fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "fig_wind_fraction.png", dpi=200, facecolor=BG)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", default=str(HERE / "out" / "sweep_speed.json"))
    ap.add_argument("--r-target", type=float, default=3.0)
    ap.add_argument("--exponent-u-ref", type=float, default=5.0)
    ap.add_argument("--out", default=str(HERE / "out" / "theory_results.json"))
    ap.add_argument("--figdir", default=str(HERE / "figs"))
    args = ap.parse_args()

    meta, rows = load(args.sweep)
    print(f"{len(rows)} cells from {Path(args.sweep).name}")

    law = wind_speed_law(rows)
    print("\n[Q1] sigma_gust ~ Uref^p        (theory: p = 1 exactly)")
    for rec in law:
        e = rec["sigma_gust"]
        print(f"   heading {rec['heading']:>3}  {rec['throw']:<16} "
              f"p = {e['exponent']:.3f}  CI [{e['ci'][0]:.3f}, {e['ci'][1]:.3f}]"
              f"  r2 {e['r2']:.4f}   (still-air floor {rec['sigma_still']:.2f} m)")
    print("\n[Q1b] mean-wind bias ~ Uref^q   (theory: q = 1 if the flight is "
          "unchanged)")
    for rec in law:
        e = rec["bias_wind"]
        print(f"   heading {rec['heading']:>3}  {rec['throw']:<16} "
              f"q = {e['exponent']:.3f}  CI [{e['ci'][0]:.3f}, {e['ci'][1]:.3f}]"
              f"  r2 {e['r2']:.4f}")

    comp = law_comparison(rows)
    print("\n[Q2] magnitude: predicted / measured")
    if comp:
        rb = np.array([c["ratio_ballistic"] for c in comp])
        rt = np.array([c["ratio_tracer"] for c in comp])
        print(f"   ballistic  median {np.median(rb):.2f}x  "
              f"[{np.percentile(rb, 10):.2f}, {np.percentile(rb, 90):.2f}]")
        print(f"   tracer     median {np.median(rt):.2f}x  "
              f"[{np.percentile(rt, 10):.2f}, {np.percentile(rt, 90):.2f}]")
        print(f"   geometric mean of the two laws: "
              f"{np.median(np.sqrt(rb * rt)):.2f}x")

    shape = shape_test(comp)
    print("\n[Q3] shape: log(measured) = a + b log(predicted); b = 1 means "
          "the law tracks the\n     data as flight and wind change, whatever "
          "its prefactor")
    sp = shape.get("_span", {})
    if sp:
        print(f"   flight span over the sweep: L_f x{sp['L_f']:.2f}  "
              f"V_d x{sp['V_d']:.2f}  T_f x{sp['T_f']:.2f}  "
              f"sigma_u x{sp['sigma_u']:.2f}   |  the two laws separate by "
              f"x{sp['law_ratio_span']:.2f}")
        if not sp["informative"]:
            print(f"   !! {sp['note'].upper()}")
    for law, d in shape.items():
        if law == "_span":
            continue
        print(f"   {law:<10} b = {d['slope']:.3f}  CI [{d['ci'][0]:.3f}, "
              f"{d['ci'][1]:.3f}]  r2 {d['r2']:.4f}  n {d['n']}   "
              f"prefactor needed {d['prefactor_needed']:.3f}x   "
              f"residual scatter {d['scatter_at_unit_slope']:.3f}x   "
              f"{'SHAPE OK' if d['shape_ok'] else 'SHAPE REJECTED'}")

    ex = exponent_tests(rows, u_ref=args.exponent_u_ref)
    print(f"\n[Q3] discriminating exponents at Uref = {args.exponent_u_ref} m/s")
    for h, d in ex.items():
        print(f"   heading {h}  ({d['n_cells']} cells, "
              f"{len(d['throws'])} throw types, mean N = {d['mean_N']:.2f})")
        for nm in ("L_f", "V_d"):
            q = d[nm]
            print(f"     sigma ~ {nm}^p:  fitted {q['fitted']:+.2f} "
                  f"CI [{q['ci'][0]:+.2f}, {q['ci'][1]:+.2f}]   "
                  f"ballistic {q['predicted_ballistic']:+.2f}   "
                  f"tracer {q['predicted_tracer']:+.2f}")
    if not ex:
        print("   (needs >= 3 throw types in the sweep; run with the full "
              "catalogue)")

    cross = crossover_speed(rows)
    print("\n[Q4a] wind vs thrower: sigma_wind = c*Uref, and where it "
          "overtakes the release scatter")
    print(f"   {'head':>5} {'throw':<16} {'c (m per m/s)':>14} "
          f"{'sigma_release':>13} {'U* (m/s)':>9} {'lin resid':>10}")
    for c in cross:
        print(f"   {c['heading']:5d} {c['throw']:<16} "
              f"{c['slope_m_per_ms']:14.4f} {c['sigma_release']:13.2f} "
              f"{c['u_star']:9.1f} {100*c['linear_fit_rel_residual']:9.1f}%")
    print("   U* scales linearly with the assumed release scatter: a thrower "
          "with half this\n   scatter reaches U* at half the wind speed. The "
          "release distribution is NOT\n   calibrated against a measured "
          "player -- see crossover_speed's docstring.")

    budget = error_budget(rows, args.r_target)
    print(f"\n[Q4] error budget at r_target = {args.r_target} m")
    print(f"   {'Uref':>5} {'head':>5} {'throw':<16} {'thrower':>8} "
          f"{'bias':>7} {'shear':>7} {'gust':>7} {'total':>7} {'wind%':>6} "
          f"{'Pi':>6}  band")
    for b in sorted(budget, key=lambda b: (b["throw"], b["heading"], b["u_ref"])):
        print(f"   {b['u_ref']:5.1f} {b['heading']:5d} {b['throw']:<16} "
              f"{b['sigma_still']:8.2f} {b['bias_wind']:7.2f} "
              f"{b['sigma_windshear']:7.2f} {b['sigma_gust']:7.2f} "
              f"{b['sigma_total']:7.2f} {100*b['wind_fraction']:6.1f} "
              f"{b['Pi']:6.2f}  {pi_band(b['Pi'])}")

    figures(rows, budget, comp, args.figdir)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(dict(meta=meta, wind_speed_law=law,
                                   law_comparison=comp, shape_test=shape,
                                   exponent_tests=ex,
                                   crossover=cross, error_budget=budget,
                                   r_target=args.r_target), indent=2))
    print(f"\nwrote {out} and figures in {args.figdir}")


if __name__ == "__main__":
    main()
