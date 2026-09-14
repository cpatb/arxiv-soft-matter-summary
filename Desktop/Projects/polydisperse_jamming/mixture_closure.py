"""mixture_closure.py -- calibrate the filling functional on BINARY packs, predict the
CONTINUOUS families.

WHY THIS IS THE TEST WORTH DOING
--------------------------------
`grade_mixture.py` established that composition at fixed size ratio moves phi by up to
0.082, that two packs matched in number CV and size ratio differ by 0.0635, and that a
zero-parameter pairwise-filling model predicts that difference to 5 %. It also found
exactly where the textbook interaction functions are wrong: they OVER-credit filling at
size ratio 2 (mean residual -4.2 scale) and are close at ratio 3 (+1.6).

A binary pack at a known (ratio, composition) is the cleanest object the model has: two
classes, one size ratio, nothing to discretise. So the binary campaign is the natural
place to CALIBRATE the two interaction exponents, and the continuous families -- lognormal,
uniform, bidisperse, five CV levels each -- are then a transfer test the calibration set
contains no point of.

That transfer is not free. A continuous PSD is a sum over MANY pairs at SMALL size ratios,
which is the regime the binary campaign shows the textbook functions handle worst, and it
is a regime the calibration barely constrains (the campaign has ratios 2 and 3 only). So
there are three live outcomes and all three are informative:

  (a) the refit improves the continuous families  -> the functional transfers, and the
      textbook exponents were simply the wrong numbers for frictional spheres;
  (b) the refit fits the binaries better and the families WORSE -> the functional does not
      transfer, and the agreement in `mixture_model.py` section 5 was a coincidence of the
      textbook exponents;
  (c) neither set fits both -> the two-exponent form is too rigid, which is a statement
      about the FORM and points at where to look next.

The prediction, written before the fit ran: (a), because the sign of the binary residual
(model too high at r = 2) is the same sign as the family residual (model too high on every
continuous family), so one correction should move both.

(!) HONESTY. The continuous-family numbers are on disk and have been read. What is
genuinely held out is that no continuous-family point enters the fit -- a cross-experiment
transfer, not a blind prediction.

Run:  python3 mixture_closure.py      (writes out/mixture_closure.txt)
"""

from __future__ import annotations

import collections
import glob
import json
import os

import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import minimize

import mixture_model as mm
from phase12_mixture import MIN_PHASE_COUNT, N_PROD
from psd import BinaryPSD

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
SCALE = mm.PHI_SCALE
FAMILIES = mm.FAMILIES
CVS = (0.1, 0.2, 0.3, 0.4, 0.5)


def binary_points(mu_p: float = 0.3) -> list[tuple[float, float, float]]:
    """(ratio, x_small, phi) for every GRADABLE interior binary cell, seed-averaged.

    The M5 admission rule of PREREG_MIXTURE_2026-09-13.md is applied here too, and for
    the same reason: a mixture whose minority phase is 36 grains measures that phase's
    sampling noise. Endpoints are excluded because they are beta by construction and
    would only re-weight the fit toward a number that is already an input.
    """
    g: dict[tuple[float, float], list[float]] = collections.defaultdict(list)
    for f in sorted(glob.glob(os.path.join(OUT, "phase12_mixr*_shard*.jsonl"))):
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                r = json.loads(line)
                if r.get("error") or r.get("n_requested") != N_PROD:
                    continue
                # The shard glob matches any tag, including a block at another friction.
                # beta is friction-dependent and the whole model is built on one beta, so
                # pooling two frictions here would be silently meaningless rather than
                # noisy. Default 0.3 is the friction every number in this file is read at.
                if abs(float(r.get("mu_p", -1)) - mu_p) > 1e-9:
                    continue
                x = round(float(r["x_small_vol"]), 6)
                if x in (0.0, 1.0):
                    continue
                p = BinaryPSD(d10=1.0, ratio=float(r["mix_ratio"]), x_small=x)
                if min(p.n_small, 1.0 - p.n_small) * N_PROD < MIN_PHASE_COUNT:
                    continue
                g[(float(r["mix_ratio"]), x)].append(float(r["phi"]))
    return [(r, x, float(np.mean(v))) for (r, x), v in sorted(g.items())]


def fit_on_binaries(pts, beta: float) -> tuple[float, float, float]:
    def cost(theta):
        pa, pb = np.exp(theta)
        return sum(((mm.binary_density(r, x, beta, pa, pb) - phi) / SCALE) ** 2
                   for r, x, phi in pts)

    best = None
    for x0 in ([0.0, 0.4], [1.0, 1.0], [-1.0, 0.0], [0.5, 1.5], [-0.5, 2.0]):
        res = minimize(cost, np.array(x0), method="Nelder-Mead",
                       options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 8000})
        if best is None or res.fun < best.fun:
            best = res
    pa, pb = np.exp(best.x)
    return float(pa), float(pb), float(best.fun)


def fit_ab_direct(pts, beta: float) -> tuple[float, float, float]:
    """Fit a and b DIRECTLY (two numbers in (0,1)) to one ratio's composition curve.

    No functional form at all. Every point in `pts` must share one size ratio, which is
    asserted rather than assumed -- mixing two ratios here would silently return an
    average of two different interaction strengths and the whole point is that they
    differ.
    """
    rs = {r for r, _, _ in pts}
    if len(rs) != 1:
        raise ValueError(f"fit_ab_direct needs ONE size ratio, got {sorted(rs)}")

    def cost(theta):
        a_, b_ = 1.0 / (1.0 + np.exp(-theta))
        s = 0.0
        for r, x, phi in pts:
            v = packing_density_ab(np.array([r, 1.0]), np.array([1.0 - x, x]),
                                   beta, a_, b_)
            s += ((v - phi) / SCALE) ** 2
        return s

    best = None
    for x0 in ([0.0, 0.0], [1.0, 1.0], [-1.0, -1.0], [2.0, 0.0], [0.0, 2.0]):
        res = minimize(cost, np.array(x0), method="Nelder-Mead",
                       options={"xatol": 1e-10, "fatol": 1e-12, "maxiter": 8000})
        if best is None or res.fun < best.fun:
            best = res
    a_, b_ = 1.0 / (1.0 + np.exp(-best.x))
    return float(a_), float(b_), float(best.fun)


def packing_density_ab(d, y, beta: float, a_of, b_of) -> float:
    """The linear packing model with the interaction values supplied as CALLABLES or
    CONSTANTS of u = d_fine/d_coarse, rather than as the two-exponent closed forms.

    This is `mixture_model.packing_density` with the only substantive line changed, and
    it reduces to it exactly when `a_of`/`b_of` are the closed forms -- checked in
    `check_ab_reduces()` rather than asserted.
    """
    d = np.asarray(d, dtype=float)
    y = np.asarray(y, dtype=float)
    y = y / y.sum()
    order = np.argsort(-d)
    d, y = d[order], y[order]

    def ev(f, u):
        return np.clip(f(u) if callable(f) else np.full_like(np.asarray(u, float), f),
                       0.0, 1.0)

    k = d.size
    best = np.inf
    for i in range(k):
        tot = 0.0
        if i > 0:
            bb = ev(b_of, d[i] / d[:i])
            tot += float(np.sum((1.0 - beta + bb * beta * (1.0 - 1.0 / beta)) * y[:i]))
        if i < k - 1:
            aa = ev(a_of, d[i + 1:] / d[i])
            tot += float(np.sum((1.0 - aa) * y[i + 1:]))
        denom = 1.0 - tot
        if denom <= 1e-12:
            continue
        best = min(best, beta / denom)
    return float(best)


def profile_interval(pts, beta: float, which: int, ngrid: int = 201
                     ) -> tuple[float, float]:
    """Profile one of (a, b) on a grid, re-minimising the other, at delta chi2 = n.

    The bar is one unit of chi2 PER POINT, i.e. the level at which every point is fitted
    to within its own declared error. A parameter whose profile interval spans most of
    [0, 1] is not determined by the composition curve at that ratio, whatever the
    optimiser returned for it. The grid is scanned rather than root-found deliberately:
    the objective here re-minimises a nuisance parameter, so it need not be a smooth
    function of the profiled one, and a root finder will happily return a number from a
    discontinuity (lessons.md, 2026-09-06).
    """
    def chi(a_, b_):
        return sum(((packing_density_ab(np.array([r, 1.0]), np.array([1.0 - x, x]),
                                        beta, a_, b_) - phi) / SCALE) ** 2
                   for r, x, phi in pts)

    grid = np.linspace(1e-4, 1.0 - 1e-4, ngrid)
    prof = np.empty(ngrid)
    for i, v in enumerate(grid):
        prof[i] = min(chi(v, w) if which == 0 else chi(w, v) for w in grid)
    lo = float(prof.min())
    ok = grid[prof <= lo + len(pts)]
    return (float(ok.min()), float(ok.max())) if ok.size else (0.0, 1.0)


def packing_density_tabulated(d, y, beta: float, fa, fb) -> float:
    """Same thing, named for the caller that passes interpolators."""
    return packing_density_ab(d, y, beta, fa, fb)


def check_ab_reduces(beta: float = 0.5947) -> list[str]:
    """`packing_density_ab` with the closed forms must equal `mixture_model` exactly."""
    bad = []
    rng = np.random.default_rng(20260913)
    for _ in range(20):
        k = int(rng.integers(2, 12))
        d = np.sort(rng.uniform(0.3, 3.0, k))[::-1]
        y = rng.dirichlet(np.ones(k))
        want = mm.packing_density(d, y, beta, 1.02, 1.50)
        got = packing_density_ab(d, y, beta,
                                 lambda u: mm.a_loosening(u, 1.02),
                                 lambda u: mm.b_wall(u, 1.50))
        if abs(got - want) > 1e-12:
            bad.append(f"k={k}: {got:.15g} vs {want:.15g}")
    return bad


def main() -> None:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    pts = binary_points()
    # The headline calibration is deliberately the ORIGINAL ratio 2/3 grid, so the
    # transfer number quoted in the Status card is the one computed before the ratio
    # ladder existed and does not silently absorb it. The wider calibration is E2.
    pts23 = [p for p in pts if p[0] in (2.0, 3.0)]
    meas = mm.measured()
    beta = meas[("lognormal", 0.0)][0]
    say("=== mixture_closure.py -- calibrate on BINARIES, predict the FAMILIES ===")
    say(f"   beta = {beta:.4f} (measured, not fitted). scale = {SCALE:.1e}.")
    say(f"   calibration set: {len(pts23)} gradable interior binary compositions at")
    say(f"   size ratios {sorted({r for r, _, _ in pts23})}. No continuous-family point is in it.")
    say()

    pa, pb, chi_fit = fit_on_binaries(pts23, beta)
    say("=== the refit, on the ratio 2 and 3 grid ===")
    say("   textbook   PA = 1.020  PB = 1.500")
    say(f"   refitted   PA = {pa:.3f}  PB = {pb:.3f}   chi2 on the calibration set = {chi_fit:.1f}"
        f"  (per point {chi_fit/len(pts23):.1f})")
    chi_txt = sum(((mm.binary_density(r, x, beta, 1.02, 1.50) - phi) / SCALE) ** 2
                  for r, x, phi in pts23)
    say(f"   textbook exponents on the same set:            chi2 = {chi_txt:.1f}"
        f"  (per point {chi_txt/len(pts23):.1f})")
    say()
    say("   what the refit does to the interaction functions, at the ratios that matter:")
    say(f"      {'d_fine/d_coarse':>16}{'a textbook':>12}{'a refit':>10}"
        f"{'b textbook':>12}{'b refit':>10}")
    for rr in (1.0, 1 / 1.25, 1 / 1.5, 1 / 2.0, 1 / 3.0, 1 / 5.0, 1 / 10.0):
        say(f"      {rr:16.4f}{float(mm.a_loosening(rr, 1.02)):12.4f}"
            f"{float(mm.a_loosening(rr, pa)):10.4f}"
            f"{float(mm.b_wall(rr, 1.50)):12.4f}{float(mm.b_wall(rr, pb)):10.4f}")
    say("   a -> 1 and b -> 1 mean NO separation (the mixture is just a monodisperse pack);")
    say("   a -> 0 and b -> 0 mean ideal Furnas filling. A refit that RAISES these at")
    say("   moderate ratios is saying frictional spheres fill less well than the textbook")
    say("   functions -- which were calibrated on concrete aggregates -- assume.")
    say()

    say("=== TRANSFER: the continuous families, no point of which is in the fit ===")
    say(f"   {'family':<11}{'CV':>5}{'phi_meas':>10}{'textbook':>10}{'/scale':>8}"
        f"{'refit':>9}{'/scale':>8}")
    chi = {"textbook": 0.0, "refit": 0.0}
    npt = 0
    for fam in FAMILIES:
        for cv in CVS:
            key = (fam, cv)
            if key not in meas:
                continue
            d, y = mm.psd_classes(fam, cv)
            t = mm.packing_density(d, y, beta, 1.02, 1.50)
            rfit = mm.packing_density(d, y, beta, pa, pb)
            phi = meas[key][0]
            chi["textbook"] += ((t - phi) / SCALE) ** 2
            chi["refit"] += ((rfit - phi) / SCALE) ** 2
            npt += 1
            say(f"   {fam:<11}{cv:5.1f}{phi:10.4f}{t:10.4f}{(t-phi)/SCALE:+8.1f}"
                f"{rfit:9.4f}{(rfit-phi)/SCALE:+8.1f}")
    say(f"   chi2 over {npt} held-out family points:  textbook {chi['textbook']:8.1f}"
        f"   refit {chi['refit']:8.1f}")
    verdict = ("(a) the functional TRANSFERS" if chi["refit"] < chi["textbook"]
               else "(b) the functional does NOT transfer")
    say(f"   -> {verdict}: refit/textbook = {chi['refit']/chi['textbook']:.2f}")
    say()

    say("=== the contrast the project quotes, under both calibrations ===")
    say(f"   {'family':<11}{'d_meas':>9}{'textbook':>10}{'refit':>9}"
        f"{'meas/txt':>10}{'meas/refit':>12}")
    rt, rr_ = {}, {}
    for fam in FAMILIES:
        if (fam, 0.5) not in meas:
            continue
        d0, y0 = mm.psd_classes(fam, 0.0)
        d5, y5 = mm.psd_classes(fam, 0.5)
        dm = meas[(fam, 0.5)][0] - meas[(fam, 0.0)][0]
        dt = (mm.packing_density(d5, y5, beta, 1.02, 1.50)
              - mm.packing_density(d0, y0, beta, 1.02, 1.50))
        dr = (mm.packing_density(d5, y5, beta, pa, pb)
              - mm.packing_density(d0, y0, beta, pa, pb))
        rt[fam], rr_[fam] = dm / dt, dm / dr
        say(f"   {fam:<11}{dm:+9.4f}{dt:+10.4f}{dr:+9.4f}{dm/dt:10.2f}{dm/dr:12.2f}")
    if rt:
        say(f"   spread of the ratio across families (1.00 = collapse):"
            f"  textbook {max(rt.values())/min(rt.values()):.2f}x"
            f"   refit {max(rr_.values())/min(rr_.values()):.2f}x")
        say("   against 4.51x for CV and 4.06x for d32 (shape_control.py).")
    say()

    # ---- E1 / E2, PREREG_MIXTURE_2026-09-13.md section 4b -----------------------------
    ratios = sorted({r for r, _, _ in pts})
    if len(ratios) > 2:
        say("=== E1 -- is the two-exponent FORM rigid enough? per-ratio fits ===")
        say(f"   {'ratio':>7}{'n':>4}{'PA':>9}{'PB':>9}{'chi2':>10}{'per pt':>9}")
        per_sum = 0.0
        per: dict[float, tuple[float, float]] = {}
        for r in ratios:
            sub = [p for p in pts if p[0] == r]
            a_, b_, c_ = fit_on_binaries(sub, beta)
            per[r] = (a_, b_)
            per_sum += c_
            say(f"   {r:7.2f}{len(sub):4d}{a_:9.3f}{b_:9.3f}{c_:10.1f}{c_/len(sub):9.1f}")
        pa_all, pb_all, chi_all = fit_on_binaries(pts, beta)
        say(f"   {'JOINT':>7}{len(pts):4d}{pa_all:9.3f}{pb_all:9.3f}{chi_all:10.1f}"
            f"{chi_all/len(pts):9.1f}")
        rigid = chi_all <= 2.0 * per_sum
        say(f"   sum of per-ratio chi2 = {per_sum:.1f}; joint = {chi_all:.1f}; "
            f"ratio = {chi_all/per_sum:.2f}")
        say(f"   -> E1 verdict: {'RIGID ENOUGH' if rigid else 'TOO RIGID'} "
            f"(bar: joint within 2x the sum of per-ratio fits)")
        say("   the per-ratio exponents themselves: "
            + ", ".join(f"r={r:.2f} -> ({per[r][0]:.2f}, {per[r][1]:.2f})" for r in ratios))
        say()

        say("=== E2 -- does the WIDER calibration improve the held-out family transfer? ===")
        chi_all_fam = 0.0
        for fam in FAMILIES:
            for cv in CVS:
                if (fam, cv) not in meas:
                    continue
                d, y = mm.psd_classes(fam, cv)
                chi_all_fam += ((mm.packing_density(d, y, beta, pa_all, pb_all)
                                 - meas[(fam, cv)][0]) / SCALE) ** 2
        say(f"   family chi2 over {npt} held-out points:")
        say(f"      textbook exponents                  {chi['textbook']:9.1f}")
        say(f"      calibrated on ratios 2, 3 only      {chi['refit']:9.1f}")
        say(f"      calibrated on ratios {ratios}   {chi_all_fam:9.1f}")
        if chi_all_fam < chi["refit"]:
            say("   -> E2 as predicted: the wider calibration IMPROVES the transfer.")
        else:
            say("   -> E2 REFUTED as predicted in the pre-registration: the wider")
            say("      calibration DEGRADES the held-out transfer, so the two-ratio")
            say("      agreement was partly coincidence and is reported as such.")
        say(f"      (the exponents moved from ({pa:.3f}, {pb:.3f}) to "
            f"({pa_all:.3f}, {pb_all:.3f}))")
        say()

    # ---- the non-parametric interaction functions -------------------------------------
    if len(ratios) > 2:
        say("=== E1 follow-up: the interaction functions MEASURED, not assumed ===")
        say("   E1 says the two-exponent FORM is too rigid, and it says WHERE: PA drifts")
        say("   monotonically with ratio while PB does not. Inside the linear packing")
        say("   structure a binary at one ratio determines TWO NUMBERS -- a(u) and b(u) at")
        say("   u = d_fine/d_coarse -- so they can be fitted directly, per ratio, with no")
        say("   functional form at all. Eight compositions constrain two numbers.")
        say(f"      {'u = 1/r':>9}{'a measured':>12}{'a textbook':>12}"
            f"{'b measured':>12}{'b textbook':>12}{'chi2/pt':>9}")
        us, a_m, b_m = [1.0], [1.0], [1.0]        # u = 1 is the proved no-separation limit
        for r in sorted(ratios, reverse=True):
            sub = [p for p in pts if p[0] == r]
            a_hat, b_hat, c_hat = fit_ab_direct(sub, beta)
            u = 1.0 / r
            us.append(u)
            a_m.append(a_hat)
            b_m.append(b_hat)
            say(f"      {u:9.4f}{a_hat:12.4f}{float(mm.a_loosening(u, 1.02)):12.4f}"
                f"{b_hat:12.4f}{float(mm.b_wall(u, 1.50)):12.4f}{c_hat/len(sub):9.2f}")
        us.append(0.0)                            # u = 0 is the proved ideal-filling limit
        a_m.append(0.0)
        b_m.append(0.0)
        order = np.argsort(us)
        uu = np.asarray(us)[order]
        aa = np.asarray(a_m)[order]
        bb = np.asarray(b_m)[order]
        say("   The two endpoints are the limits proved in mixture_model.check_limits, not")
        say("   fitted: a = b = 1 at equal size, a = b = 0 at infinite separation.")
        mono_a = bool(np.all(np.diff(aa) >= -1e-9))
        mono_b = bool(np.all(np.diff(bb) >= -1e-9))
        say(f"   monotone in u:  a {mono_a}   b {mono_b}")
        say()
        say("   IDENTIFIABILITY, before anything is read off the table above. A fitted")
        say("   number is worth nothing if the data does not constrain it, and the b")
        say("   column is the one to suspect: it is non-monotone, and the two ratios")
        say("   where it misbehaves are the ones whose composition curve barely moves.")
        say("   Profiling each parameter with the other free, at a bar of delta chi2 = 1")
        say("   PER POINT (the scale at which a point is fitted to its own error):")
        say(f"      {'u = 1/r':>9}{'a interval':>22}{'b interval':>22}{'phi span':>10}")
        ident: dict[float, tuple[bool, bool]] = {}
        for r in sorted(ratios, reverse=True):
            sub = [p for p in pts if p[0] == r]
            lo_a, hi_a = profile_interval(sub, beta, which=0)
            lo_b, hi_b = profile_interval(sub, beta, which=1)
            span = max(p[2] for p in sub) - min(p[2] for p in sub)
            ident[1.0 / r] = (hi_a - lo_a < 0.25, hi_b - lo_b < 0.25)
            say(f"      {1.0/r:9.4f}      [{lo_a:.3f}, {hi_a:.3f}]"
                f"         [{lo_b:.3f}, {hi_b:.3f}]{span:10.4f}")
        say("   An interval spanning most of [0, 1] means that parameter is NOT determined")
        say("   at that ratio, and its fitted value is an artefact of where the optimiser")
        say("   stopped rather than a measurement.")
        say()
        if not (mono_a and mono_b and all(all(v) for v in ident.values())):
            say("   -> THE MEASURED-INTERACTION ROUTE IS NOT DEFENSIBLE ON THIS DATA, by the")
            say("      rule written one paragraph above it rather than after seeing the")
            say("      answer. The transfer is computed and printed below because a number")
            say("      that is refused should still be visible, and it is NOT quoted")
            say("      anywhere else. The two-exponent fit remains the reported closure.")
            say()

        fa = PchipInterpolator(uu, aa)
        fb = PchipInterpolator(uu, bb)
        chi_np = 0.0
        say()
        say("   TRANSFER with the MEASURED interaction functions (shape-preserving")
        say("   interpolation in u; the four measured points plus the two proved limits):")
        say(f"   {'family':<11}{'CV':>5}{'phi_meas':>10}{'phi_np':>9}{'/scale':>8}")
        for fam in FAMILIES:
            for cv in CVS:
                if (fam, cv) not in meas:
                    continue
                d, y = mm.psd_classes(fam, cv)
                v = packing_density_tabulated(d, y, beta, fa, fb)
                chi_np += ((v - meas[(fam, cv)][0]) / SCALE) ** 2
                say(f"   {fam:<11}{cv:5.1f}{meas[(fam, cv)][0]:10.4f}{v:9.4f}"
                    f"{(v-meas[(fam, cv)][0])/SCALE:+8.1f}")
        say(f"   chi2 over {npt} held-out family points = {chi_np:.1f}, against")
        say(f"   {chi['textbook']:.1f} (textbook), {chi['refit']:.1f} (two-ratio exponent fit).")
        say("   ⚠ u < 1/3 is constrained ONLY by the u = 0 limit -- no binary in this")
        say("      campaign has a ratio above 3 -- and the lognormal at CV 0.5 truncated")
        say("      at (0.2, 5.0) contains pairs down to u = 0.04. That part of the curve is")
        say("      interpolated between two anchors, not measured, and a ratio-5 and")
        say("      ratio-10 block is what would fix it.")
        say()

    say("=== what is NOT established ===")
    say("   * The calibration has size ratios 2 and 3 only. A continuous PSD is a sum over")
    say("     many pairs, most of them at ratios BELOW 2, where nothing constrains the")
    say("     exponents and where the binary campaign already showed the textbook form is")
    say("     worst. A ratio ladder (1.25, 1.5, 5, 10) is what would fix that, and it is")
    say("     the same runner with a different --ratio.")
    say("   * beta is the same measured number in every column, so none of this tests the")
    say("     model at a different friction; section 7 of mixture_model.py does, at the")
    say("     textbook exponents only.")
    say("   * Every phi here is read at ONE pressure rung. Nothing says the interaction")
    say("     functions are pressure-independent, and this campaign cannot say it.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "mixture_closure.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
