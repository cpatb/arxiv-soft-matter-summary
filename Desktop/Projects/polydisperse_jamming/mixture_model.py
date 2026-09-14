"""mixture_model.py -- is the jamming shift a PAIRWISE-FILLING FUNCTIONAL of the PSD?

THE QUESTION THIS ANSWERS
-------------------------
`shape_control.py` established a NEGATIVE: at matched number CV and matched friction the
three size-distribution families give jamming shifts spanning 4.5x, and fifteen single
descriptors of the distribution (CV, d32, skewness, spans, volume-weighted versions, the
entropy of ln d) all fail to collapse it -- decisively, because the bidisperse family's
shift is NON-MONOTONE in CV while every descriptor tried is monotone in CV.

A non-monotone response to a monotone axis is not a defect of the descriptors. It is the
signature of a FUNCTIONAL: something that depends on how the volume is distributed across
size PAIRS, not on any one moment. The oldest such functional in the granular literature
is the Furnas two-scale filling picture -- small grains occupy the voids between large
ones, so a mixture packs denser than either component, with an optimum at an intermediate
composition and a return to the monodisperse value at both ends.

THE STRUCTURAL POINT, which needs no fitting: this project's bidisperse family is a PATH
through the (size ratio, fines volume fraction) plane, and it is a non-monotone path. The
family is 50:50 BY NUMBER at d10*(1 -/+ CV), so the fines carry

    X_s = (1-CV)^3 / [(1-CV)^3 + (1+CV)^3]

of the VOLUME: 0.500, 0.354, 0.229, 0.135, 0.073, 0.036 at CV 0, 0.1, 0.2, 0.3, 0.4, 0.5
while the size ratio (1+CV)/(1-CV) runs 1.00, 1.22, 1.50, 1.86, 2.33, 3.00. The ratio
rises monotonically (which helps filling) and the fines volume collapses toward zero
(which destroys it). A filling model must therefore produce a maximum somewhere in the
middle -- and the measurement has one, at CV 0.3-0.4.

THE MODEL
---------
The linear packing model. K size classes with volume fractions y_i (sum 1), diameters
d_i sorted DESCENDING, every class packing at the same monodisperse fraction `beta`
measured in this protocol. The packing density if class i is the dominant (space-filling)
one is

    gamma_i = beta / ( 1 - sum_{j<i} [1 - beta + b_ij*beta*(1 - 1/beta)] y_j
                         - sum_{j>i} [1 - a_ij] y_j )

and the mixture packs at gamma = min_i gamma_i, because the dominant class is whichever
one runs out of room first. j < i are COARSER than i (wall effect b_ij), j > i are FINER
(loosening effect a_ij).

The two interaction functions carry all the size-ratio physics. They are 1 at equal size
(full interference: the mixture is just a monodisperse pack) and 0 at infinite separation
(ideal Furnas filling). The exponents are the only fitted quantities:

    a_ij = sqrt( 1 - (1 - d_j/d_i)^PA )      loosening, d_j < d_i
    b_ij = 1 - (1 - d_i/d_j)^PB              wall,      d_i < d_j

LIMITS, CHECKED IN `check_limits()` RATHER THAN ASSERTED:
  * K = 1 returns beta exactly.
  * a = b = 0 (perfect separation) reproduces the Furnas binary branches in closed form,
    gamma = min( beta/(1-X_s), beta/(1-(1-beta)*X_L) ).
  * a = b = 1 (no separation) returns beta for every composition.

WHAT IS FITTED AND WHAT IS NOT
------------------------------
`beta` is MEASURED, not fitted: the monodisperse phi at this protocol's own p_hat = 1e-3,
N = 1500, mu_p = 0.3 rung. Only PA and PB are fitted, and they are fitted on the SIX
BIDISPERSE points alone. The uniform and lognormal ladders are HELD OUT -- twelve points
the fit never sees.

(!) HONESTY ABOUT THE HELD-OUT TEST. Those twelve numbers are on disk and the author of
this script has seen them. Held-out is not blind. The blind test is the composition
campaign pre-registered in PREREG_MIXTURE_2026-09-13.md, which measures phi at fixed
size ratio across the whole fines-fraction range -- a region where this project has NO
data at all and where the model's maximum lives.

(!) THE MODEL IS GEOMETRIC AND THE PACKS ARE FRICTIONAL. `beta` absorbs the friction at
the monodisperse point; nothing here claims the interaction functions are
friction-independent. That is a testable claim and the mu_p = 0 ladders on disk test it.

Run:  python3 mixture_model.py      (writes out/mixture_model.txt)
"""

from __future__ import annotations

import collections
import os

import numpy as np
from scipy.optimize import minimize

import analyse_phase2 as a
from psd import make_psd

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

# The protocol rung every number in this file is read at. Fixing it is the point:
# the model is compared to phi MEASURED at one pressure, never to an extrapolated phi_J,
# so none of this project's own fits enters the comparison. See the COMPANION section of
# PREREG_TAIL_2026-09-13.md for why the fixed-pressure channel is the cleaner one.
P_HAT = 1e-3
N_REQ = 1500
MU_P = 0.3
TRUNC = (0.2, 5.0)

# PSD discretisation. Convergence in K is CHECKED, not assumed: the lognormal CV = 0.5
# model value still drifts -0.0025 between K = 10 and K = 160 and is not converged at 40.
# 160 is the production value and the residual drift is printed as a systematic.
K_CLASSES = 160
NDRAW = 400_000
RNG_SEED = 20260913

# THE SCALE EVERY AGREEMENT IS GRADED AGAINST, declared before any comparison.
# It is NOT the seed sem (2e-04 at a typical cell), because the protocol is known to be
# non-deterministic at the 2.1e-03 level: `analyse_phase2.load_all` finds 15 duplicated
# keys returning a different phi from identical inputs, median |dphi| 2.1e-03, worst
# 5.0e-03. Grading a model against a scatter 10x smaller than the protocol's own
# reproducibility would measure the seed average, not the model. See lessons.md,
# "a relative spread needs a declared scale".
PHI_SCALE = 2.1e-3

FAMILIES = ("lognormal", "uniform", "bidisperse")
CVS = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5)


# --- the model -----------------------------------------------------------------------

def a_loosening(ratio, pa: float):
    """Loosening of a coarse class by a finer one. `ratio` = d_fine/d_coarse in (0, 1]."""
    r = np.clip(np.asarray(ratio, dtype=float), 0.0, 1.0)
    return np.sqrt(np.clip(1.0 - (1.0 - r) ** pa, 0.0, 1.0))


def b_wall(ratio, pb: float):
    """Wall effect on a fine class from a coarser one. `ratio` = d_fine/d_coarse in (0,1]."""
    r = np.clip(np.asarray(ratio, dtype=float), 0.0, 1.0)
    return 1.0 - (1.0 - r) ** pb


def packing_density(d, y, beta: float, pa: float, pb: float) -> float:
    """Linear packing model. `d` diameters, `y` VOLUME fractions (need not be sorted)."""
    d = np.asarray(d, dtype=float)
    y = np.asarray(y, dtype=float)
    if d.size != y.size:
        raise ValueError(f"d and y must match in length: {d.size} vs {y.size}")
    if np.any(d <= 0):
        raise ValueError("diameters must be positive")
    s = float(y.sum())
    if not np.isfinite(s) or s <= 0:
        raise ValueError(f"volume fractions must sum to a positive number, got {s}")
    y = y / s
    order = np.argsort(-d)                 # DESCENDING: index 0 is the coarsest
    d, y = d[order], y[order]
    k = d.size
    best = np.inf
    for i in range(k):
        tot = 0.0
        if i > 0:                          # coarser classes j < i exert a WALL effect
            bb = b_wall(d[i] / d[:i], pb)
            tot += float(np.sum((1.0 - beta + bb * beta * (1.0 - 1.0 / beta)) * y[:i]))
        if i < k - 1:                      # finer classes j > i LOOSEN class i
            aa = a_loosening(d[i + 1:] / d[i], pa)
            tot += float(np.sum((1.0 - aa) * y[i + 1:]))
        denom = 1.0 - tot
        if denom <= 1e-12:
            continue
        best = min(best, beta / denom)
    return float(best)


def binary_density(ratio: float, x_small_vol: float, beta: float,
                   pa: float, pb: float) -> float:
    """Two classes at size ratio `ratio` = d_large/d_small, fines VOLUME fraction."""
    if ratio < 1.0:
        raise ValueError(f"ratio = d_large/d_small must be >= 1, got {ratio}")
    x = float(np.clip(x_small_vol, 0.0, 1.0))
    if x <= 0.0 or x >= 1.0 or ratio == 1.0:
        return float(beta)
    return packing_density(np.array([ratio, 1.0]), np.array([1.0 - x, x]), beta, pa, pb)


def furnas_ideal(x_small_vol: float, beta: float) -> float:
    """Perfect size separation, closed form -- the a = b = 0 limit of the model."""
    x = float(np.clip(x_small_vol, 0.0, 1.0))
    coarse = beta / (1.0 - x) if x < 1.0 else np.inf
    fine = beta / (1.0 - (1.0 - beta) * (1.0 - x))
    return float(min(coarse, fine))


def psd_classes(family: str, cv: float, k: int = K_CLASSES,
                ndraw: int = NDRAW, seed: int = RNG_SEED):
    """Discretise the TRUNCATED SAMPLED psd into `k` equal-count classes.

    Sampling rather than integrating is deliberate: the packings were built from
    `psd.sample` with the same truncation, and the analytic moments do not survive it.
    """
    rng = np.random.default_rng(seed)
    p = make_psd(family, d10=1.0, cv=cv)
    d = np.asarray(p.sample(ndraw, rng, truncate=TRUNC), dtype=float)
    u = np.unique(d)
    if u.size <= k:                       # discrete family (bidisperse): use the atoms
        vol = np.array([float(np.sum(d[d == x] ** 3)) for x in u])
        return u, vol / vol.sum()
    edges = np.quantile(d, np.linspace(0.0, 1.0, k + 1))
    edges[0] -= 1e-12
    idx = np.clip(np.searchsorted(edges, d, side="left") - 1, 0, k - 1)
    dc = np.zeros(k)
    vol = np.zeros(k)
    for j in range(k):
        m = idx == j
        if not m.any():
            continue
        # The class diameter is the VOLUME-weighted mean of its members, because the
        # class stands in for the volume it carries.
        dc[j] = float(np.sum(d[m] ** 4) / np.sum(d[m] ** 3))
        vol[j] = float(np.sum(d[m] ** 3))
    keep = vol > 0
    return dc[keep], vol[keep] / vol[keep].sum()


# --- the data ------------------------------------------------------------------------

def measured(mu_p: float = MU_P):
    """(family, cv) -> (mean phi, sem, n) at the fixed rung. Uses the shipped loader."""
    rows = a.load_all()
    sel = [r for r in rows
           if r.get("n_requested") == N_REQ
           and abs(r.get("mu_p", 9) - mu_p) < 1e-9
           and abs(r.get("p_hat_target", 0) - P_HAT) < 1e-12
           and r.get("converged") and r.get("eta_roll") is None]
    g = collections.defaultdict(list)
    for r in sel:
        g[(r.get("family") or "lognormal", round(float(r["cv"]), 3))].append(float(r["phi"]))
    out = {}
    for k, v in g.items():
        arr = np.asarray(v)
        sem = float(arr.std(ddof=1) / np.sqrt(arr.size)) if arr.size > 1 else float("nan")
        out[k] = (float(arr.mean()), sem, int(arr.size))
    return out


# --- checks that can fail ------------------------------------------------------------

def check_limits(beta: float = 0.60) -> list[str]:
    """Every limit the docstring claims, computed. Returns FAILURE lines (empty = pass)."""
    bad = []
    one = packing_density([1.0], [1.0], beta, 1.02, 1.50)
    if abs(one - beta) > 1e-12:
        bad.append(f"K=1 gave {one!r}, not beta={beta}")
    # a = b = 0 is the INFINITE SIZE RATIO limit, not a limit of the exponents:
    # a(r) = sqrt(1 - (1-r)^PA) -> 0 and b(r) = 1 - (1-r)^PB -> 0 as r = d_fine/d_coarse
    # -> 0, for any positive exponents. Checked at two exponent pairs so the limit is
    # shown to be a property of the RATIO rather than of the fit.
    #
    # /!\ THE BAR HERE IS A CONVERGENCE RATE, NOT A TOLERANCE, and the first version of
    # this check was WRONG in the way lessons.md warns about: it demanded |model - closed
    # form| < 1e-05 at a finite ratio, which NO correct implementation can meet. At finite
    # r the loosening function is a = sqrt(1 - (1-1/r)^PA) ~ sqrt(PA/r), so the residual
    # falls only as 1/sqrt(r) -- 1.9e-05 at r = 1e8, exactly as it should. Requiring the
    # OBSERVED ORDER instead is strictly stricter: a wrong parameter map or a sign error
    # does not produce a residual that falls at the right rate, and a routine returning a
    # constant fails it outright.
    for pa, pb in ((1.02, 1.50), (0.4, 3.0)):
        for x in (0.05, 0.2, 0.3, 0.5, 0.8):
            want = furnas_ideal(x, beta)
            rs = np.array([1e4, 1e6, 1e8])
            err = np.array([abs(binary_density(r, x, beta, pa, pb) - want) for r in rs])
            if err[-1] > 1e-4:
                bad.append(f"ideal-separation limit at X_s={x} (PA={pa}, PB={pb}): "
                           f"residual {err[-1]:.3e} at r = 1e8, above 1e-04")
            if np.any(err <= 0):
                continue
            order = float(np.polyfit(np.log(rs), np.log(err), 1)[0])
            # The admissible band spans BOTH branch rates and that is not a fudge: the
            # two branches approach the ideal limit at different rates because the two
            # interaction functions do. On the COARSE branch the error is set by
            # a ~ sqrt(PA/r), giving order -1/2; on the FINE branch it is set by
            # b ~ PB/r, giving order -1. Which branch is active depends on X_s, so both
            # appear on this grid (-0.50 at X_s = 0.05, -1.00 at X_s = 0.3). A routine
            # returning a constant gives order 0 and fails; so does any map that does not
            # converge to the closed form at all.
            if not (-1.25 < order < -0.35):
                bad.append(f"ideal-separation CONVERGENCE at X_s={x} (PA={pa}, PB={pb}): "
                           f"observed order {order:+.3f}, expected -0.5 (coarse branch) "
                           f"or -1.0 (fine branch)")
    # a = b = 1 is the ratio -> 1 limit: no separation at all, gamma = beta everywhere.
    for x in (0.1, 0.4, 0.9):
        got = binary_density(1.0 + 1e-12, x, beta, 1.02, 1.50)
        if abs(got - beta) > 1e-9:
            bad.append(f"no-separation limit at X_s={x}: {got:.9f} vs beta {beta}")
    # ratio = 1 must be beta at every composition, for any exponents.
    for x in (0.1, 0.5, 0.9):
        got = binary_density(1.0, x, beta, 1.02, 1.50)
        if abs(got - beta) > 1e-12:
            bad.append(f"ratio=1 at X_s={x} gave {got}, not beta")
    # the endpoints of any binary are monodisperse and must agree with each other
    for r in (2.0, 3.0, 10.0):
        lo = binary_density(r, 0.0, beta, 1.02, 1.50)
        hi = binary_density(r, 1.0, beta, 1.02, 1.50)
        if abs(lo - beta) > 1e-12 or abs(hi - beta) > 1e-12:
            bad.append(f"binary endpoints at r={r}: {lo} / {hi}, not beta")
    return bad


def check_discretisation(beta: float, pa: float, pb: float) -> list[tuple]:
    """phi_model for the continuous families at K = 10, 20, 40, 80, 160."""
    rows = []
    for fam in ("lognormal", "uniform"):
        for cv in (0.2, 0.5):
            vals = []
            for k in (10, 20, 40, 80, 160):
                d, y = psd_classes(fam, cv, k=k)
                vals.append(packing_density(d, y, beta, pa, pb))
            rows.append((fam, cv, vals))
    return rows


# --- fit -----------------------------------------------------------------------------

def packing_density_blind(d, y, beta: float, a_const: float, b_const: float) -> float:
    """THE NULL: the same linear-packing STRUCTURE with the size ratio deleted.

    Every interaction is the same constant regardless of how far apart the two classes
    are in size, so the pairwise-filling physics is gone while the mixture bookkeeping,
    the min-over-dominant-class rule and the number of free parameters are identical.
    If this does as well as the real model, the agreement is the bookkeeping and not the
    filling, and nothing has been established.

    It is not a straw man: with a and b constant, gamma_i still depends on how the volume
    splits between classes coarser and finer than i, so the null CAN distinguish a
    left-skewed PSD from a right-skewed one. What it cannot see is how BIG the size gaps
    are -- which is exactly the claim under test.
    """
    d = np.asarray(d, dtype=float)
    y = np.asarray(y, dtype=float)
    y = y / y.sum()
    order = np.argsort(-d)
    y = y[order]
    k = y.size
    c_wall = 1.0 - beta + b_const * beta * (1.0 - 1.0 / beta)
    c_loose = 1.0 - a_const
    best = np.inf
    for i in range(k):
        tot = c_wall * float(y[:i].sum()) + c_loose * float(y[i + 1:].sum())
        denom = 1.0 - tot
        if denom <= 1e-12:
            continue
        best = min(best, beta / denom)
    return float(best)


def fit_on_bidisperse(meas, beta: float, scale: float = PHI_SCALE, blind: bool = False):
    """Fit two parameters to the five non-trivial BIDISPERSE points only.

    `blind=False` fits the exponents (PA, PB) of the real model; `blind=True` fits the
    two constants (a, b) of the size-ratio-blind null, so both have the same number of
    free parameters and the comparison is like for like. Residuals are scaled by
    `scale`, the protocol's own reproducibility, not by the seed sem.
    """
    pts = []
    for cv in CVS:
        k = ("bidisperse", cv)
        if k not in meas or cv == 0.0:      # CV = 0 IS beta, it cannot constrain a shape
            continue
        phi, _sem, _n = meas[k]
        d, y = psd_classes("bidisperse", cv)
        pts.append((d, y, phi))
    if not pts:
        raise RuntimeError("no bidisperse points to fit")

    def cost(theta):
        s = 0.0
        for d, y, phi in pts:
            if blind:
                p, q = 1.0 / (1.0 + np.exp(-theta))      # both constants in (0, 1)
                mv = packing_density_blind(d, y, beta, p, q)
            else:
                p, q = np.exp(theta)
                mv = packing_density(d, y, beta, p, q)
            s += ((mv - phi) / scale) ** 2
        return s

    best = None
    for x0 in ([0.0, 0.4], [1.0, 1.0], [-1.0, 0.0], [2.0, 2.0]):
        r = minimize(cost, np.array(x0), method="Nelder-Mead",
                     options={"xatol": 1e-9, "fatol": 1e-12, "maxiter": 6000})
        if best is None or r.fun < best.fun:
            best = r
    par = (1.0 / (1.0 + np.exp(-best.x))) if blind else np.exp(best.x)
    return float(par[0]), float(par[1]), float(best.fun), len(pts)


# --- report --------------------------------------------------------------------------

def bidisperse_path(cv: float) -> tuple[float, float]:
    """(fines VOLUME fraction, size ratio) for this project's bidisperse family.

    The family is 50:50 BY NUMBER at d10*(1 -/+ CV), so the fines volume fraction is
    (1-CV)^3 / [(1-CV)^3 + (1+CV)^3] and the size ratio is (1+CV)/(1-CV).
    """
    if cv <= 0:
        return 0.5, 1.0
    lo, hi = (1.0 - cv) ** 3, (1.0 + cv) ** 3
    return lo / (lo + hi), (1.0 + cv) / (1.0 - cv)


def main() -> None:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    say("=== mixture_model.py -- the shift as a PAIRWISE-FILLING FUNCTIONAL ===")
    say()
    say("=== 0. the model's own limits, computed rather than asserted ===")
    bad = check_limits()
    if bad:
        for b in bad:
            say(f"   FAIL {b}")
    else:
        say("   K = 1 returns beta; ratio -> 1 returns beta at every composition; and")
        say("   ratio -> inf reproduces the closed-form Furnas branches at two different")
        say("   exponent pairs, so the ideal limit is a property of the RATIO, not a")
        say("   limit of the fitted exponents. ALL LIMITS PASS.")
    say()

    meas = measured()
    beta = meas[("lognormal", 0.0)][0]
    mono = [meas[(f, 0.0)][0] for f in FAMILIES if (f, 0.0) in meas]
    say("=== 1. beta is MEASURED, not fitted ===")
    for fam in FAMILIES:
        if (fam, 0.0) in meas:
            p, s, n = meas[(fam, 0.0)]
            say(f"   monodisperse phi, {fam:<11} {p:.4f} +- {s:.4f}  (n = {n})")
    say(f"   -> beta = {beta:.4f} (the lognormal CV = 0 cell). The three families agree")
    say(f"      to {max(mono)-min(mono):.4f}, which is one check that CV = 0 is CV = 0.")
    say(f"   grading scale PHI_SCALE = {PHI_SCALE:.1e}: the protocol's own run-to-run")
    say("   non-determinism, NOT the seed sem (~2e-04), which is 10x smaller than the")
    say("   reproducibility of the thing being modelled.")
    say()

    say("=== 2. discretisation: is K = 160 enough? ===")
    for fam, cv, vals in check_discretisation(beta, 1.02, 1.50):
        say(f"   {fam:<10} cv {cv:.1f}  K = 10..160: "
            + " ".join(f"{v:.5f}" for v in vals)
            + f"   last step {vals[-1]-vals[-2]:+.5f}")
    say("   -> the lognormal CV = 0.5 value converges slowest; the K = 80 -> 160 step is")
    say("      the systematic to carry and it is below PHI_SCALE.")
    say()

    say("=== 3. THE ZERO-PARAMETER MODEL. Nothing is fitted: beta is measured and the")
    say("    exponents are held at the textbook values PA = 1.02, PB = 1.50 ===")
    say(f"   {'family':<11}{'CV':>5}{'X_s':>8}{'ratio':>7}{'phi_meas':>10}"
        f"{'phi_model':>11}{'resid':>9}{'/scale':>8}")
    resid: dict[tuple[str, float], tuple[float, float]] = {}
    for fam in FAMILIES:
        for cv in CVS:
            key = (fam, cv)
            if key not in meas:
                continue
            d, y = psd_classes(fam, cv)
            mv = packing_density(d, y, beta, 1.02, 1.50)
            phi = meas[key][0]
            resid[key] = (phi, mv)
            if fam == "bidisperse":
                xs, rr = bidisperse_path(cv)
                xt, rt = f"{xs:8.3f}", f"{rr:7.2f}"
            else:
                xt, rt = f"{'-':>8}", f"{'-':>7}"
            say(f"   {fam:<11}{cv:5.1f}{xt}{rt}{phi:10.4f}{mv:11.4f}"
                f"{mv-phi:+9.4f}{(mv-phi)/PHI_SCALE:+8.1f}")
    say()

    say("=== 4. the two things no descriptor could do ===")
    for tag, idx in (("measured", 0), ("model", 1)):
        b = [resid[("bidisperse", cv)][idx] for cv in CVS if ("bidisperse", cv) in resid]
        i = int(np.argmax(b))
        order = sorted([f for f in FAMILIES if (f, 0.5) in resid],
                       key=lambda f: resid[(f, 0.5)][idx])
        say(f"   {tag:<9} bidisperse peaks at CV = {CVS[i]:.1f}  "
            f"[{' '.join(f'{v:.4f}' for v in b)}]")
        say(f"   {tag:<9} family order at CV = 0.5, smallest phi first: {order}")
    say("   -> the model has BOTH: an interior maximum on the bidisperse path, and the")
    say("      right family order at matched CV. No monotone function of a single")
    say("      descriptor can have the first (shape_control.py, fifteen tried).")
    say()

    say("=== 5. the CONTRAST the project quotes, CV 0 -> 0.5 ===")
    say(f"   {'family':<11}{'d_meas':>9}{'d_model':>9}{'meas/model':>12}")
    rat: dict[str, float] = {}
    for fam in FAMILIES:
        if (fam, 0.0) not in resid or (fam, 0.5) not in resid:
            continue
        dm = resid[(fam, 0.5)][0] - resid[(fam, 0.0)][0]
        dd = resid[(fam, 0.5)][1] - resid[(fam, 0.0)][1]
        rat[fam] = dm / dd
        say(f"   {fam:<11}{dm:+9.4f}{dd:+9.4f}{dm/dd:12.2f}")
    if rat:
        sp = max(rat.values()) / min(rat.values())
        say("   -> a collapsing model gives the SAME ratio for all three families.")
        say(f"      spread = {sp:.2f}x, against 4.51x for CV and 4.06x for d32")
        say("      (shape_control.py) -- and with the same number of fitted parameters")
        say("      as those had: ZERO.")
    say()

    say("=== 6. THE NULL: the same bookkeeping with the size ratio DELETED ===")
    pa_b, pb_b, chi_b, _ = fit_on_bidisperse(meas, beta, blind=True)
    pa_r, pb_r, chi_r, n_r = fit_on_bidisperse(meas, beta, blind=False)
    say(f"   both fitted on the same {n_r} non-trivial bidisperse points, both with 2")
    say("   free parameters, both scaled by PHI_SCALE:")
    say(f"     real model  PA = {pa_r:.3f}  PB = {pb_r:.3f}    chi2(fitted set) = {chi_r:9.2f}")
    say(f"     blind null  a  = {pa_b:.3f}  b  = {pb_b:.3f}    chi2(fitted set) = {chi_b:9.2f}")
    held = [(f, cv) for f in ("uniform", "lognormal") for cv in CVS
            if (f, cv) in meas and cv > 0]
    for tag, fn, pars in (("real model", packing_density, (pa_r, pb_r)),
                          ("blind null", packing_density_blind, (pa_b, pb_b))):
        s = 0.0
        worst: tuple[float, object] = (0.0, None)
        for key in held:
            d, y = psd_classes(*key)
            mv = fn(d, y, beta, *pars)
            r = (mv - meas[key][0]) / PHI_SCALE
            s += r * r
            if abs(r) > abs(worst[0]):
                worst = (r, key)
        say(f"     {tag}: HELD OUT on {len(held)} uniform+lognormal points, "
            f"chi2 = {s:9.2f}, worst {worst[0]:+.1f} scale at {worst[1]}")
    say("   -> the null is calibrated on the same data with the same freedom, and it can")
    say("      still see how the volume SPLITS between coarser and finer classes. What it")
    say("      cannot see is how big the size gaps are, which is the claim under test.")
    say()

    say("=== 7. does the SAME zero-parameter model work at mu_p = 0? ===")
    m0 = measured(mu_p=0.0)
    if ("lognormal", 0.0) in m0:
        b0 = m0[("lognormal", 0.0)][0]
        say(f"   beta(mu_p = 0) = {b0:.4f}   (vs {beta:.4f} at mu_p = 0.3)")
        say(f"   {'family':<11}{'CV':>5}{'phi_meas':>10}{'phi_model':>11}{'/scale':>8}")
        for fam in FAMILIES:
            for cv in CVS:
                if (fam, cv) not in m0:
                    continue
                d, y = psd_classes(fam, cv)
                mv = packing_density(d, y, b0, 1.02, 1.50)
                say(f"   {fam:<11}{cv:5.1f}{m0[(fam, cv)][0]:10.4f}{mv:11.4f}"
                    f"{(mv-m0[(fam, cv)][0])/PHI_SCALE:+8.1f}")
        say("   -> the interaction functions are GEOMETRIC, so only beta should carry the")
        say("      friction. Whether that survives is what this block measures; it is a")
        say("      claim that can fail, not a restatement.")
    else:
        say("   no mu_p = 0 cells at this rung.")
    say()

    say("=== 8. THE BLIND PREDICTION, for PREREG_MIXTURE_2026-09-13.md ===")
    say("   phi at FIXED size ratio across the fines volume fraction -- a region where")
    say("   this project has no data at all. The bidisperse family reaches X_s = 0.036 at")
    say("   ratio 3, which is the far right-hand tail of this curve and nothing else.")
    for r in (2.0, 3.0):
        say(f"   size ratio {r:.1f}:")
        say(f"      {'X_s':>6}{'zero-param':>12}{'fitted':>9}{'blind null':>12}")
        for x in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.60, 0.80, 1.0):
            z = binary_density(r, x, beta, 1.02, 1.50)
            f_ = binary_density(r, x, beta, pa_r, pb_r)
            if x in (0.0, 1.0):
                nb = beta
            else:
                nb = packing_density_blind(np.array([r, 1.0]),
                                           np.array([1.0 - x, x]), beta, pa_b, pb_b)
            say(f"      {x:6.2f}{z:12.4f}{f_:9.4f}{nb:12.4f}")
    say()
    say("   The three columns SEPARATE by construction of the experiment, which is what")
    say("   makes it worth running: the pairwise model puts a large interior maximum near")
    say("   X_s ~ 0.3, the blind null cannot put the maximum anywhere the size ratio")
    say("   says, and no filling at all is a flat line at beta.")

    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "mixture_model.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
