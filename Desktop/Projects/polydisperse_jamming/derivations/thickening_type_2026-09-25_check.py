"""thickening_type_2026-09-25_check.py -- where a Wyart-Cates suspension is CST, DST or shear-jamming.

Closure = the one wc_drop_impact.c integrates (steady state, fhat = fss):
    f(sigma) = exp(-(sigma*/sigma)^delta),  phi_J = phi_0 [1 - f (1 - chi)],  eta = eta_s (1 - phi/phi_J)^-n
In x = phi/phi_0 and u = (sigma*/sigma)^delta, the steady flow curve is gammadot(sigma) = (sigma/eta_s) g^n with
    g = 1 - x/J,  J = 1 - (1 - chi) e^{-u}.
DST (an S-shaped flow curve, d gammadot/d sigma < 0 somewhere) <=> max_u H(u) > 1 with
    H(u) = n delta x (1 - chi) u e^{-u} / [J (J - x)]            (derived in thickening_type_2026-09-25.md)
Sections: C1 SymPy derivation of H from d ln gammadot / d ln sigma; C2 brute-force check of the DST criterion
against a direct numerical flow curve at 4 parameter points; C3 x_DST(chi, delta, n) on the project's measured
(chi, delta) for spheres, doublets and CV; C4 limits. Writes thickening_type_2026-09-25_check.txt.
"""
from __future__ import annotations


import os

import numpy as np
import sympy as sp
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
L = []


def say(s=""):
    print(s)
    L.append(s)


# ---------------------------------------------------------------- C1
say("C1  SymPy: d ln(gammadot)/d ln(sigma) = 1 - H(u)")
s, x, chi, n, d = sp.symbols("s x chi n delta", positive=True)
u = s ** (-d)
J = 1 - (1 - chi) * sp.exp(-u)
lng = sp.log(s) + n * sp.log(1 - x / J)
slope = sp.simplify(sp.diff(lng, s) * s)
U = sp.symbols("u", positive=True)
H = n * d * x * (1 - chi) * U * sp.exp(-U) / ((1 - (1 - chi) * sp.exp(-U)) * (1 - (1 - chi) * sp.exp(-U) - x))
check = sp.simplify((slope - (1 - H.subs(U, u))))
say(f"   slope - (1 - H) simplifies to: {check}")
assert check == 0


# ---------------------------------------------------------------- helpers
def Hmax(xv, chiv, dv, nv):
    uu = np.logspace(-4, 3, 20001)
    Jv = 1 - (1 - chiv) * np.exp(-uu)
    ok = Jv > xv
    h = np.where(ok, nv * dv * xv * (1 - chiv) * uu * np.exp(-uu) / (Jv * np.maximum(Jv - xv, 1e-300)), np.inf)
    return float(h.max())


def x_dst(chiv, dv, nv):
    """Smallest x < chi with max_u H = 1 (DST onset). H rises monotonically with x at fixed u."""
    f = lambda xv: Hmax(xv, chiv, dv, nv) - 1.0  # noqa: E731
    lo, hi = 1e-6, chiv * (1 - 1e-9)
    if f(lo) > 0:
        return 0.0
    if f(hi) < 0:
        return float("nan")
    return brentq(f, lo, hi, xtol=1e-10)


def flowcurve_nonmonotone(xv, chiv, dv, nv):
    ss = np.logspace(-3, 3, 200001)
    uu = ss ** (-dv)
    Jv = 1 - (1 - chiv) * np.exp(-uu)
    g = 1 - xv / Jv
    gd = ss * g ** nv
    return bool(np.any(np.diff(gd) < 0))


# ---------------------------------------------------------------- C2
say("\nC2  brute force: criterion vs direct flow curve (non-monotone gammadot(sigma))")
for chiv, dv, nv in ((0.908, 1.254, 2.0), (0.863, 1.063, 2.0), (0.908, 1.254, 5.0), (0.80, 2.0, 2.0)):
    xd = x_dst(chiv, dv, nv)
    below, above = flowcurve_nonmonotone(xd - 2e-3, chiv, dv, nv), flowcurve_nonmonotone(xd + 2e-3, chiv, dv, nv)
    say(f"   chi {chiv} delta {dv} n {nv}: x_DST = {xd:.5f}; flow curve non-monotone at x_DST - 0.002: {below}, + 0.002: {above}")
    assert (not below) and above

# ---------------------------------------------------------------- C3
say("\nC3  the CST / DST / SJ bands in x = phi/phi_0 for measured (chi, delta)  [chi, delta sources in the .md]")
cases = [
    ("spheres, mu_p 0.3, CV 0 (sheared delta)", 0.908, 1.254),
    ("pure doublets, mu_p 0.3 (sheared delta)", 0.863, 1.063),
    ("spheres, CV 0.5 (delta x0.81 of CV 0)", 0.917, 1.254 * 0.81),
    ("clumps alpha 1.1-2.0, p 1e-3 chi (delta not measured: sphere value)", 0.894, 1.254),
]
for nv in (2.0, 5.0):
    say(f"   n = {nv}:")
    for lab, chiv, dv in cases:
        xd = x_dst(chiv, dv, nv)
        say(f"      {lab:62s} x_DST {xd:.4f}  DST band {chiv - xd:.4f}  SJ band {1 - chiv:.4f}")

# ---------------------------------------------------------------- C4
say("\nC4  limits")
say(f"   delta -> 0 (no activation contrast): no DST at all (nan = never crosses 1): delta 0.05: {x_dst(0.9, 0.05, 2.0):.4f} (chi 0.9)")
say(f"   delta large: x_DST at delta 20: {x_dst(0.9, 20.0, 2.0):.4f}")
say(f"   chi -> 1 (no frictional drop): H -> 0, no DST: Hmax(x 0.99, chi 0.999999) = {Hmax(0.99, 0.999999, 1.25, 2.0):.2e}")
say("   x -> chi: for small u, J - x ~ (chi - x) + (1 - chi) u, so sup_u H -> n*delta*x/J -> n*delta. Numerically:")
for dv, nv in ((0.4, 2.0), (0.6, 2.0), (1.254, 2.0), (0.3, 5.0)):
    uu = np.logspace(-8, 3, 200001)
    Jv = 1 - 0.1 * np.exp(-uu)
    hm = float((nv * dv * (0.9 - 1e-7) * 0.1 * uu * np.exp(-uu) / (Jv * (Jv - (0.9 - 1e-7)))).max())
    say(f"      n*delta {nv*dv:.3f}: sup H at x = chi - 1e-7 (chi 0.9) = {hm:.4f}")
say("   => a DST band below chi exists IFF n*delta > 1; for n*delta < 1 the suspension goes CST -> SJ directly")
say(f"   scaling: H is LINEAR in n*delta, so x_DST depends on n and delta only through n*delta:"
    f" x_DST(chi .9, delta 2, n 2) = {x_dst(0.9, 2.0, 2.0):.6f} vs (delta 1, n 4) = {x_dst(0.9, 1.0, 4.0):.6f}")

with open(os.path.join(HERE, "thickening_type_2026-09-25_check.txt"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
