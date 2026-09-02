"""Closed forms for the quasi-steady shortfall, and where they stop being true.

`SCALING_2026-08-31.md` §4 reports the shortfall as a table of integrated numbers.
A table is not understanding, and it gives nothing to check the integrator
against. This derives the same curve analytically under one stated approximation,
verifies it symbolically and numerically, and then measures the error of the
approximation rather than assuming it away.

Run:  python derivations/shortfall_asymptotics.py

Result summary lives in `derivations/SHORTFALL_ASYMPTOTICS.md`, written by hand
from this script's output.

The approximation
-----------------
The exact problem is nonlinear because the stress feeds back on its own
activation: f grows, viscosity grows, stress grows, f_ss grows. Freeze that
feedback by evaluating f_ss once, at the *steady* stress. Then the strain clock
integrates exactly:

    f(A) = f_ss * (1 - exp(-A)),          A = gamma_dot tau / gamma_0

and the shortfall follows in closed form. Because the true stress starts at the
lubricated value and climbs, the true f_ss is smaller than the frozen one
throughout, so this is an UPPER BOUND on the delivered stress and therefore on
the ratio. §4 of the output checks that the bound is never violated.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "layer2"))

from activation import STFMaterial, quasi_steady_ratio, steady_stress  # noqa: E402


# ---------------------------------------------------------------- symbolic ---

def symbolic() -> dict:
    """Derive eta(f), its logarithmic slope, and the crossover activation number."""
    f, fs, A, phi, phi_0, phi_m, eta_s, n, eps = sp.symbols(
        "f f_s A phi phi_0 phi_m eta_s n epsilon", positive=True
    )
    D = phi_0 - phi_m                      # width of the jamming interval
    phi_J = phi_0 - f * D                  # (1-f) phi_0 + f phi_m, expanded
    eta = eta_s * ((phi_J - phi) / phi_J) ** (-n)

    # Check the expanded form really is the interpolation dst_fit uses.
    assert sp.simplify(phi_J - ((1 - f) * phi_0 + f * phi_m)) == 0

    # d ln eta / df, the amplification of a change in frictional fraction.
    dlneta = sp.simplify(sp.diff(sp.log(eta), f))
    dlneta_closed = n * D * phi / ((phi_0 - f * D) * (phi_0 - f * D - phi))
    assert sp.simplify(dlneta - dlneta_closed) == 0, sp.simplify(dlneta - dlneta_closed)

    # Frozen-feedback trajectory and the resulting stress ratio.
    f_of_A = fs * (1 - sp.exp(-A))
    ratio = (eta.subs(f, f_of_A) / eta.subs(f, fs))

    # Small-A expansion about the plateau eta(0)/eta(f_s).
    plateau = sp.simplify(eta.subs(f, 0) / eta.subs(f, fs))
    slope0 = sp.simplify(sp.diff(ratio, A).subs(A, 0))
    slope0_closed = sp.simplify(plateau * n * D * phi * fs / (phi_0 * (phi_0 - phi)))
    assert sp.simplify(slope0 - slope0_closed) == 0, sp.simplify(slope0 - slope0_closed)

    # Crossover: the f at which the ratio reaches 1 - eps, then invert the clock.
    g = (phi_0 - f * D) / (phi_0 - f * D - phi)
    # eta = eta_s g^n, so an eta ratio of (1-eps) is a g ratio of (1-eps)^(1/n).
    G = (1 - eps) ** (sp.Rational(1, 1) / n) * g.subs(f, fs)
    f_star = sp.simplify((phi_0 - G * phi / (G - 1)) / D)
    A_star = -sp.log(1 - f_star / fs)

    return {
        "eta": eta, "ratio": ratio, "plateau": plateau,
        "dlneta": dlneta_closed, "slope0": slope0_closed,
        "f_star": f_star, "A_star": A_star,
        "syms": (f, fs, A, phi, phi_0, phi_m, eta_s, n, eps),
    }


# ---------------------------------------------------------------- numeric ----

def ratio_frozen(A: np.ndarray, mat: STFMaterial, gamma_dot: float) -> np.ndarray:
    """Closed-form shortfall under the frozen-feedback approximation."""
    sigma_ss = steady_stress(gamma_dot, mat)
    f_s = float(mat.f_steady(sigma_ss))
    f = f_s * (1.0 - np.exp(-np.asarray(A, dtype=float)))
    return np.array([float(mat.viscosity(fi)) for fi in f]) / float(mat.viscosity(f_s))


def a_star_closed(mat: STFMaterial, gamma_dot: float, eps: float) -> float:
    """Activation number at which the frozen-feedback ratio reaches 1 - eps."""
    sigma_ss = steady_stress(gamma_dot, mat)
    f_s = float(mat.f_steady(sigma_ss))
    D = mat.phi_0 - mat.phi_m
    g_s = (mat.phi_0 - f_s * D) / (mat.phi_0 - f_s * D - mat.phi)
    G = (1.0 - eps) ** (1.0 / mat.n) * g_s
    f_star = (mat.phi_0 - G * mat.phi / (G - 1.0)) / D
    if not 0.0 < f_star < f_s:
        return float("nan")
    return float(-np.log(1.0 - f_star / f_s))


def main() -> int:
    S = symbolic()
    print("=" * 74)
    print("1. SYMBOLIC IDENTITIES  (each asserted inside symbolic(), so a failure raises)")
    print("=" * 74)
    print("   phi_J(f)          = phi_0 - f (phi_0 - phi_m)     [= (1-f)phi_0 + f phi_m]")
    print("   d ln eta / d f    =", sp.simplify(S["dlneta"]))
    print("   plateau (A -> 0)  = eta(0) / eta(f_s)")
    print("   d ratio / d A |_0 = plateau * n D phi f_s / (phi_0 (phi_0 - phi))")
    print("   f*(eps)           =", sp.simplify(S["f_star"]))

    SS, PHI_M, PHI_0, N, GD = 50.0, 0.54, 0.64, 2.0, 100.0
    ETA_S = SS / ((1 - 0.40 / PHI_0) ** -N * GD)
    phis = (0.40, 0.46, 0.50, 0.52)
    mats = {p: STFMaterial(phi=p, phi_m=PHI_M, sigma_star=SS, eta_s=ETA_S,
                           gamma_0=1.0, phi_0=PHI_0, n=N) for p in phis}

    print()
    print("=" * 74)
    print("2. PLATEAU: closed form vs the integrator, at four packing fractions")
    print("=" * 74)
    print("   phi    integrated      eta(0)/eta(f_s)      rel.diff")
    for p in phis:
        m = mats[p]
        integ, _, _ = quasi_steady_ratio(GD, 1e-6 * m.gamma_0 / GD, m)
        closed = ratio_frozen(np.array([0.0]), m, GD)[0]
        print(f"   {p:.2f}  {integ:12.6f}  {closed:19.6f}   {abs(integ-closed)/closed:9.2e}")

    print()
    print("=" * 74)
    print("3. SMALL-A SLOPE: symbolic expression vs a finite difference of the closed form")
    print("=" * 74)
    f_, fs_, A_, phi_, phi0_, phim_, etas_, n_, eps_ = S["syms"]
    print("   phi     symbolic        finite-diff        rel.diff")
    for p in phis:
        m = mats[p]
        f_s = float(m.f_steady(steady_stress(GD, m)))
        sym = float(S["slope0"].subs({phi_: p, phi0_: PHI_0, phim_: PHI_M,
                                      n_: N, fs_: f_s}))
        h = 1e-6
        fd = (ratio_frozen(np.array([h]), m, GD)[0] - ratio_frozen(np.array([0.0]), m, GD)[0]) / h
        print(f"   {p:.2f}  {sym:12.6f}  {fd:16.6f}   {abs(sym-fd)/abs(sym):9.2e}")

    print()
    print("=" * 74)
    print("4. THE BOUND: frozen feedback must never sit below the full integration")
    print("=" * 74)
    print("   phi     A      integrated     frozen      frozen-integ   bound ok")
    worst = 0.0
    for p in phis:
        m = mats[p]
        for A in (0.1, 0.3, 1.0, 3.0, 10.0):
            integ, _, _ = quasi_steady_ratio(GD, A * m.gamma_0 / GD, m)
            fr = ratio_frozen(np.array([A]), m, GD)[0]
            gap = fr - integ
            worst = max(worst, gap)
            print(f"   {p:.2f} {A:6.1f} {integ:12.6f} {fr:12.6f} {gap:14.6f}   {gap >= -1e-9}")
    print(f"\n   worst overshoot of the bound: {worst:.4f}")
    print("   (a large gap means the stress feedback is strong, i.e. the frozen")
    print("    approximation is a loose bound there, not that either is wrong)")

    print()
    print("=" * 74)
    print("5. CROSSOVER A*(eps = 0.01): closed form vs bisection on the integrator")
    print("=" * 74)
    print("   phi    A* closed    A* integrated    note")
    for p in phis:
        m = mats[p]
        ac = a_star_closed(m, GD, 0.01)
        lo, hi = 1e-3, 1e3
        for _ in range(60):
            mid = np.sqrt(lo * hi)
            r, _, _ = quasi_steady_ratio(GD, mid * m.gamma_0 / GD, m)
            lo, hi = (mid, hi) if r < 0.99 else (lo, mid)
        ai = np.sqrt(lo * hi)
        print(f"   {p:.2f}  {ac:10.3f}  {ai:14.3f}    "
              f"{'closed form is the bound, so A* closed < A* true' if ac < ai else 'CHECK'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
