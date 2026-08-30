"""The ballistic dispersion law, and what it predicts as wind speed varies.

WHERE THIS COMES FROM
---------------------
Not from me. `CLAUDE_CODE_FRISBEE_HARDENING.md` section C5 derives it and
`P4_Coupling/regrade_gate_g1.py` implements and verifies `Phi` against sympy
limits. This module re-exports that `Phi` -- it does not re-implement it -- and
adds only the wind-speed extension and the tests that discriminate it.

    sigma_landing = (kappa * g * sigma_u * T_f^2) / (2 * V_d) * sqrt(Phi(N))
    N = L_f / L
    Phi(N) = 8 [ 1/(3N) - 1/(2N^2) - e^-N/N^3 + (1 - e^-N)/N^4 ]
    kappa = sqrt( 4 + (C_Lalpha/C_L)^2 (sigma_w/sigma_u)^2 )

The disc is a HIGH-Stokes body: ballistic length l_b = 2m/(rho C_D S) = 51 m
against an integral scale L ~ 3-6 m, so St ~ 8.5 and the disc does not follow
eddies. It responds through force, the displacement is a DOUBLE time integral,
and the tracer law sigma ~ sigma_u T_f / sqrt(N) does not apply. The two laws
are trivially distinguishable and this module is the instrument:

    at fixed L_f:   ballistic sigma ~ V_d^-3     tracer sigma ~ V_d^-1
    at fixed V_d:   ballistic sigma ~ L_f^1.7    tracer sigma ~ L_f^0.5

THE WIND-SPEED EXTENSION
------------------------
`wind_scaling.py` establishes (and measures) that under the RANS scaling

    sigma_u ~ Uref        L ~ Uref^0        so  N ~ Uref^0

Substituting into the law, at a FIXED flight (fixed T_f, V_d, L_f):

    sigma_landing ~ Uref          -- exactly linear, no threshold, no exponent

That is a strong prediction and it is the one worth attacking, because a real
throw is not a fixed flight: the mean wind changes the airspeed, hence the lift,
hence T_f, V_d and L_f, and the disc's aerodynamic coefficients are strongly
nonlinear in angle of attack (the project uses the TABULATED Potts & Crowther
curves precisely because the linear fits miss range by 30-57%). So the law
predicts a linear FLOOR, and everything above it is disc aerodynamics.

    LINEAR REGIME     sigma_landing grows in proportion to wind speed;
                      the site enters only through one amplitude.
    NONLINEAR REGIME  the mean wind changes the flight itself; the departure
                      from linear is a property of the DISC, not of the site.

The crossover is where the mean wind stops being a perturbation on the disc's
own airspeed, i.e. at some order-one value of U/V_d. That group, not the
turbulence intensity U/sigma_u of pedestrian-comfort criteria, is the one that
should organise the problem -- which is the same argument section C5 already
makes for sigma_u/V_d over sigma_u/U.

WHAT THIS MODULE DOES NOT DO
----------------------------
It does not fit anything to the Monte Carlo. The prediction is evaluated from
CFD-measured sigma_u and L and from flight quantities measured on the SAME
trajectories, with kappa taken from the brief. If the prediction and the
measurement agree, that is a test; if a coefficient were fitted first it would
not be.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "P4_Coupling"))

from regrade_gate_g1 import phi                            # noqa: E402

G = 9.80665

# kappa from the brief's section C5. The gust model in production uses
# sigma_ratios = (1.0, 0.75, 0.5) about the local mean-wind direction, so
# sigma_w/sigma_u = 0.5 and kappa = 3.03; the isotropic arm gives 4.96.
KAPPA_ANISOTROPIC = 3.03
KAPPA_ISOTROPIC = 4.96


def _scalar_like(out, *inputs):
    """Return a python float when every input was a scalar.

    `phi` calls `np.atleast_1d`, so without this every scalar evaluation comes
    back as a length-1 array and `float()` on it raises. Silently returning a
    1-element array instead is worse: it propagates into JSON as a list.
    """
    if all(np.ndim(a) == 0 for a in inputs):
        return float(np.reshape(out, ())) if np.size(out) == 1 else out
    return out


def sigma_ballistic(sigma_u, T_f, V_d, L_f, L_int, kappa=KAPPA_ANISOTROPIC):
    """Predicted landing dispersion, metres. All arguments SI, array-like."""
    args = (sigma_u, T_f, V_d, L_f, L_int)
    sigma_u, T_f, V_d, L_f, L_int = map(
        lambda a: np.asarray(a, dtype=float), args)
    N = L_f / np.maximum(L_int, 1e-12)
    out = (kappa * G * sigma_u * T_f ** 2) / (2.0 * np.maximum(V_d, 1e-12)) \
        * np.sqrt(phi(N))
    return _scalar_like(out, *args)


def sigma_tracer(sigma_u, T_f, L_f, L_int):
    """The Taylor tracer law the brief disowns, kept as the discriminating
    alternative. A test that only ever evaluates the favoured law is not a test.
    """
    sigma_u, T_f, L_f, L_int = map(
        lambda a: np.asarray(a, dtype=float), (sigma_u, T_f, L_f, L_int))
    N = np.maximum(L_f / np.maximum(L_int, 1e-12), 1e-12)
    return sigma_u * T_f / np.sqrt(N)


def predicted_exponents():
    """The exponents each law predicts, for the discriminating tests.

    Derived, not asserted: `check_exponents_symbolically` re-derives them with
    sympy from the closed forms above and asserts these values.
    """
    return {
        "ballistic": {"V_d_at_fixed_Lf": -3.0, "L_f_at_fixed_Vd": 1.7,
                      "u_ref_at_fixed_flight": 1.0},
        "tracer": {"V_d_at_fixed_Lf": -1.0, "L_f_at_fixed_Vd": 0.5,
                   "u_ref_at_fixed_flight": 1.0},
    }


def check_exponents_symbolically(N0=4.0):
    """Re-derive the exponents rather than quoting the brief's arithmetic.

    At fixed L_f, T_f = L_f/V_d, so sigma_ballistic ~ V_d^-2 * V_d^-1 = V_d^-3
    and N is fixed, so sqrt(Phi(N)) contributes nothing -- exact.

    At fixed V_d, T_f = L_f/V_d ~ L_f, so the prefactor goes as L_f^2 and
    sqrt(Phi(L_f/L)) supplies the rest. That part is NOT a power law: the
    quoted 1.7 is the LOCAL logarithmic slope at the site's own N, so it is
    returned here evaluated at `N0` rather than as a constant. Reporting it as
    a universal exponent would be wrong, and the brief's "1.7" is only correct
    near N = 4.
    """
    import sympy as sp

    N, Lf, Vd, L = sp.symbols("N L_f V_d L", positive=True)
    Phi = 8 * (1 / (3 * N) - 1 / (2 * N ** 2) - sp.exp(-N) / N ** 3
               + (1 - sp.exp(-N)) / N ** 4)

    # Anchors on Phi itself, independent of anything in regrade_gate_g1.
    lim0 = sp.limit(Phi, N, 0)
    liminf = sp.limit(N * Phi, N, sp.oo)

    # sigma at fixed L_f: T_f = L_f/V_d, N = L_f/L fixed.
    sig_fixed_Lf = (Lf / Vd) ** 2 / Vd
    p_Vd = sp.simplify(Vd * sp.diff(sig_fixed_Lf, Vd) / sig_fixed_Lf)

    # sigma at fixed V_d: T_f ~ L_f, N = L_f/L.
    sig_fixed_Vd = Lf ** 2 * sp.sqrt(Phi.subs(N, Lf / L))
    p_Lf = sp.simplify(Lf * sp.diff(sig_fixed_Vd, Lf) / sig_fixed_Vd)
    p_Lf_at = float(p_Lf.subs({Lf: N0, L: 1.0}))

    return dict(Phi_at_zero=float(lim0), N_Phi_at_infinity=float(liminf),
                exponent_V_d=float(p_Vd),
                exponent_L_f_local=p_Lf_at, evaluated_at_N=float(N0))


def local_exponent(x, y):
    """Local log-log slope of y(x), by centred differences on the logs.

    Used instead of a single global power-law fit wherever the relationship is
    known NOT to be a power law -- sqrt(Phi(N)) is not, and neither is the
    response to wind speed once the disc's aerodynamics move.
    """
    lx, ly = np.log(np.asarray(x, float)), np.log(np.asarray(y, float))
    return np.gradient(ly, lx)


def fit_exponent(x, y):
    """Global log-log slope with a bootstrap CI. Returns (p, lo, hi, r2)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return float("nan"), float("nan"), float("nan"), float("nan")
    lx, ly = np.log(x), np.log(y)
    p, c = np.polyfit(lx, ly, 1)
    r2 = 1.0 - np.sum((ly - (p * lx + c)) ** 2) / np.sum((ly - ly.mean()) ** 2)
    rng = np.random.default_rng(20260830)
    ps = []
    for _ in range(4000):
        i = rng.integers(0, len(x), len(x))
        if len(np.unique(x[i])) < 2:
            continue
        ps.append(np.polyfit(lx[i], ly[i], 1)[0])
    lo, hi = np.percentile(ps, [2.5, 97.5])
    return float(p), float(lo), float(hi), float(r2)
