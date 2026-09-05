"""Needle injectability: pressure drop and plunger force for a generalised-Newtonian
fluid in laminar, fully developed flow through a circular needle.

Closed form
-----------
For Herschel-Bulkley (``sigma = tau_y + k gd^n``) the flow rate at pressure drop ``dp``
through a needle of radius ``R`` and length ``L`` is (``derivations/buckingham_reiner.md``)

    tau_w = dp R / (2 L),   S = tau_w - tau_y,
    Q = pi R^3 / tau_w^3 * k^(-1/n) * [ S^(1/n+3)/(1/n+3) + 2 tau_y S^(1/n+2)/(1/n+2)
                                        + tau_y^2 S^(1/n+1)/(1/n+1) ],

``Q = 0`` when ``tau_w <= tau_y`` (the plug ``r_p = 2 tau_y L / dp`` fills the needle).
Its ``n = 1`` case is the Buckingham-Reiner equation and its ``tau_y = 0`` case the
power-law result; both, and Poiseuille, are the anchors in ``tests/test_injectability.py``.

For any other registry model (Cross, Carreau-Yasuda, Casson, ...) :func:`flow_rate_gnf`
integrates ``Q = pi int_0^R r^2 gd(r) dr`` numerically, inverting ``sigma(gd) = tau(r)``
at each quadrature node. The two routes are checked against each other on HB.

What is and is not in the force
-------------------------------
``F = dp * pi R_s^2 + F_friction``: the needle's viscous pressure drop acting on the
plunger of barrel radius ``R_s``, plus a friction term the caller supplies (it is a
property of the syringe, not the fluid). NOT included: the barrel's own pressure drop
(smaller by ``(R/R_s)^4`` for a Newtonian fluid), entrance and exit losses, the kinetic
term ``rho U^2``, and the transient at plunger start. These are all inertial or
geometric corrections that grow with Reynolds number; pass ``rho`` to get the wall
Reynolds number reported beside the force and judge them. ALSO NOT included, and not
inertial: wall slip. The closed form assumes no slip at the needle wall; microgel
suspensions (Carbopol) and many biological fluids slip at wall stresses of the kPa
order this geometry produces, which lowers the true force below the number here and
is not covered by the extrapolation flag either (it is a boundary condition, not a
constitutive extrapolation). A slip-vs-no-slip test needs two needle radii at the
same wall stress (Mooney's method); none of the registry data has it.

``reynolds`` is ``rho U (2R) / eta(gd_w)``, the apparent viscosity at the true wall
shear rate, NOT the Metzner-Reed generalised Reynolds number built on ``8U/D``; for a
power-law fluid ``Re_here = Re_MR (3n+1)/(4n)``, 1.1-1.3x at n = 0.4-0.65. Compare
with a literature transition threshold only after converting.

Where the flow curve is being extrapolated
------------------------------------------
The shear rate in the needle runs from 0 at the plug edge (or axis) to ``gd_w`` at the
wall, and the wall region carries most of ``Q``. :func:`force_with_uq` compares ``gd_w``
with the identified model's fitted envelope and flags the force as extrapolated when the
wall rate is outside it. That flag is the answer to "did your rheometer measure the rate
the needle imposes", which is the first question a reviewer asks.

Uncertainty
-----------
:func:`force_with_uq` reports two intervals and says which is which. (1) Profile
push-through: each parameter is moved to the ends of its profile-likelihood interval
with the others held at the fit, and the force envelope over those moves is reported.
This honours the (often asymmetric, often open) profile intervals but ignores parameter
correlation. (2) Covariance sampling: ``n_samples`` draws from the Wald normal on the
parameters, forced by rejection to stay inside the model bounds, with force quantiles
at the identified level. This carries the correlations but assumes a symmetric normal
on parameters that live on decades. When the two disagree the profile interval is
the one to quote; when the covariance is singular (degenerate fit) only the profile
interval exists and the sampled interval is NaN.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

from ..constitutive.base import RateModel
from ..constitutive.gnf import HERSCHEL_BULKLEY
from ..identify.fit import IdentifiedModel

__all__ = [
    "Needle", "Syringe", "InjectionForce", "ForceUQ",
    "wall_stress", "plug_radius", "wall_shear_rate",
    "flow_rate_hb", "flow_rate_gnf", "flow_rate", "pressure_for_flow_rate",
    "injection_force", "sensitivities", "force_with_uq",
]

_GD_BRACKET = (1e-12, 1e12)        # shear-rate bracket for inverting sigma(gd) = tau, 1/s
_DP_BRACKET = (1e-6, 1e12)         # pressure bracket for the flow-rate inversion, Pa
_SENS_REL_STEP = 1e-4              # central-difference step for d ln F / d ln p
_QUAD_LIMIT = 200                  # subdivision limit for the profile quadrature
_MONOTONE_RTOL = 1e-9              # relative decrease in sigma(gd) that counts as non-monotone
_DP_START = 1e3                    # Pa, first upper bracket tried by the pressure inversion


@dataclass(frozen=True)
class Needle:
    """Inner radius and length, metres. Inner diameter is what the datasheet gives, and
    a nominal gauge is not a radius: ISO 9626 allows a range per gauge, so the radius
    must come from the manufacturer's stated inner diameter, not the gauge number."""
    radius: float
    length: float

    def __post_init__(self):
        if not (self.radius > 0 and self.length > 0):
            raise ValueError(f"needle radius and length must be > 0, got {self}")


@dataclass(frozen=True)
class Syringe:
    """Barrel inner radius, metres. Force = pressure x plunger area."""
    barrel_radius: float

    def __post_init__(self):
        if not self.barrel_radius > 0:
            raise ValueError(f"barrel radius must be > 0, got {self.barrel_radius}")

    @property
    def area(self) -> float:
        return float(np.pi * self.barrel_radius ** 2)


# --- primitives --------------------------------------------------------------------

def wall_stress(dp: float, R: float, L: float) -> float:
    """``tau_w = dp R / (2 L)``, from the momentum balance alone (any fluid)."""
    return dp * R / (2.0 * L)


def plug_radius(dp: float, L: float, tau_y: float) -> float:
    """``r_p = 2 tau_y L / dp``; compare with ``R`` (>= R means no flow)."""
    return 2.0 * tau_y * L / dp


def flow_rate_hb(dp: float, R: float, L: float, tau_y: float, k: float, n: float) -> float:
    """Herschel-Bulkley flow rate (m^3/s) through a needle; closed form for every n.

    Returns 0.0 when ``tau_w <= tau_y``. Raises on non-positive geometry, pressure,
    ``k`` or ``n``; ``tau_y`` may be 0 (power law) and ``n`` may be 1 (Bingham).
    """
    if not (dp > 0 and R > 0 and L > 0 and k > 0 and n > 0 and tau_y >= 0):
        raise ValueError(f"flow_rate_hb: bad inputs dp={dp}, R={R}, L={L}, "
                         f"tau_y={tau_y}, k={k}, n={n}")
    tw = wall_stress(dp, R, L)
    S = tw - tau_y
    if S <= 0.0:
        return 0.0
    m = 1.0 / n
    bracket = (S ** (m + 3) / (m + 3) + 2.0 * tau_y * S ** (m + 2) / (m + 2)
               + tau_y ** 2 * S ** (m + 1) / (m + 1))
    return float(np.pi * R ** 3 / tw ** 3 * k ** (-m) * bracket)


def _shear_rate_at_stress(model: RateModel, params, tau: float) -> float:
    """Invert ``model.stress(gd) = tau`` for gd; 0 if the stress at the lower bracket
    already exceeds ``tau`` (below yield, or effectively so)."""
    lo, hi = _GD_BRACKET

    def f(log_gd):
        return float(model.stress(np.exp(log_gd), *params)) - tau

    if tau <= 0.0 or f(np.log(lo)) >= 0.0:
        return 0.0
    if f(np.log(hi)) < 0.0:
        raise ValueError(f"stress {tau:g} Pa not reached below gd = {hi:g} 1/s")
    return float(np.exp(brentq(f, np.log(lo), np.log(hi), xtol=1e-13, rtol=1e-13)))


def _assert_monotone_stress(model: RateModel, params, tau_w: float) -> None:
    """The inversion needs sigma(gd) increasing on the rates the needle imposes;
    a Cross model with m > 1 can violate that (a non-monotone flow curve)."""
    gd_w = _shear_rate_at_stress(model, params, tau_w)
    if gd_w <= 0.0:
        return
    gd = np.logspace(np.log10(_GD_BRACKET[0]), np.log10(gd_w), 400)
    sig = model.stress(gd, *params)
    # a genuine stress maximum is a decrease of order sig itself; a yield plateau at
    # tiny rates gives diffs that round to exactly 0, which is not a decrease
    if not np.all(np.diff(sig) >= -_MONOTONE_RTOL * sig[1:]):
        raise ValueError(f"{model.name}: stress is not monotone in rate up to gd_w = {gd_w:.3g}; "
                         "the needle-flow inversion is not defined for this flow curve")


def flow_rate_gnf(dp: float, R: float, L: float, model: RateModel, params) -> float:
    """Flow rate (m^3/s) for any registry :class:`RateModel` by quadrature of
    ``pi int_0^R r^2 gd(tau(r)) dr``, inverting the flow curve at each node.

    Independent of :func:`flow_rate_hb` (no closed form used); the two agree to
    1e-8 on HB in the tests. ``quad`` handles the integrable ``(r - r_p)^(1/n)``
    singularity at the plug edge.
    """
    if not (dp > 0 and R > 0 and L > 0):
        raise ValueError(f"flow_rate_gnf: bad geometry/pressure dp={dp}, R={R}, L={L}")
    params = np.asarray(params, dtype=float)
    model.check_params(params)
    tw = wall_stress(dp, R, L)
    _assert_monotone_stress(model, params, tw)

    def integrand(r):
        return r * r * _shear_rate_at_stress(model, params, tw * r / R)

    val, err = quad(integrand, 0.0, R, limit=_QUAD_LIMIT, epsabs=0.0, epsrel=1e-10)
    if val > 0 and err / val > 1e-6:
        raise RuntimeError(f"flow_rate_gnf: quadrature error {err / val:.2e} relative")
    return float(np.pi * val)


def flow_rate(dp: float, needle: Needle, model: RateModel, params) -> float:
    """Closed form when ``model`` is Herschel-Bulkley, quadrature otherwise."""
    if model is HERSCHEL_BULKLEY:
        return flow_rate_hb(dp, needle.radius, needle.length, *np.asarray(params, float))
    return flow_rate_gnf(dp, needle.radius, needle.length, model, params)


def pressure_for_flow_rate(Q: float, needle: Needle, model: RateModel, params) -> float:
    """Pressure drop (Pa) giving flow rate ``Q``; ``brentq`` in log(dp) on the
    monotone ``Q(dp)``. Raises if ``Q`` is outside the bracket ``_DP_BRACKET``."""
    if not Q > 0:
        raise ValueError(f"Q must be > 0, got {Q}")

    def f(log_dp):
        return flow_rate(np.exp(log_dp), needle, model, params) - Q

    # grow the upper bracket by decades from _DP_START rather than evaluating at an
    # absurd pressure, where the flow-curve inversion itself would run out of range
    lo, hi = _DP_BRACKET[0], _DP_START
    while f(np.log(hi)) < 0.0:
        lo, hi = hi, hi * 10.0
        if hi > _DP_BRACKET[1]:
            raise ValueError(f"Q = {Q:g} m^3/s needs dp > {_DP_BRACKET[1]:g} Pa")
    if not f(np.log(lo)) < 0.0:
        raise ValueError(f"Q = {Q:g} m^3/s is below the flow rate at dp = {lo:g} Pa")
    return float(np.exp(brentq(f, np.log(lo), np.log(hi), xtol=1e-12, rtol=1e-12)))


def wall_shear_rate(dp: float, needle: Needle, model: RateModel, params) -> float:
    """Shear rate at the needle wall, the largest rate the fluid sees."""
    return _shear_rate_at_stress(model, np.asarray(params, float),
                                 wall_stress(dp, needle.radius, needle.length))


# --- force -------------------------------------------------------------------------

@dataclass(frozen=True)
class InjectionForce:
    force: float            # N, on the plunger
    dp: float               # Pa, across the needle
    wall_stress: float      # Pa
    wall_shear_rate: float  # 1/s
    plug_fraction: float    # r_p / R for a yield-stress model, 0 otherwise
    reynolds: float         # rho U D / eta(gd_w); NaN if rho not given
    friction: float         # N, as supplied

    def __str__(self) -> str:
        return (f"F = {self.force:.3f} N (dp = {self.dp / 1e5:.3f} bar, tau_w = {self.wall_stress:.1f} Pa, "
                f"gd_w = {self.wall_shear_rate:.3g} 1/s, plug r_p/R = {self.plug_fraction:.3f}, "
                f"Re = {self.reynolds:.3g})")


def injection_force(Q: float, needle: Needle, syringe: Syringe, model: RateModel, params,
                    *, friction: float = 0.0, rho: float | None = None) -> InjectionForce:
    """Plunger force for flow rate ``Q`` (m^3/s): needle pressure drop times plunger
    area, plus ``friction`` (N). See the module docstring for what is left out."""
    params = np.asarray(params, dtype=float)
    dp = pressure_for_flow_rate(Q, needle, model, params)
    tw = wall_stress(dp, needle.radius, needle.length)
    gd_w = _shear_rate_at_stress(model, params, tw)
    pd = model.params_dict(params)
    plug = plug_radius(dp, needle.length, pd["tau_y"]) / needle.radius if "tau_y" in pd else 0.0
    if rho is None:
        re = float("nan")
    else:
        U = Q / (np.pi * needle.radius ** 2)
        eta_w = float(model(gd_w, *params)) if gd_w > 0 else float("inf")
        re = float(rho * U * 2 * needle.radius / eta_w)
    return InjectionForce(force=float(dp * syringe.area + friction), dp=dp, wall_stress=tw,
                          wall_shear_rate=gd_w, plug_fraction=float(plug), reynolds=re,
                          friction=friction)


def sensitivities(Q: float, needle: Needle, syringe: Syringe, model: RateModel, params,
                  *, friction: float = 0.0) -> dict[str, float]:
    """``d ln F / d ln p_j`` at ``params`` by central differences (relative step 1e-4).
    A value near 0 says the force does not care about that parameter at this
    operating point; the report states these explicitly."""
    params = np.asarray(params, dtype=float)
    out = {}
    for j, name in enumerate(model.param_names):
        if params[j] == 0.0:
            out[name] = 0.0  # a parameter at exactly 0 has no log-sensitivity; report inert
            continue
        h = _SENS_REL_STEP * abs(params[j])
        up, dn = params.copy(), params.copy()
        up[j] += h
        dn[j] -= h
        f_up = injection_force(Q, needle, syringe, model, up, friction=friction).force
        f_dn = injection_force(Q, needle, syringe, model, dn, friction=friction).force
        f_0 = injection_force(Q, needle, syringe, model, params, friction=friction).force
        out[name] = float((f_up - f_dn) / (2 * h) * params[j] / f_0)
    return out


@dataclass(frozen=True)
class ForceUQ:
    nominal: InjectionForce
    level: float
    profile_lo: float               # N, min force over one-at-a-time profile pushes
    profile_hi: float               # N, max
    profile_open: tuple[str, ...]   # parameters whose profile interval hit a bound
    per_parameter: dict[str, tuple[float, float]]  # F at (lo, hi) of each parameter
    sampled_lo: float               # N, quantile of the covariance sample; NaN if singular
    sampled_hi: float
    n_samples: int
    n_rejected: int                 # draws outside the model bounds, discarded
    seed: int
    sensitivity: dict[str, float]   # d ln F / d ln p at the fit
    insensitive_to: tuple[str, ...]  # |sensitivity| < INSENSITIVE
    extrapolated: bool              # wall shear rate outside the envelope, at the fit or any profile end
    envelope: tuple[float, float]   # fitted (gammadot_min, gammadot_max)
    notes: tuple[str, ...] = field(default_factory=tuple)

    def summary(self) -> str:
        n = self.nominal
        lines = [str(n),
                 f"profile interval ({self.level:.0%}, one-at-a-time): [{self.profile_lo:.3f}, {self.profile_hi:.3f}] N"
                 + (f"  OPEN in {', '.join(self.profile_open)}" if self.profile_open else ""),
                 f"covariance sample ({self.n_samples} draws, {self.n_rejected} rejected, seed {self.seed}): "
                 f"[{self.sampled_lo:.3f}, {self.sampled_hi:.3f}] N",
                 "sensitivity d ln F / d ln p: " + ", ".join(f"{k} {v:+.3f}" for k, v in self.sensitivity.items()),
                 f"insensitive to: {', '.join(self.insensitive_to) or 'none'}",
                 f"wall shear rate {n.wall_shear_rate:.3g} 1/s vs fitted envelope "
                 f"[{self.envelope[0]:.3g}, {self.envelope[1]:.3g}]: "
                 + ("EXTRAPOLATED" if self.extrapolated else "inside")]
        lines += list(self.notes)
        return "\n".join(lines)


INSENSITIVE = 0.05  # |d ln F / d ln p| below which the report calls the force insensitive


def force_with_uq(identified: IdentifiedModel, Q: float, needle: Needle, syringe: Syringe,
                  *, friction: float = 0.0, rho: float | None = None,
                  n_samples: int = 2000, seed: int = 0) -> ForceUQ:
    """Force with both uncertainty routes described in the module docstring."""
    model, values = identified.model, identified.values
    nominal = injection_force(Q, needle, syringe, model, values, friction=friction, rho=rho)

    # (1) profile push-through, one parameter at a time
    per_param: dict[str, tuple[float, float]] = {}
    forces = [nominal.force]
    open_params = []
    for j, name in enumerate(model.param_names):
        est = identified.params[name]
        pair = []
        for edge in (est.lo, est.hi):
            p = values.copy()
            p[j] = edge
            pair.append(injection_force(Q, needle, syringe, model, p, friction=friction).force)
        per_param[name] = (pair[0], pair[1])
        forces += pair
        if est.lo_at_bound or est.hi_at_bound:
            open_params.append(name)

    # (2) covariance sampling with rejection to the model bounds
    cov = np.asarray(identified.fit.cov, dtype=float)
    notes = []
    if cov.shape == (values.size, values.size) and np.all(np.isfinite(cov)):
        rng = np.random.default_rng(seed)
        draws = rng.multivariate_normal(values, cov, size=n_samples)
        lo_b, hi_b = np.array(model.bounds[0], float), np.array(model.bounds[1], float)
        keep = np.all((draws > lo_b) & (draws < hi_b), axis=1)
        n_rej = int(n_samples - keep.sum())
        f = np.array([injection_force(Q, needle, syringe, model, d, friction=friction).force
                      for d in draws[keep]])
        if f.size < 10:
            s_lo = s_hi = float("nan")
            notes.append(f"covariance sample unusable: {f.size} of {n_samples} draws inside the bounds")
        else:
            a = (1.0 - identified.level) / 2.0
            s_lo, s_hi = (float(x) for x in np.quantile(f, [a, 1.0 - a]))
            if n_rej > 0.1 * n_samples:
                notes.append(f"{n_rej / n_samples:.0%} of covariance draws fell outside the bounds: "
                             "the Wald normal is a poor description of this posterior")
    else:
        n_rej = 0
        s_lo = s_hi = float("nan")
        notes.append("covariance is singular or undefined (degenerate fit): sampled interval not available")

    sens = sensitivities(Q, needle, syringe, model, values, friction=friction)
    insens = tuple(k for k, v in sens.items() if abs(v) < INSENSITIVE)
    env = (float(identified.envelope.gammadot_min), float(identified.envelope.gammadot_max))
    # judged at the nominal fit AND at every one-at-a-time profile end: an interval whose
    # far end sits outside the envelope is partly extrapolation even if the fit is not
    wall_rates = [nominal.wall_shear_rate]
    for j, name in enumerate(model.param_names):
        est = identified.params[name]
        for edge in (est.lo, est.hi):
            p = values.copy()
            p[j] = edge
            wall_rates.append(injection_force(Q, needle, syringe, model, p, friction=friction).wall_shear_rate)
    extrap = any(not (env[0] <= w <= env[1]) for w in wall_rates if np.isfinite(w))
    return ForceUQ(nominal=nominal, level=identified.level,
                   profile_lo=float(min(forces)), profile_hi=float(max(forces)),
                   profile_open=tuple(open_params), per_parameter=per_param,
                   sampled_lo=s_lo, sampled_hi=s_hi, n_samples=n_samples, n_rejected=n_rej,
                   seed=seed, sensitivity=sens, insensitive_to=insens,
                   extrapolated=extrap, envelope=env, notes=tuple(notes))
