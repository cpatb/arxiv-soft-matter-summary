"""Constitutive models, all written as viscosity as a function of *stress*.

Why stress and not shear rate
-----------------------------
For a discontinuously shear-thickening suspension the flow curve sigma(gamma_dot)
is S-shaped and multivalued: at one shear rate there can be three stresses. The
inverse, gamma_dot(sigma), is single-valued. Wyart-Cates is written stress-
controlled for exactly this reason, and forcing it into rate space would break
it on the very case the package exists to fit.

So every model here is ``eta(sigma; params)``, and the measured stress is taken
from ``FlowSweep.stress`` (recomputed as eta * gamma_dot, an identity for steady
simple shear). The thinning baselines are given in their stress forms, Ellis
rather than Cross, so that no model needs a root find.

What each model can and cannot do
---------------------------------
Herschel-Bulkley and Ellis describe yielding and thinning. Wyart-Cates describes
thickening and jamming and has **no thinning branch at all**: its viscosity is
monotonically increasing in stress. Real corn starch thins and then thickens, so
fitting Wyart-Cates to a whole sweep is a category error. Fit it to the
thickening branch and say where the branch was cut, because that cut moves
sigma*. :func:`thickening_branch` picks the cut and returns the index so it can
be reported rather than buried.

Fitted parameters versus fixed constants
----------------------------------------
Wyart-Cates as written here has three *fitted* parameters (phi, phi_m, sigma*)
and four *fixed* constants (eta_s, phi_0, n, beta). The constants are properties
of the suspending fluid, the packing and the chosen functional form, not of the
flow curve, so they are supplied rather than inferred; freeing eta_s alongside
phi and phi_m is one of the degeneracies N3 exists to catch.

Fixed does not mean universal. Every one of the four is material- and
solvent-specific, and until 2026-08-27 all four were module-level constants
carrying values appropriate to corn starch in water. Applying those to a
published PMMA-in-organic-solvent dataset silently mis-states the answer, and
the error lands in phi_m, which is the number the package is for. So they are
now bound per dataset with :meth:`Model.with_fixed`, and the module defaults are
unchanged only so that results obtained before that change still reproduce::

    guy = WYART_CATES.with_fixed(eta_s=2.4e-3, phi_0=0.611, n=2.2)
    fit_model(guy, stress, viscosity, branch=br)

``with_fixed`` refuses a constant the model's function does not accept, so a
typo raises instead of being ignored.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass, field, replace
from typing import Callable, Mapping, Sequence

import numpy as np

# Solvent viscosity of water at 25 degC, Pa.s. The *default* only; see the module
# docstring. Guy et al. 2015 is 2.4e-3, Richards et al. 2019 is 6.0e-3.
ETA_S_WATER = 8.9e-4

# Random close packing of frictionless monodisperse spheres. Again a default:
# the frictionless divergence of Guy et al.'s d = 3770 nm series sits at 0.611,
# and Richards et al. use 0.546 for corn starch in weight-fraction units.
PHI_0_RCP = 0.64

# Exponent of the divergence, eta ~ (1 - phi/phi_J)^-n. n = 2 is the
# Maron-Pierce / Krieger-Dougherty value and what most closures assume. Guy et
# al. fit n = 2.2(2) to their own upper branch, and a phi_m measured at one
# exponent cannot be compared with a phi_m measured at another.
WC_EXPONENT = 2.0

# Sharpness of the stress activation, f = exp(-(sigma*/sigma)^beta). beta = 1 is
# the original Wyart-Cates form; Richards et al. fit 0.94 to corn starch.
WC_BETA = 1.0

# phi and phi_m must both stay below the frictionless packing phi_0, otherwise
# "thickening" runs backwards. Bounds are held this far below it so the fitter
# never evaluates the divergence exactly at the pole.
PHI_BOUND_MARGIN = 1e-3


@dataclass(frozen=True, eq=False)
class Model:
    """A constitutive model expressed as eta(sigma).

    ``func`` takes ``(sigma, *params, **fixed)`` in SI and returns viscosity in
    Pa.s. ``bounds`` is ``(lower, upper)`` per *fitted* parameter, used by the
    fitter and by the synthetic-data generator. ``fixed`` holds the constants
    that are supplied rather than fitted; see :meth:`with_fixed`.

    ``eq=False`` because ``fixed`` is a dict: the default frozen-dataclass
    ``__hash__`` would raise on it. Nothing compares models, and identity
    semantics are the honest default for something holding a closure.
    """

    name: str
    param_names: tuple[str, ...]
    func: Callable[..., np.ndarray]
    bounds: tuple[tuple[float, ...], tuple[float, ...]]
    describes_thickening: bool
    notes: str = ""
    fixed: Mapping[str, float] = field(default_factory=dict)

    def __call__(self, sigma: np.ndarray, *params: float) -> np.ndarray:
        return self.func(np.asarray(sigma, dtype=float), *params, **self.fixed)

    @property
    def n_params(self) -> int:
        return len(self.param_names)

    def check_params(self, params: Sequence[float]) -> None:
        lo, hi = self.bounds
        if len(params) != self.n_params:
            raise ValueError(
                f"{self.name}: expected {self.n_params} params "
                f"{self.param_names}, got {len(params)}"
            )
        for name, p, a, b in zip(self.param_names, params, lo, hi):
            if not (a <= p <= b):
                raise ValueError(f"{self.name}: {name}={p:g} outside bounds [{a:g}, {b:g}]")

    def accepts_fixed(self) -> tuple[str, ...]:
        """Constant names ``func`` will accept, i.e. its keyword-only arguments."""
        sig = inspect.signature(self.func)
        return tuple(
            n for n, p in sig.parameters.items()
            if p.kind is inspect.Parameter.KEYWORD_ONLY
        )

    def with_fixed(self, **overrides: float) -> "Model":
        """Bind non-fitted constants for one dataset, returning a new model.

        The name is left unchanged so the fitter's per-model initial guess and
        every existing report keep working; what changed travels in ``fixed``
        and is printed by :meth:`describe_fixed`.

        Setting ``phi_0`` also moves the upper bounds on ``phi`` and ``phi_m``,
        because the model requires both to sit below the frictionless packing
        fraction. That is not a tuning knob: leaving the bounds at their
        water/RCP values while lowering phi_0 would let the fitter return
        phi_m > phi_0, which inverts the sign of the thickening.

        Raises
        ------
        ValueError
            If a name is not a keyword-only argument of ``func``, or if a
            supplied value is not a positive finite number.
        """
        allowed = self.accepts_fixed()
        unknown = sorted(set(overrides) - set(allowed))
        if unknown:
            raise ValueError(
                f"{self.name}: cannot fix {unknown}; it takes "
                f"{list(allowed) if allowed else 'no fixed constants'}"
            )
        for k, v in overrides.items():
            if not np.isfinite(v) or v <= 0:
                raise ValueError(
                    f"{self.name}: fixed constant {k}={v!r} must be finite and > 0"
                )

        new_fixed = {**self.fixed, **overrides}
        lo = np.array(self.bounds[0], dtype=float)
        hi = np.array(self.bounds[1], dtype=float)
        if "phi_0" in overrides:
            cap = float(overrides["phi_0"]) - PHI_BOUND_MARGIN
            for pname in ("phi", "phi_m"):
                if pname in self.param_names:
                    i = self.param_names.index(pname)
                    if lo[i] >= cap:
                        raise ValueError(
                            f"{self.name}: phi_0={overrides['phi_0']:g} leaves no room "
                            f"for {pname} above its lower bound {lo[i]:g}"
                        )
                    hi[i] = cap
        return replace(
            self,
            fixed=new_fixed,
            bounds=(tuple(lo.tolist()), tuple(hi.tolist())),
        )

    def describe_fixed(self) -> str:
        """The constants in force, for a report header. Defaults are shown too."""
        defaults = {
            n: p.default
            for n, p in inspect.signature(self.func).parameters.items()
            if p.kind is inspect.Parameter.KEYWORD_ONLY
        }
        if not defaults:
            return "no fixed constants"
        parts = []
        for n, d in defaults.items():
            v = self.fixed.get(n, d)
            parts.append(f"{n}={v:g}" + ("" if n in self.fixed else " (default)"))
        return ", ".join(parts)


# --------------------------------------------------------------------------
# thinning / yielding baselines
# --------------------------------------------------------------------------
def _herschel_bulkley(sigma, tau_y, k, n):
    """eta(sigma) for Herschel-Bulkley, exact inverse of sigma = tau_y + k*gdot^n.

    Below the yield stress the material does not flow, so the viscosity is
    unbounded. Returned as inf rather than a large finite number: a regularized
    value here would silently become a fitting parameter. Regularization is a
    solver concern and belongs at the solver boundary, not in the constitutive
    definition.
    """
    sigma = np.asarray(sigma, dtype=float)
    out = np.full(sigma.shape, np.inf)
    flowing = sigma > tau_y
    excess = np.maximum(sigma[flowing] - tau_y, 0.0)
    gdot = (excess / k) ** (1.0 / n)
    with np.errstate(divide="ignore", invalid="ignore"):
        out[flowing] = np.where(gdot > 0, sigma[flowing] / gdot, np.inf)
    return out


def _ellis(sigma, eta_0, eta_inf, sigma_c, a):
    """Ellis model, the stress-space form of Cross. Purely thinning."""
    sigma = np.asarray(sigma, dtype=float)
    return eta_inf + (eta_0 - eta_inf) / (1.0 + (sigma / sigma_c) ** a)


# --------------------------------------------------------------------------
# thickening
# --------------------------------------------------------------------------
def frictional_fraction(sigma, sigma_star, beta=WC_BETA):
    """Wyart-Cates stress-activated fraction of frictional contacts.

    f -> 0 as sigma -> 0 (lubricated), f -> 1 as sigma -> inf (frictional).
    ``beta`` sets how sharply the crossover happens; 1 is the original form.
    """
    sigma = np.asarray(sigma, dtype=float)
    with np.errstate(divide="ignore", over="ignore"):
        return np.exp(-((sigma_star / np.maximum(sigma, 1e-300)) ** beta))


def jamming_fraction(f, phi_0=PHI_0_RCP, phi_m=0.58):
    """Jamming packing fraction interpolating frictionless phi_0 to frictional phi_m."""
    return (1.0 - f) * phi_0 + f * phi_m


def divergence(phi, phi_j, eta_s=ETA_S_WATER, n=WC_EXPONENT):
    """eta = eta_s * (1 - phi/phi_j)^-n, with inf at and beyond the pole.

    The divergence is returned as inf, not clipped: clipping it would turn a
    jammed state into a large-but-finite viscosity that a least-squares fit
    would happily chase.
    """
    ratio = np.asarray(np.asarray(phi, dtype=float) / np.asarray(phi_j, dtype=float))
    out = np.full(ratio.shape, np.inf, dtype=float)
    free = ratio < 1.0
    out[free] = eta_s * (1.0 - ratio[free]) ** -float(n)
    return out


def _wyart_cates(
    sigma,
    phi,
    phi_m,
    sigma_star,
    *,
    eta_s=ETA_S_WATER,
    phi_0=PHI_0_RCP,
    n=WC_EXPONENT,
    beta=WC_BETA,
):
    """Wyart-Cates viscosity, eta = eta_s * (1 - phi/phi_J(sigma))^-n.

    Diverges where phi reaches phi_J, which is the shear-jamming condition.
    ``eta_s``, ``phi_0``, ``n`` and ``beta`` are fixed constants of the material
    and the chosen form, not fitted; bind them with ``Model.with_fixed``.
    """
    f = frictional_fraction(sigma, sigma_star, beta=beta)
    phi_j = jamming_fraction(f, phi_0=phi_0, phi_m=phi_m)
    phi_full = np.broadcast_to(float(phi), np.shape(phi_j))
    return divergence(phi_full, phi_j, eta_s=eta_s, n=n)


def _wyart_cates_known_phi(
    sigma,
    phi_m,
    sigma_star,
    *,
    phi,
    eta_s=ETA_S_WATER,
    phi_0=PHI_0_RCP,
    n=WC_EXPONENT,
    beta=WC_BETA,
):
    """Wyart-Cates with the volume fraction supplied instead of fitted.

    Same physics as :func:`_wyart_cates`; the difference is entirely in what is
    being asked of the data. ``phi`` and ``phi_m`` enter the model only through
    the ratio ``phi/phi_J``, and near jamming a flow curve constrains that ratio
    far better than it constrains either number alone. Fitting both from one
    curve therefore returns a ``phi_m`` that tracks the fitted ``phi`` rather
    than a property of the particles.

    When the volume fraction is known independently, which for a made-up
    suspension it is, pinning it removes the trade-off and leaves ``phi_m``
    actually measured. ``phi`` has no default: a wrong one is worse than none,
    so it must be supplied.
    """
    f = frictional_fraction(sigma, sigma_star, beta=beta)
    phi_j = jamming_fraction(f, phi_0=phi_0, phi_m=phi_m)
    return divergence(np.broadcast_to(float(phi), np.shape(phi_j)), phi_j, eta_s=eta_s, n=n)


def _maron_pierce(sigma, phi, phi_m, *, eta_s=ETA_S_WATER, n=WC_EXPONENT):
    """Rate-independent Krieger-Dougherty/Maron-Pierce plateau. No stress dependence.

    Included as the null model: if a sweep is fit as well by this as by
    Wyart-Cates, the data contains no evidence of stress-activated thickening.
    """
    val = eta_s * (1.0 - phi / phi_m) ** -float(n) if phi < phi_m else np.inf
    return np.full(np.shape(sigma), val, dtype=float)


# --------------------------------------------------------------------------
# registry
# --------------------------------------------------------------------------
HERSCHEL_BULKLEY = Model(
    name="herschel_bulkley",
    param_names=("tau_y", "k", "n"),
    func=_herschel_bulkley,
    bounds=((0.0, 1e-6, 0.05), (1e5, 1e6, 2.0)),
    describes_thickening=False,
    notes="yield stress plus power law; n > 1 can mimic thickening but has no phi",
)

ELLIS = Model(
    name="ellis",
    param_names=("eta_0", "eta_inf", "sigma_c", "a"),
    func=_ellis,
    bounds=((1e-6, 1e-9, 1e-6, 0.05), (1e6, 1e3, 1e6, 10.0)),
    describes_thickening=False,
    notes="stress-space Cross; thinning only",
)

_PHI_CAP = PHI_0_RCP - PHI_BOUND_MARGIN  # 0.639, unchanged from before the refactor

WYART_CATES = Model(
    name="wyart_cates",
    param_names=("phi", "phi_m", "sigma_star"),
    func=_wyart_cates,
    bounds=((0.10, 0.30, 1e-4), (_PHI_CAP, _PHI_CAP, 1e5)),
    describes_thickening=True,
    notes="eta_s, phi_0, n and beta are fixed constants; bind them with with_fixed()",
)

WYART_CATES_KNOWN_PHI = Model(
    name="wyart_cates_known_phi",
    param_names=("phi_m", "sigma_star"),
    func=_wyart_cates_known_phi,
    bounds=((0.30, 1e-4), (_PHI_CAP, 1e5)),
    describes_thickening=True,
    notes="phi supplied, not fitted; with_fixed(phi=...) is required before use",
)

MARON_PIERCE = Model(
    name="maron_pierce",
    param_names=("phi", "phi_m"),
    func=_maron_pierce,
    bounds=((0.10, 0.30), (_PHI_CAP, 0.70)),
    describes_thickening=False,
    notes="null model: no stress dependence at all",
)

MODELS: dict[str, Model] = {
    m.name: m
    for m in (HERSCHEL_BULKLEY, ELLIS, WYART_CATES, WYART_CATES_KNOWN_PHI, MARON_PIERCE)
}


def get_model(name: str) -> Model:
    try:
        return MODELS[name]
    except KeyError:
        raise KeyError(f"unknown model {name!r}; available: {sorted(MODELS)}") from None


# --------------------------------------------------------------------------
# branch selection
# --------------------------------------------------------------------------
def thickening_branch(viscosity: np.ndarray, pad: int = 0) -> slice:
    """Indices from the viscosity minimum onward, i.e. the thickening branch.

    Corn starch thins before it thickens and Wyart-Cates has no thinning branch,
    so a whole-sweep fit is a category error. Where this cut falls changes the
    recovered sigma*, so the fitter records it and the caller should report it.

    ``pad`` moves the cut later (positive) or earlier (negative) so the
    sensitivity of a result to the cut can be measured rather than assumed.
    """
    eta = np.asarray(viscosity, dtype=float)
    if eta.size == 0:
        raise ValueError("empty viscosity array")
    start = int(np.argmin(eta)) + pad
    start = max(0, min(start, eta.size - 1))
    return slice(start, eta.size)
