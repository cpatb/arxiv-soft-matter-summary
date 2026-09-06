"""Fit a rate model to a flow curve and say, per parameter, whether it was measured.

The verdict rules are fixed here, before any data is looked at, and applied in code
(the ``enforce-gates-in-code`` lesson: a rule applied by hand at read time gets
skipped). They are read off ``dst_fit.profiling.profile_likelihood``, which walks the
profile likelihood outward from the best fit until the F statistic crosses its
critical value or the parameter's bound is hit:

measurable
    the profile interval closes on both sides, and neither end sits on a bound;
bounded_only
    the interval closes on one side only (the other end ran into a bound), so the data
    gives a one-sided bound and not a value;
not_identifiable
    both ends are at bounds: the profile never crossed the critical value, the data
    says nothing about this parameter on its own;
degenerate
    at least one end is open AND the parameter is traded against another one (a ratio
    or product enters the model, the ``parameters-that-enter-as-a-ratio`` lesson),
    detected by either of two signatures: the fit correlation has ``|corr| > 0.99``,
    or the log-log sensitivity Jacobian ``d ln eta / d ln p`` at the fit has a null
    direction (``s_min/s_max < DEGENERATE_COND``) in which this parameter and a
    partner both carry weight (``|v_j| > NULL_WEIGHT``). Checked before
    ``bounded_only``/``not_identifiable``, so a degenerate pair is named as such
    rather than reported as two unrelated bounds.

    The second signature was added 2026-09-05 after the structural gate
    (``tests/test_structural.py``) missed its pre-registered bar at 16/20: on 4 of 20
    seeds the fitter sat exactly on the symmetric point ``alpha = t_r``, where the
    Jacobian is rank-deficient, the covariance is all-NaN, and the corr rule read the
    strongest possible degeneracy signal as no evidence. The bar was not moved; the
    detector got a second, independent channel that a singular covariance cannot hide
    from.

"Bounded" here means the likelihood-ratio interval at ``level``, not the Wald error;
``dst_fit`` reports both, and the profile is the one that is honest near a bound. A
``measurable`` verdict is a statement about closure, not precision; ``ParamEstimate.decades``
carries the width and must be quoted beside it.

**CAVEAT ON THE CORRELATION SIGNATURE, and it is a real one.** ``|corr| > 0.99`` is read
off the Wald covariance -- a local quadratic approximation to the likelihood at the
optimum -- and it is consulted ONLY when ``open_side`` is true, i.e. only when the profile
ran into a bound. That is exactly the regime where the quadratic approximation is least
trustworthy: at or near a bound the likelihood is not locally quadratic, the covariance
linearises a problem the optimiser solved under a constraint, and the curvature the
correlation is computed from need not describe the surface the profile actually walked.
So the first signature is being used where it is weakest, by construction.

Three things follow, and they are why the code is shaped as it is rather than being a
reason to distrust the verdicts:

* the SECOND signature (a null direction of the log-log sensitivity Jacobian) does not go
  through the covariance at all, so it survives the case where the covariance is singular
  or all-NaN -- which is precisely how it was found (see above);
* the two signatures are ORed, never ANDed, so a failure of the Wald channel can only
  cost a ``degenerate`` label, never manufacture one out of nothing;
* ``degenerate`` and ``bounded_only`` and ``not_identifiable`` all refuse an export
  the same way, so a mislabel among the three does not change what ships.

What the caveat DOES forbid is quoting ``ParamEstimate.corr_with`` as a measured
correlation between two parameters when either sits at a bound. It is a detector input
there, not an estimate. Measured 2026-09-06 (``data/horner2020/multistart_scatter.md``):
the verdicts nevertheless predict which records move when only the optimiser's seed
changes, on 20 donors x 4 families, with no all-``measurable`` record moving more than
1e-3 -- an independent check that does not go through the covariance either.

Every verdict above is conditional on the model describing the curve. A parameter can
be ``measurable`` in a model that is wrong: the interval then says how precisely the
wrong model's constant is pinned, not how well the fluid is described. So the fit
also reports the STRUCTURE of its log residuals ordered by shear rate, as a
Wald-Wolfowitz runs test on their signs (``ResidualStructure``): a model that misses
the curve's shape leaves long runs of one sign, far fewer runs than the
``1 + 2 n+ n- / n`` a random sequence has. ``structured`` is set when the standard
score is below ``RUNS_Z = -2.0`` (one-sided, fixed 2026-09-05 before the blood
records were fitted). Hand anchor: 16 points with signs ``--+++++++++-----`` have 3
runs where 8.875 are expected, z = -3.09. The flag is a statement about the model
form and is written into the sidecar beside the verdicts; it changes no verdict.

Read it WITH ``rmse_log``, never alone: it measures whether the misfit is systematic,
not how large it is. On the 15 Carbopol ETD2020 records (54 points each) it fires on
13, at rmse_log from 0.002 (0.5 wt%, z = -3.6: a misfit of 0.2 % that a 54-point
series still resolves as non-random) to 0.17 (0.1 wt%, z = -6.6, signs
``++++++++-----...+++++-------``: Herschel-Bulkley is the wrong form for that curve
within the data). Synthetic curves drawn from the fitted model trip it at the
expected rate (``tests/test_identify.py``).

This module is the only writer of a record's ``identifiability`` block.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from dst_fit import FitResult, fit_model
from dst_fit.profiling import Profile, profile_likelihood

from ..constitutive import RateModel
from ..registry.schema import Envelope

__all__ = ["VERDICTS", "ParamEstimate", "IdentifiedModel", "ResidualStructure",
           "fit_rate_curve", "verdict_for", "null_direction", "runs_test"]

VERDICTS = ("measurable", "bounded_only", "degenerate", "not_identifiable")
DEGENERATE_CORR = 0.99  # |corr| above which two parameters count as one
DEGENERATE_COND = 1e-8  # s_min/s_max of the log-log Jacobian below which a null direction exists
NULL_WEIGHT = 0.3       # |component| in the null vector for a parameter to "carry weight"
_FD_REL_STEP = 1e-6     # central-difference step, relative to |p|
RUNS_Z = -2.0           # runs-test score below which the residuals are called structured


@dataclass(frozen=True)
class ParamEstimate:
    name: str
    value: float
    lo: float  # profile-likelihood interval; equals the bound where at_bound
    hi: float
    verdict: str
    lo_at_bound: bool
    hi_at_bound: bool
    wald_stderr: float
    partner: str | None = None  # the parameter it is degenerate with, if any

    def __post_init__(self):
        if self.verdict not in VERDICTS:
            raise ValueError(f"verdict {self.verdict!r} not in {VERDICTS}")

    @property
    def decades(self) -> float:
        """Width of the interval in decades, ``log10(hi/lo)``; inf if ``lo <= 0``.

        "measurable" says the interval CLOSES, not that it is narrow: on a curve with
        no yield plateau the fitter has returned a closed tau_y interval spanning
        [0.002, 9] Pa (3.7 decades) around a truth of 1 Pa. Quote the width with the
        verdict.
        """
        if self.lo <= 0 or not np.isfinite(self.hi):
            return float("inf")
        return float(np.log10(self.hi / self.lo))


@dataclass(frozen=True)
class ResidualStructure:
    """Wald-Wolfowitz runs test on the signs of the log residuals ordered by rate."""
    n: int
    n_runs: int
    expected_runs: float
    z: float
    structured: bool
    signs: str  # e.g. "--+++++++++-----", data - model, in order of increasing rate

    def as_dict(self) -> dict:
        return {"n": self.n, "n_runs": self.n_runs, "expected_runs": self.expected_runs,
                "z": self.z, "structured": self.structured, "signs": self.signs}


def runs_test(residuals: np.ndarray) -> ResidualStructure:
    """Runs test on ``sign(residuals)`` in the order given (caller sorts by rate).

    Under the null of exchangeable signs, with ``n1`` positive and ``n2`` negative,
    the number of runs has mean ``1 + 2 n1 n2 / n`` and variance
    ``2 n1 n2 (2 n1 n2 - n) / (n^2 (n - 1))`` (Wald & Wolfowitz 1940). Exact zeros are
    dropped rather than assigned a sign. If either sign is absent the score is
    undefined (nan) and ``structured`` is False: one-signed residuals are a bias,
    which a runs test cannot see, and rmse carries it.
    """
    r = np.asarray(residuals, dtype=float)
    r = r[np.isfinite(r) & (r != 0)]
    n = int(r.size)
    signs = "".join("+" if x > 0 else "-" for x in r)
    n1 = int((r > 0).sum())
    n2 = n - n1
    n_runs = 1 + int((np.sign(r[1:]) != np.sign(r[:-1])).sum()) if n else 0
    if n1 == 0 or n2 == 0 or n < 2:
        return ResidualStructure(n, n_runs, float("nan"), float("nan"), False, signs)
    mu = 1.0 + 2.0 * n1 * n2 / n
    var = 2.0 * n1 * n2 * (2.0 * n1 * n2 - n) / (n * n * (n - 1.0))
    z = (n_runs - mu) / np.sqrt(var) if var > 0 else float("nan")
    return ResidualStructure(n, n_runs, float(mu), float(z), bool(np.isfinite(z) and z < RUNS_Z), signs)


@dataclass(frozen=True)
class IdentifiedModel:
    model: RateModel
    params: dict[str, ParamEstimate]
    fit: FitResult
    envelope: Envelope
    level: float
    residual_structure: ResidualStructure | None = None

    @property
    def values(self) -> np.ndarray:
        return np.array([self.params[p].value for p in self.model.param_names])

    @property
    def verdicts(self) -> dict[str, str]:
        return {k: v.verdict for k, v in self.params.items()}

    def identifiability_block(self) -> dict[str, str]:
        """What the registry stores under ``identifiability``; nothing else writes it."""
        return {k: (v.verdict if v.partner is None else f"{v.verdict} (with {v.partner})")
                for k, v in self.params.items()}


def null_direction(model: RateModel, gd: np.ndarray, params: np.ndarray
                   ) -> tuple[float, np.ndarray]:
    """``(s_min/s_max, v)`` of the log-log sensitivity Jacobian ``d ln eta_i / d ln p_j``
    at ``params``, with ``v`` the unit right-singular vector of the smallest singular
    value. A ratio near machine zero means the data cannot move along ``v`` at all;
    the parameters with large ``|v_j|`` are the ones traded against each other.

    Central differences with a relative step; a parameter at exactly 0 is stepped
    absolutely (its log-sensitivity is then a plain sensitivity, which is fine for
    the rank question). If one side of the step leaves the model's domain (eta <= 0,
    typically a parameter sitting on its lower bound 0) the other side is used alone,
    so the channel is not switched off by a bound. If the Jacobian is still
    non-finite the result is ``(nan, nan-vector)`` and ``verdict_for`` falls back to
    the correlation rule.
    """
    p = np.asarray(params, dtype=float)
    cols = []
    for j in range(p.size):
        scale = abs(p[j]) if p[j] != 0 else 1.0
        h = _FD_REL_STEP * scale
        up, dn = p.copy(), p.copy()
        up[j] += h
        dn[j] -= h
        with np.errstate(all="ignore"):
            f_up, f_dn = np.log(model(gd, *up)), np.log(model(gd, *dn))
            f_0 = None
            if not np.all(np.isfinite(f_dn)):   # e.g. p_j = 0 stepped negative
                f_0 = np.log(model(gd, *p))
                d = (f_up - f_0) / h
            elif not np.all(np.isfinite(f_up)):
                f_0 = np.log(model(gd, *p))
                d = (f_0 - f_dn) / h
            else:
                d = (f_up - f_dn) / (2 * h)
        cols.append(d * scale)
    J = np.column_stack(cols)
    if not np.all(np.isfinite(J)):
        return float("nan"), np.full(p.size, np.nan)
    _, s, vt = np.linalg.svd(J, full_matrices=False)
    return float(s[-1] / s[0]) if s[0] > 0 else 0.0, vt[-1]


def verdict_for(profile: Profile, corr_row: np.ndarray, index: int,
                param_names: tuple[str, ...], *, cond: float = np.inf,
                null_vec: np.ndarray | None = None) -> tuple[str, str | None]:
    """Apply the rules in the module docstring to one parameter's profile."""
    open_side = profile.lower_at_bound or profile.upper_at_bound
    if open_side:
        others = [(abs(float(corr_row[j])), param_names[j])
                  for j in range(len(param_names)) if j != index]
        if others:
            c, partner = max(others)
            if np.isfinite(c) and c > DEGENERATE_CORR:
                return "degenerate", partner
        # second signature: a null direction of the sensitivity Jacobian in which this
        # parameter and another both carry weight
        if null_vec is not None and np.isfinite(cond) and cond < DEGENERATE_COND                 and abs(float(null_vec[index])) > NULL_WEIGHT:
            weights = [(abs(float(null_vec[j])), param_names[j])
                       for j in range(len(param_names)) if j != index]
            if weights:
                w, partner = max(weights)
                if w > NULL_WEIGHT:
                    return "degenerate", partner
    if profile.lower_at_bound and profile.upper_at_bound:
        return "not_identifiable", None
    if open_side:
        return "bounded_only", None
    return "measurable", None


def fit_rate_curve(model: RateModel, gammadot, eta, *, level: float = 0.95,
                   seed: int = 0, n_starts: int = 8, temperature_C: float = float("nan"),
                   geometry: str = "unspecified") -> IdentifiedModel:
    """Fit ``model`` to a measured ``(gammadot, eta)`` curve and profile every parameter.

    The envelope is the measured rate range of the data passed in; ``eta_at`` on a
    record built from this result will refuse to extrapolate beyond it.

    ``n_starts`` IS NOT A COUNT OF DISTINCT OPTIMA. ``dst_fit.fit_model`` builds its start
    list as ``[start, heuristic, *(n_starts - 1) log-uniform draws]``, and three things
    make the effective count smaller than the nominal one: ``start`` (this function always
    passes ``model.guess``) can coincide with ``dst_fit``'s own heuristic; several random
    draws can fall into the same basin; and for a 5-parameter model, 8 draws from a box
    spanning decades is sparse. So a fit is not "well searched" because ``n_starts`` is 8.

    What that costs was MEASURED rather than left as a caveat
    (``scripts/multistart_seed_scatter.py`` -> ``data/horner2020/multistart_scatter.md``,
    2026-09-06): refitting all 20 blood curves under four models with five different seeds
    moves the OBJECTIVE by at most ~1e-9 relative, so no start is landing on a worse
    optimum -- the search is fine. Where the PARAMETERS move, they move on a flat
    direction, and no record with all parameters ``measurable`` moved by more than 1e-3
    (54 such records, worst 6.0e-5).

    **But "the flagged set catches every mover" is FALSE, and the pre-registered gate that
    said so failed.** Herschel-Bulkley on donor G moves 4.149 in ``tau_y`` while that
    parameter is graded ``bounded_only`` -- which is neither ``degenerate`` nor
    ``not_identifiable``, and so is NOT one of the verdicts
    ``solver.openfoam.check_exportable`` refuses. Do not read ``bounded_only`` as
    "determined, with a one-sided interval": it can be a parameter the optimiser places
    anywhere along a flat direction. The seed is recorded on the result, so such a value
    is reproducible; it is not, on its own, a measurement.
    """
    gd = np.asarray(gammadot, dtype=float)
    et = np.asarray(eta, dtype=float)
    if gd.shape != et.shape or gd.ndim != 1:
        raise ValueError(f"gammadot and eta must be equal-length 1-d, got {gd.shape} and {et.shape}")
    good = np.isfinite(gd) & np.isfinite(et) & (gd > 0) & (et > 0)
    n_dropped = int((~good).sum())
    if n_dropped:
        print(f"fit_rate_curve: dropping {n_dropped} of {gd.size} points (non-finite or <= 0)")
    gd, et = gd[good], et[good]
    if gd.size <= model.n_params:
        raise ValueError(f"{gd.size} points cannot determine {model.n_params} parameters")

    # dst_fit's fitter is generic in the independent variable; its own guess is keyed on
    # model name and is not this model's, so the guess is passed explicitly.
    fit = fit_model(model, gd, et, n_starts=n_starts, seed=seed, start=model.guess(gd, et))
    if not fit.success:
        raise RuntimeError(f"fit of {model.name} failed: {fit.message}")

    cond, null_vec = null_direction(model, gd, np.asarray(fit.params, dtype=float))
    params: dict[str, ParamEstimate] = {}
    for i, name in enumerate(model.param_names):
        prof = profile_likelihood(model, gd, et, name, level=level, fit=fit, seed=seed)
        verdict, partner = verdict_for(prof, fit.corr[i], i, model.param_names,
                                       cond=cond, null_vec=null_vec)
        params[name] = ParamEstimate(
            name=name, value=float(fit.params[i]), lo=float(prof.lower), hi=float(prof.upper),
            verdict=verdict, lo_at_bound=bool(prof.lower_at_bound),
            hi_at_bound=bool(prof.upper_at_bound), wald_stderr=float(prof.wald_stderr),
            partner=partner,
        )
    envelope = Envelope(gammadot_min=float(gd.min()), gammadot_max=float(gd.max()),
                        temperature_C=temperature_C, geometry=geometry)
    # residuals recomputed here rather than read off the fitter, so the ordering is
    # by rate whatever order the data arrived in
    order = np.argsort(gd)
    with np.errstate(all="ignore"):
        resid = np.log(et[order]) - np.log(model(gd[order], *np.asarray(fit.params, dtype=float)))
    structure = runs_test(resid)
    return IdentifiedModel(model=model, params=params, fit=fit, envelope=envelope, level=level,
                           residual_structure=structure)
