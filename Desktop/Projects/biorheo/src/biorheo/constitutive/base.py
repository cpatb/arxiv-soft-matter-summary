"""Rate-space constitutive models as ``dst_fit.models.Model`` instances.

Why reuse ``dst_fit.Model`` rather than define a new interface
--------------------------------------------------------------
``dst_fit.fit_model`` and ``dst_fit.profile_likelihood`` only require
``model(x, *params)`` to return a viscosity for an independent variable ``x``;
neither cares whether ``x`` is a stress or a shear rate. ``dst_fit`` fits in stress
space because a discontinuous shear-thickening flow curve is multivalued in rate.
The generalised-Newtonian fluids here (yield stress, thinning) are single-valued in
rate, which is the variable a solver hands a viscosity model, so they are written
``eta(gammadot)`` and tagged ``control="rate"``. One fitter, one profiler, two
control variables, and every identifiability verdict in ``dst_fit`` applies here.

The one thing the ``dst_fit`` fitter keys on the model *name* is its initial guess,
with a bounds-midpoint fallback for unknown names. A midpoint of bounds spanning
decades is a poor start, so a ``RateModel`` carries its own ``guess`` and
``biorheo.identify`` passes it through ``fit_model(start=...)``.

Regularisation is not here. A constitutive definition returns ``inf`` below the
yield stress; a regularised variant is built by ``biorheo.constitutive.regularise``
with the scheme and its knob named in the model's ``fixed`` and ``name``, so a
solver export cannot receive a regularised model without saying which one.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Mapping

import numpy as np
from dst_fit.models import Model

GuessFn = Callable[[np.ndarray, np.ndarray], np.ndarray]


@dataclass(frozen=True, eq=False)
class RateModel(Model):
    """A ``dst_fit.Model`` whose independent variable is shear rate (1/s).

    ``limits`` documents the analytic limiting cases, each of which is exercised
    numerically in ``tests/test_gnf.py``; the SymPy derivations are in
    ``derivations/gnf_limits.md``. ``units`` gives the unit of each fitted
    parameter as a pint string, used by the registry to check a record's
    parameters carry the right dimension.
    """

    control: str = "rate"
    guess_fn: GuessFn | None = None
    limits: Mapping[str, str] = field(default_factory=dict)
    units: Mapping[str, str] = field(default_factory=dict)
    # Set only by ``regularise``: the scheme name and knob value, and the exact model
    # this one approximates, so the induced error can always be computed.
    regularisation: Mapping[str, float] | None = None
    base: "RateModel | None" = None

    def guess(self, gammadot: np.ndarray, eta: np.ndarray) -> np.ndarray:
        """A starting point inside the bounds, from the data's own scales."""
        lo, hi = np.array(self.bounds[0]), np.array(self.bounds[1])
        if self.guess_fn is None:
            g = np.sqrt(np.maximum(lo, 1e-12) * hi)  # log-midpoint
        else:
            g = np.asarray(self.guess_fn(np.asarray(gammadot, float), np.asarray(eta, float)), float)
        return np.clip(g, lo + 1e-12, hi - 1e-12)

    def stress(self, gammadot: np.ndarray, *params: float) -> np.ndarray:
        """sigma = eta * gammadot for gd > 0; NaN at zero rate (the 0 * inf there is
        model-specific: a yield-stress model tends to tau_y, a GNF to 0; read the
        model's ``limits`` instead of this function for that point)."""
        gd = np.asarray(gammadot, dtype=float)
        with np.errstate(invalid="ignore"):
            return np.where(gd > 0, self(gd, *params) * gd, np.nan)

    def params_dict(self, params) -> dict[str, float]:
        self.check_params(params)
        return dict(zip(self.param_names, (float(p) for p in params)))

    @property
    def is_regularised(self) -> bool:
        return self.regularisation is not None
