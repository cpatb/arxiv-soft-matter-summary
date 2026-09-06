"""The fluid record: what a registry entry is allowed to contain, and how it is used.

Rules (from ``03_ARCHITECTURE.md`` §4 data contracts):

- every parameter is ``{value, unit}``; ``si`` is computed at load, never stored;
- the envelope (rate range, temperature, geometry) is mandatory and enforced by
  ``eta_at``: evaluating outside the measured range raises rather than extrapolates;
- ``sources`` is a non-empty list of resolved DOIs;
- ``identifiability`` may only be written by ``biorheo.identify``, never by hand;
- ``licence.commercial_use`` is explicit, so a record cannot be shipped to a client
  without someone having read the licence.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Mapping

import numpy as np

from ..constitutive import RateModel
from ..units import to_si

__all__ = ["EnvelopeError", "ProvenanceError", "Parameter", "Envelope", "Source",
           "Licence", "FluidRecord", "round_out"]


def round_out(x: float, up: bool, sig: int = 6) -> float:
    """Round ``x`` to ``sig`` significant figures AWAY from zero-error: down for a lower
    bound (``up=False``), up for an upper bound (``up=True``).

    Why this exists. A record's envelope is written into YAML with ``%.6g``, and
    six-significant-figure rounding goes whichever way the digits happen to fall. When it
    rounds a MAXIMUM down or a MINIMUM up, the record's own extreme measured rate lands
    OUTSIDE its own envelope and ``FluidRecord.eta_at`` refuses to evaluate there.
    Measured 2026-09-06 before the fix: **12 of the 15 Carbopol records** excluded at
    least one of their own extreme rates (e.g. a raw maximum of 199.9612229703406 written
    as 199.961). The 60 blood records were unaffected, because 0.1 and 700 are exact at
    six figures. Rounding outward cannot exclude a measured rate, and it widens the
    envelope by at most one unit in the sixth figure.

    The generators call this; the loader cannot, because it never sees the raw curve.

    It works on the ROUNDED DECIMAL, not on pure arithmetic, because the decimal is what
    is written and read back. ``%.6g`` already rounds the safe way about half the time,
    and taking ``ceil(x/q)*q`` regardless widens those cases for nothing -- it turns an
    exact 0.1 into 0.100001, because the float 0.1 is a hair above the decimal 0.1. So:
    format, read back, and step one unit in the last figure ONLY if that landed on the
    wrong side. The value returned is exactly what ``f"{value:.6g}"`` will produce.
    """
    if not np.isfinite(x):
        raise ValueError(f"round_out needs a finite value, got {x!r}")
    if sig < 1:
        raise ValueError(f"sig must be at least 1, got {sig}")
    if x == 0.0:
        return 0.0
    v = float(f"{x:.{sig}g}")
    if (v >= x) if up else (v <= x):
        return v
    # The step is one unit in x's OWN last significant figure, not v's. They differ when
    # the rounding crossed a power of ten -- %.6g takes 999999.9 to 1e+06, whose quantum
    # is 10 where x's is 1, and stepping by v's would return 999990 instead of 999999,
    # ten times wider than necessary. Caught by the width property in test_registry.py.
    q = 10.0 ** (np.floor(np.log10(abs(x))) - (sig - 1))
    out = float(f"{v + q if up else v - q:.{sig}g}")
    if (out < x) if up else (out > x):  # cannot happen for finite x; refuse rather than lie
        raise AssertionError(f"round_out({x!r}, up={up}) produced {out!r}, which does not "
                             "bracket the input; the envelope would exclude a measured rate")
    return out


class ProvenanceError(ValueError):
    """A record cites nothing, cites something unresolved, or claims what only the
    identify layer may write."""


class EnvelopeError(ValueError):
    """The measured range is malformed, or a value was requested outside it."""


@dataclass(frozen=True)
class Parameter:
    name: str
    value: float
    unit: str
    si: float  # computed at load from (value, unit)


@dataclass(frozen=True)
class Envelope:
    gammadot_min: float  # 1/s, SI
    gammadot_max: float  # 1/s, SI
    temperature_C: float
    geometry: str

    def contains(self, gammadot: np.ndarray) -> np.ndarray:
        gd = np.asarray(gammadot, dtype=float)
        return (gd >= self.gammadot_min) & (gd <= self.gammadot_max)


@dataclass(frozen=True)
class Source:
    doi: str
    role: str  # "dataset" | "paper" | "fit" | ...
    resolved_on: str
    label: str


@dataclass(frozen=True)
class Licence:
    name: str
    commercial_use: bool


@dataclass(frozen=True)
class FluidRecord:
    id: str
    fluid: str
    model: RateModel
    parameters: Mapping[str, Parameter]
    envelope: Envelope
    sources: tuple[Source, ...]
    licence: Licence
    notes: str = ""
    identifiability: Mapping[str, str] | None = None  # written by biorheo.identify only
    path: str | None = None
    extra: Mapping[str, object] = field(default_factory=dict)

    @property
    def commercial_ok(self) -> bool:
        return self.licence.commercial_use

    @property
    def verdicts(self) -> dict[str, str] | None:
        """Per-parameter identifiability verdict, or None when nothing has been fitted."""
        if self.identifiability is None:
            return None
        return {k: v["verdict"] for k, v in self.identifiability["parameters"].items()}

    @property
    def params_si(self) -> np.ndarray:
        """Parameter values in SI, in the model's own order."""
        return np.array([self.parameters[p].si for p in self.model.param_names])

    def eta_at(self, gammadot) -> np.ndarray:
        """Viscosity (Pa s) at ``gammadot`` (1/s), refusing anything outside the envelope."""
        gd = np.asarray(gammadot, dtype=float)
        inside = self.envelope.contains(gd)
        if not np.all(inside):
            bad = gd[~inside]
            raise EnvelopeError(
                f"{self.id}: gammadot {bad.tolist()} outside measured range "
                f"[{self.envelope.gammadot_min:g}, {self.envelope.gammadot_max:g}] 1/s; "
                "the record does not extrapolate"
            )
        return self.model(gd, *self.params_si)

    def sigma_at(self, gammadot) -> np.ndarray:
        gd = np.asarray(gammadot, dtype=float)
        return self.eta_at(gd) * gd


def make_parameter(name: str, spec: object) -> Parameter:
    """Build a Parameter from the YAML ``{value, unit}`` mapping, computing ``si``."""
    if not isinstance(spec, Mapping):
        raise ValueError(f"parameter {name!r} must be a mapping {{value, unit}}, got {spec!r}")
    if "si" in spec:
        raise ValueError(
            f"parameter {name!r} stores an 'si' field; si is computed at load, never stored"
        )
    missing = {"value", "unit"} - set(spec)
    if missing:
        raise ValueError(f"parameter {name!r} is missing {sorted(missing)}: needs value and unit")
    value = float(spec["value"])
    unit = str(spec["unit"])
    if not np.isfinite(value):
        raise ValueError(f"parameter {name!r} is not finite: {value}")
    return Parameter(name=name, value=value, unit=unit, si=to_si(value, unit))
