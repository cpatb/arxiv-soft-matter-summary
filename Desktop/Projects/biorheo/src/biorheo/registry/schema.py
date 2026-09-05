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
           "Licence", "FluidRecord"]


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
