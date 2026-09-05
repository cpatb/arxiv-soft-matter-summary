"""One pint registry for the whole package.

Boundaries convert once, here; everything inside is an SI float. A registry record
stores ``value`` and ``unit`` and never ``si``; the SI magnitude is computed at load.
"""

from __future__ import annotations

import pint

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


def to_si(value: float, unit: str) -> float:
    """SI base magnitude of ``value`` in ``unit``. Raises on an unknown unit."""
    return float(Q_(value, unit).to_base_units().magnitude)


def dimensionality(unit: str) -> str:
    """Canonical dimensionality string, e.g. '[mass] / [length] / [time]' for Pa s."""
    return str(ureg.Unit(unit).dimensionality)
