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
    """Canonical dimensionality string, e.g. '[mass] / [length] / [time]' for Pa s.

    For DISPLAY and for comparison against a literal like ``"[time]"``. To ask whether
    two units share a dimension, use ``same_dimension``, not this: see its docstring."""
    return str(ureg.Unit(unit).dimensionality)


def same_dimension(a: str, b: str) -> bool:
    """Do ``a`` and ``b`` have the same physical dimension?

    Compares pint's ``UnitsContainer`` objects, not their ``str()``, and the difference is
    a SILENT ACCEPT rather than a cosmetic one. ``str()`` of a dimensionality rounds each
    exponent to 6 significant figures, and this package's models carry FITTED exponents
    (``Pa*s**n`` for Herschel-Bulkley, substituted with the record's own ``n``), so the
    exponent is whatever float the fit produced. Measured 2026-09-05::

        dimensionality("Pa*s**0.3333333")  -> '[mass] / [length] / [time] ** 1.66667'
        dimensionality("Pa*s**0.33333334") -> '[mass] / [length] / [time] ** 1.66667'

    -- equal as strings, different as dimensions. The string test is therefore strictly
    more permissive than the dimensional one and would admit a record whose unit
    disagrees with its model below the 7th digit. Added 2026-09-05, the referee item
    left open in NEXT_STEPS item 5.
    """
    return ureg.Unit(a).dimensionality == ureg.Unit(b).dimensionality
