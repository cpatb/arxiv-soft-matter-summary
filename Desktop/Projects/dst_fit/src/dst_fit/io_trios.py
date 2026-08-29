"""Reader for TA Instruments TRIOS flow-sweep exports (legacy .xls).

Why this module exists
----------------------
The TRIOS export writes ``Stress`` and ``Normal stress`` in **MPa**, which for a
dense suspension means the physically interesting numbers live around 1e-8 to
1e-4 in the stored column. Nothing is lost (the values carry six significant
figures) but anything that prints, rounds, or thresholds the raw column without
converting first will silently see zeros. Reading this file with a default
``head()`` shows a column of ``0.0``; that is a display artifact, and mistaking
it for the data is the trap this reader exists to close.

Everything is converted to SI on the way in, and stress is *also* recomputed
from the two independent columns,

    sigma = eta * gamma_dot

which is an identity for steady simple shear, not a model. Comparing the two is
a real check on the unit interpretation: if ``MPa`` were misread the two
disagree by six orders of magnitude. On the corn-starch series they agree to
about 5e-6 relative, which is the storage precision. See
:func:`check_stress_consistency`.

Sheet layout, as exported (0-indexed rows of the raw sheet):

    row 0   sheet title in col 0, NaN elsewhere      e.g. "Flow sweep - 1"
    row 1   column names                             Stress, Shear rate, ...
    row 2   units                                    MPa, 1/s, Pa.s, min, degC, MPa
    row 3+  data

Everything returned by :func:`read_flow_sweep` is in SI (Pa, 1/s, Pa.s, s, degC).
"""

from __future__ import annotations

import math
import re
import warnings
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# Unit strings TRIOS emits, mapped to the multiplier that takes them to SI.
# Unknown units raise rather than being silently assumed to be SI.
_UNIT_TO_SI: dict[str, float] = {
    "Pa": 1.0,
    "kPa": 1e3,
    "MPa": 1e6,
    "1/s": 1.0,
    "s^-1": 1.0,
    "Pa.s": 1.0,
    "Pa-s": 1.0,
    "mPa.s": 1e-3,
    "s": 1.0,
    "min": 60.0,
    "degC": 1.0,
    "C": 1.0,
}

# TRIOS column name -> our field name. Matched case-insensitively after
# collapsing whitespace.
_COLUMN_ALIASES: dict[str, str] = {
    "stress": "stress",
    "shear stress": "stress",
    "shear rate": "shear_rate",
    "viscosity": "viscosity",
    "step time": "step_time",
    "temperature": "temperature",
    "normal stress": "normal_stress",
}

_REQUIRED = ("shear_rate", "viscosity")


@dataclass(frozen=True)
class FlowSweep:
    """One steady flow sweep, in SI units.

    ``stress`` is the reconstructed ``viscosity * shear_rate``, never the
    exported stress column. ``stress_reported`` keeps the exported column so a
    caller can audit the unit interpretation; it is quantised to the export's
    display precision and is usually mostly zeros.
    """

    shear_rate: np.ndarray  # 1/s
    viscosity: np.ndarray  # Pa.s
    stress: np.ndarray  # Pa, reconstructed
    stress_reported: np.ndarray | None  # Pa, as exported (low precision)
    temperature: np.ndarray | None  # degC
    step_time: np.ndarray | None  # s
    normal_stress: np.ndarray | None  # Pa
    source: Path
    sheet: str
    metadata: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return int(self.shear_rate.size)

    @property
    def phi(self) -> float | None:
        """Volume fraction parsed from the filename, as a fraction (not %)."""
        return self.metadata.get("phi")

    def summary(self) -> str:
        eta = self.viscosity
        i_min = int(np.argmin(eta))
        thickening = float(eta[-1] / eta[i_min])
        return (
            f"{self.source.name:>10s}  n={len(self):>3d}  "
            f"gdot {self.shear_rate.min():9.4f} to {self.shear_rate.max():9.2f} 1/s  "
            f"eta {eta.min():.5f} to {eta.max():.5f} Pa.s  "
            f"eta_max/eta_min(after min) {thickening:6.2f}x"
        )


def _clean(text: object) -> str:
    """Normalise a header cell: strip non-ascii (TRIOS writes a raw degree sign)."""
    s = "" if text is None or (isinstance(text, float) and math.isnan(text)) else str(text)
    s = s.encode("ascii", "ignore").decode("ascii")
    return re.sub(r"\s+", " ", s).strip()


def _find_header_row(raw: pd.DataFrame, max_scan: int = 10) -> int:
    """Return the index of the row holding the column names.

    Located by content, not by a hardcoded offset, so an export with an extra
    banner row still reads correctly.
    """
    for i in range(min(max_scan, len(raw))):
        cells = {_clean(c).lower() for c in raw.iloc[i]}
        if "shear rate" in cells and "viscosity" in cells:
            return i
    raise ValueError(
        "no header row containing both 'Shear rate' and 'Viscosity' in the "
        f"first {max_scan} rows; is this a TRIOS flow-sweep export?"
    )


def _to_si(values: np.ndarray, unit: str, column: str) -> np.ndarray:
    key = _clean(unit)
    if key not in _UNIT_TO_SI:
        raise ValueError(
            f"column {column!r} has unrecognised unit {key!r}. Refusing to guess; "
            f"add it to _UNIT_TO_SI with an explicit factor. Known: "
            f"{sorted(_UNIT_TO_SI)}"
        )
    return values * _UNIT_TO_SI[key]


def _phi_from_name(path: Path) -> float | None:
    """Parse a leading integer percent from names like ``47_e1.xls`` -> 0.47."""
    m = re.match(r"^(\d{1,3})(?:[_\-.]|$)", path.stem)
    if not m:
        return None
    pct = int(m.group(1))
    return pct / 100.0 if 0 < pct <= 100 else None


def read_flow_sweep(path: str | Path, sheet: str | None = None) -> FlowSweep:
    """Read one TRIOS flow sweep into SI units.

    Parameters
    ----------
    path
        The exported ``.xls`` file.
    sheet
        Sheet name. Defaults to the first sheet whose name starts with
        "Flow sweep".

    Raises
    ------
    ValueError
        If the header cannot be located, a unit is unrecognised, required
        columns are missing, or the data violate a physical invariant
        (non-positive shear rate or viscosity, NaNs in required columns).
    """
    path = Path(path)
    with warnings.catch_warnings():
        # TRIOS writes a legacy .xls whose styling xlrd grumbles about.
        warnings.simplefilter("ignore")
        book = pd.read_excel(path, sheet_name=None, header=None)

    if sheet is None:
        candidates = [k for k in book if _clean(k).lower().startswith("flow sweep")]
        if not candidates:
            raise ValueError(
                f"{path.name}: no sheet starting with 'Flow sweep'; found {list(book)}"
            )
        sheet = candidates[0]
    raw = book[sheet]

    h = _find_header_row(raw)
    names = [_clean(c).lower() for c in raw.iloc[h]]
    units = [_clean(c) for c in raw.iloc[h + 1]]
    body = raw.iloc[h + 2 :].reset_index(drop=True)

    columns: dict[str, np.ndarray] = {}
    for j, name in enumerate(names):
        field_name = _COLUMN_ALIASES.get(name)
        if field_name is None:
            continue
        values = np.asarray(pd.to_numeric(body.iloc[:, j], errors="coerce"), dtype=float)
        columns[field_name] = _to_si(values, units[j], field_name)

    missing = [c for c in _REQUIRED if c not in columns]
    if missing:
        raise ValueError(
            f"{path.name}[{sheet}]: missing required column(s) {missing}; "
            f"parsed {sorted(columns)} from header {names}"
        )

    # Drop trailing all-NaN rows that some exports leave behind, then validate.
    keep = np.isfinite(columns["shear_rate"]) & np.isfinite(columns["viscosity"])
    if not keep.any():
        raise ValueError(f"{path.name}[{sheet}]: no finite (shear_rate, viscosity) rows")
    for k in list(columns):
        columns[k] = columns[k][keep]

    gdot = columns["shear_rate"]
    eta = columns["viscosity"]
    if not np.all(gdot > 0):
        raise ValueError(f"{path.name}[{sheet}]: non-positive shear rate present")
    if not np.all(eta > 0):
        raise ValueError(f"{path.name}[{sheet}]: non-positive viscosity present")

    # The identity TRIOS itself used to build the viscosity column.
    stress = eta * gdot

    metadata = {"phi": _phi_from_name(path), "units_row": dict(zip(names, units))}
    metadata.update(_read_details(book))

    return FlowSweep(
        shear_rate=gdot,
        viscosity=eta,
        stress=stress,
        stress_reported=columns.get("stress"),
        temperature=columns.get("temperature"),
        step_time=columns.get("step_time"),
        normal_stress=columns.get("normal_stress"),
        source=path,
        sheet=str(sheet),
        metadata=metadata,
    )


def _read_details(book: dict) -> dict:
    """Pull instrument/operator/rundate out of the 'Details' sheet if present."""
    for key in book:
        if _clean(key).lower() != "details":
            continue
        df = book[key]
        if df.shape[1] < 2:
            return {}
        wanted = {"instrument name", "operator", "rundate", "proceduresegments"}
        out = {}
        for _, row in df.iterrows():
            label = _clean(row.iloc[0]).lower()
            if label in wanted:
                out[label.replace(" ", "_")] = _clean(row.iloc[1])
        return out
    return {}


def check_stress_consistency(sweep: FlowSweep, rtol: float = 1e-2) -> dict:
    """Cross-check ``eta * gamma_dot`` against the exported stress column.

    This is a test of the unit interpretation, not of the physics. If ``MPa``
    were misread as ``Pa`` the two would differ by a factor of 1e6, so agreement
    at any sane tolerance pins the conversion. The storage precision is six
    significant figures, so the achievable agreement is around 5e-6 relative and
    the default 1% tolerance leaves three orders of margin.

    Points with a non-positive reported stress are excluded rather than counted
    as agreement; the count is reported so a caller can see how much of the
    curve the check actually covered.
    """
    reported = sweep.stress_reported
    n = len(sweep)
    if reported is None:
        return {"comparable": 0, "n": n, "ok": False, "reason": "no stress column"}

    informative = reported > 0
    n_excluded = int((~informative).sum())

    if not informative.any():
        return {
            "comparable": 0,
            "n": n,
            "n_excluded": n_excluded,
            "ok": False,
            "reason": "exported stress is non-positive at every point",
        }

    rel = np.abs(sweep.stress[informative] - reported[informative]) / reported[informative]
    worst = float(rel.max())
    return {
        "comparable": int(informative.sum()),
        "n": n,
        "n_excluded": n_excluded,
        "max_rel_resid": worst,
        "rtol": rtol,
        "ok": bool(worst <= rtol),
        "reason": "",
    }


def read_series(paths) -> list[FlowSweep]:
    """Read several sweeps and return them sorted by volume fraction."""
    sweeps = [read_flow_sweep(p) for p in paths]
    return sorted(sweeps, key=lambda s: (s.phi is None, s.phi))
