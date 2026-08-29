"""Tests for the TRIOS reader.

The synthetic fixtures deliberately reproduce the awkward parts of a real
export rather than a clean table: a sheet-title row above the header, a separate
units row, stress in MPa, and a non-ascii degree sign in the temperature unit. A
fixture that omits those would pass while the reader failed on real files.

One real-data test is included and skips if the corn-starch directory is absent,
so the suite still runs on a machine without it.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from dst_fit import check_stress_consistency, read_flow_sweep, read_series
from dst_fit.io_trios import _find_header_row, _phi_from_name

REAL_DATA_DIR = Path(r"C:\Users\compl\Desktop\Projects\Corn Starch Bulk Rheology AM")


# --------------------------------------------------------------------------
# fixture builder
# --------------------------------------------------------------------------
def write_sweep(
    path: Path,
    gdot: np.ndarray,
    eta: np.ndarray,
    *,
    stress_unit: str = "MPa",
    rate_unit: str = "1/s",
    banner_rows: int = 1,
    sheet: str = "Flow sweep - 1",
    drop_column: str | None = None,
) -> Path:
    """Write an .xlsx laid out like a TRIOS export."""
    stress_si = eta * gdot
    scale = {"MPa": 1e-6, "kPa": 1e-3, "Pa": 1.0}[stress_unit]

    names = ["Stress", "Shear rate", "Viscosity", "Step time", "Temperature"]
    units = [stress_unit, rate_unit, "Pa.s", "min", "\u00b0C"]
    cols = [
        stress_si * scale,
        gdot,
        eta,
        np.arange(1, gdot.size + 1) * 0.5,
        np.full(gdot.size, 25.0),
    ]
    if drop_column is not None:
        j = names.index(drop_column)
        names.pop(j)
        units.pop(j)
        cols.pop(j)

    rows: list[list] = [[sheet] + [None] * (len(names) - 1) for _ in range(banner_rows)]
    rows.append(list(names))
    rows.append(list(units))
    rows.extend([list(v) for v in np.array(cols, dtype=object).T])

    with pd.ExcelWriter(path, engine="openpyxl") as xl:
        pd.DataFrame(rows).to_excel(xl, sheet_name=sheet, header=False, index=False)
    return path


@pytest.fixture
def gdot() -> np.ndarray:
    return np.logspace(-1, 3, 25)


@pytest.fixture
def eta(gdot: np.ndarray) -> np.ndarray:
    # thinning then thickening, so the curve has the shape the reader will meet
    return 0.1 * gdot**-0.5 + 1e-4 * gdot


# --------------------------------------------------------------------------
# happy path
# --------------------------------------------------------------------------
def test_reads_si_units(tmp_path, gdot, eta):
    p = write_sweep(tmp_path / "47_e1.xlsx", gdot, eta)
    s = read_flow_sweep(p)

    assert len(s) == gdot.size
    np.testing.assert_allclose(s.shear_rate, gdot, rtol=1e-9)
    np.testing.assert_allclose(s.viscosity, eta, rtol=1e-9)
    # stress is reconstructed, and lands in Pa not MPa
    np.testing.assert_allclose(s.stress, eta * gdot, rtol=1e-9)
    # step time converted min -> s
    np.testing.assert_allclose(s.step_time[0], 30.0, rtol=1e-9)


def test_mpa_is_converted_not_assumed(tmp_path, gdot, eta):
    """The same numbers labelled Pa and MPa must differ by exactly 1e6."""
    a = read_flow_sweep(write_sweep(tmp_path / "a.xlsx", gdot, eta, stress_unit="MPa"))
    b = read_flow_sweep(write_sweep(tmp_path / "b.xlsx", gdot, eta, stress_unit="Pa"))
    np.testing.assert_allclose(a.stress_reported, b.stress_reported, rtol=1e-9)


def test_header_found_despite_extra_banner(tmp_path, gdot, eta):
    p = write_sweep(tmp_path / "c.xlsx", gdot, eta, banner_rows=3)
    s = read_flow_sweep(p)
    assert len(s) == gdot.size


def test_phi_parsed_from_filename():
    assert _phi_from_name(Path("47_e1.xls")) == 0.47
    assert _phi_from_name(Path("40.xls")) == 0.40
    assert _phi_from_name(Path("sample_e1.xls")) is None
    assert _phi_from_name(Path("999_e1.xls")) is None


def test_series_sorted_by_phi(tmp_path, gdot, eta):
    for name in ("52_e1.xlsx", "40_e1.xlsx", "47_e1.xlsx"):
        write_sweep(tmp_path / name, gdot, eta)
    sweeps = read_series(sorted(tmp_path.glob("*.xlsx")))
    assert [s.phi for s in sweeps] == [0.40, 0.47, 0.52]


# --------------------------------------------------------------------------
# failure modes: the reader must refuse, not guess
# --------------------------------------------------------------------------
def test_unknown_unit_raises(tmp_path, gdot, eta):
    p = write_sweep(tmp_path / "d.xlsx", gdot, eta, rate_unit="rpm")
    with pytest.raises(ValueError, match="unrecognised unit"):
        read_flow_sweep(p)


def test_missing_required_column_raises(tmp_path, gdot, eta):
    p = write_sweep(tmp_path / "e.xlsx", gdot, eta, drop_column="Viscosity")
    with pytest.raises(ValueError, match="no header row"):
        read_flow_sweep(p)


def test_non_positive_viscosity_raises(tmp_path, gdot, eta):
    bad = eta.copy()
    bad[3] = 0.0
    p = write_sweep(tmp_path / "f.xlsx", gdot, bad)
    with pytest.raises(ValueError, match="non-positive viscosity"):
        read_flow_sweep(p)


def test_non_positive_shear_rate_raises(tmp_path, gdot, eta):
    bad = gdot.copy()
    bad[0] = 0.0
    p = write_sweep(tmp_path / "g.xlsx", bad, eta)
    with pytest.raises(ValueError, match="non-positive shear rate"):
        read_flow_sweep(p)


def test_no_flow_sweep_sheet_raises(tmp_path, gdot, eta):
    p = write_sweep(tmp_path / "h.xlsx", gdot, eta, sheet="Frequency sweep - 1")
    with pytest.raises(ValueError, match="no sheet starting with"):
        read_flow_sweep(p)


def test_find_header_row_reports_clearly():
    junk = pd.DataFrame([["nothing", "useful"], ["here", "either"]])
    with pytest.raises(ValueError, match="no header row"):
        _find_header_row(junk)


# --------------------------------------------------------------------------
# consistency check
# --------------------------------------------------------------------------
def test_consistency_passes_on_clean_data(tmp_path, gdot, eta):
    s = read_flow_sweep(write_sweep(tmp_path / "i.xlsx", gdot, eta))
    chk = check_stress_consistency(s)
    assert chk["ok"]
    assert chk["comparable"] == gdot.size
    assert chk["max_rel_resid"] < 1e-9


def test_consistency_catches_a_wrong_unit(tmp_path, gdot, eta):
    """If the exported column were really kPa and we read it as MPa, the check
    must fail rather than quietly accept a factor of 1000."""
    s = read_flow_sweep(write_sweep(tmp_path / "j.xlsx", gdot, eta))
    corrupted = s.__class__(
        **{
            **{k: getattr(s, k) for k in s.__dataclass_fields__},
            "stress_reported": s.stress_reported * 1e3,
        }
    )
    assert not check_stress_consistency(corrupted)["ok"]


def test_zero_stress_points_excluded_not_counted_as_agreement(tmp_path, gdot, eta):
    s = read_flow_sweep(write_sweep(tmp_path / "k.xlsx", gdot, eta))
    zeroed = s.stress_reported.copy()
    zeroed[:5] = 0.0
    corrupted = s.__class__(
        **{
            **{k: getattr(s, k) for k in s.__dataclass_fields__},
            "stress_reported": zeroed,
        }
    )
    chk = check_stress_consistency(corrupted)
    assert chk["n_excluded"] == 5
    assert chk["comparable"] == gdot.size - 5


# --------------------------------------------------------------------------
# real data
# --------------------------------------------------------------------------
@pytest.mark.skipif(not REAL_DATA_DIR.is_dir(), reason="corn-starch data not present")
def test_real_corn_starch_series():
    sweeps = read_series(sorted(REAL_DATA_DIR.glob("*_e1.xls")))
    assert [s.phi for s in sweeps] == [0.40, 0.45, 0.47, 0.50, 0.52, 0.54]
    for s in sweeps:
        chk = check_stress_consistency(s)
        assert chk["ok"], f"{s.source.name}: {chk}"
        assert np.all(s.shear_rate > 0) and np.all(s.viscosity > 0)
        # stress lands in Pa: a corn-starch sweep spans roughly 1e-2 to 1e3 Pa
        assert 1e-3 < s.stress.max() < 1e5, s.stress.max()
