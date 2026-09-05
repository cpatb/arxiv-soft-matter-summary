import pytest

from biorheo.units import dimensionality, to_si


def test_pascal_second_to_si_is_identity():
    assert to_si(0.056, "Pa*s") == pytest.approx(0.056)


def test_centipoise_converts():
    # 1 cP = 1e-3 Pa s exactly.
    assert to_si(3.45, "cP") == pytest.approx(3.45e-3)


def test_unknown_unit_raises():
    with pytest.raises(Exception):
        to_si(1.0, "furlongs_per_fortnight_squared")


def test_viscosity_dimensionality():
    assert dimensionality("Pa*s") == dimensionality("cP")
    assert dimensionality("Pa*s") != dimensionality("Pa")
