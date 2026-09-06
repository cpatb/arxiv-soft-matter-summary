import pytest

from biorheo.units import dimensionality, same_dimension, to_si


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


def test_same_dimension_is_not_a_string_comparison():
    """`str()` of a pint dimensionality rounds each exponent to 6 significant figures.
    Herschel-Bulkley's `k` has unit `Pa*s**n` with a FITTED n, so the loader's parameter
    check runs on exponents that differ in the far digits. The string comparison the
    loader used until 2026-09-05 ACCEPTS two such units as equal; the dimensional one
    does not. Recorded because the failure direction is a silent accept."""
    a, b = "Pa*s**0.3333333", "Pa*s**0.33333334"
    assert dimensionality(a) == dimensionality(b)      # the old test, still true
    assert not same_dimension(a, b)                    # what the loader now asks
    # and it still says yes to the things it must
    assert same_dimension("cP", "Pa*s") and same_dimension("N/m**2", "Pa")
    assert not same_dimension("Pa*s", "Pa")


def test_the_loader_rejects_a_unit_that_only_stringifies_the_same():
    """End to end: a record whose `k` unit differs from its own `n` below the 6th
    significant figure must be refused. It passed before 2026-09-05."""
    from biorheo.constitutive import get_model
    from biorheo.registry.loader import _expected_unit

    model = get_model("herschel_bulkley")
    si = {"tau_y": 1.0, "k": 1.0, "n": 0.3333333}
    expected = _expected_unit(model.units["k"], si)
    assert expected == "Pa*s**0.3333333"
    assert dimensionality(expected) == dimensionality("Pa*s**0.33333334")
    assert not same_dimension(expected, "Pa*s**0.33333334")
