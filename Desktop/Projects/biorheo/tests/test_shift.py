"""Item 21: the temperature-shift layer.

Two layers, as everywhere else here. The first tests `biorheo.constitutive.shift` against
the identity it exists to implement, at parameter values chosen in this file rather than
read off a record -- the identity is algebraic and a check on one fit is not a check on it.
The second pins the gate verdicts in `data/temperature_shift.json`, INCLUDING the one that
failed.

The load-bearing test here is not that the map is right for the six families. It is that
the map is not the IDENTITY (`test_the_map_actually_moves_the_curve`) and that the
collapse residual can tell a real shift from a shape mismatch. Everything else is passed
by a routine that returns its input.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from biorheo.apps.injectability import Needle, Syringe, injection_force
from biorheo.constitutive.gnf import MODELS, NEWTONIAN
from biorheo.constitutive.shift import (
    GAS_CONSTANT, SHIFTABLE, ThermalShiftNotRepresentable, arrhenius_shift, collapse_shift,
    shift_envelope, shift_from_two_viscosities, shift_parameters, shift_record,
)
from biorheo.constitutive.structural import STRUCTURAL_STEADY, STRUCTURAL_STEADY_SEPARATE
from biorheo.registry import load_record
from biorheo.registry.schema import EnvelopeError

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "data" / "temperature_shift.json"
MD = REPO / "data" / "temperature_shift.md"
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"

PROBE = {
    "newtonian": np.array([3.5e-3]),
    "power_law": np.array([2.0, 0.45]),
    "herschel_bulkley": np.array([23.75, 8.9, 0.4526]),
    "casson": np.array([2.18e-3, 3.85e-3]),
    "cross": np.array([0.03, 2.4e-3, 8.0, 0.7]),
    "carreau_yasuda": np.array([0.03, 2.4e-3, 8.0, 2.0, 0.586]),
}
GD = np.logspace(-2, 3, 200)
FACTORS = (0.25, 0.5, 2.0, 4.0)


@pytest.fixture(scope="module")
def art() -> dict:
    assert JSON.exists(), "run scripts/temperature_shift.py"
    return json.loads(JSON.read_text(encoding="utf-8"))


# --- the identity the module exists to implement --------------------------------------

@pytest.mark.parametrize("name", sorted(PROBE))
@pytest.mark.parametrize("a", FACTORS)
def test_the_map_reproduces_the_shifted_curve(name, a):
    m, p = MODELS[name], PROBE[name]
    lhs = a * m(a * GD, *p)
    rhs = m(GD, *shift_parameters(m, p, a))
    assert np.max(np.abs(lhs / rhs - 1.0)) < 1e-12


@pytest.mark.parametrize("name", sorted(PROBE))
def test_the_map_actually_moves_the_curve(name):
    """**The discriminating test.** A `shift_parameters` that returned its argument would
    pass every identity test in this file, because `a*eta(a*gd) == eta(gd)` is exactly what
    it would be compared against -- no. It would fail THIS one, which compares the shifted
    curve against the UNSHIFTED one at the same rates."""
    m, p = MODELS[name], PROBE[name]
    mid = np.array([float(np.sqrt(GD[0] * GD[-1]))])
    base = float(m(mid, *p)[0])
    moved = float(m(mid, *shift_parameters(m, p, 2.0))[0])
    assert abs(moved / base - 1.0) > 0.10, (name, base, moved)


@pytest.mark.parametrize("name", sorted(PROBE))
def test_a_equals_one_is_bitwise_identity(name):
    """Bitwise, not approximate: an identity held only to a tolerance would hide a defect
    in the map behind that tolerance."""
    p = PROBE[name]
    assert np.array_equal(shift_parameters(MODELS[name], p, 1.0), p)


@pytest.mark.parametrize("name", sorted(PROBE))
def test_the_shift_composes(name):
    m, p = MODELS[name], PROBE[name]
    for a in FACTORS:
        for b in FACTORS:
            two = shift_parameters(m, shift_parameters(m, p, a), b)
            one = shift_parameters(m, p, a * b)
            assert np.allclose(two, one, rtol=1e-12), (name, a, b)


def test_the_yield_stress_is_invariant_and_the_viscous_parameter_is_not():
    """Derived in the pre-registration BEFORE this code existed, not imposed. The second
    half is the discriminating one: a map that changed nothing would report invariance by
    doing nothing at all."""
    for name, visc in (("casson", "eta_c"), ("herschel_bulkley", "k")):
        m, p = MODELS[name], PROBE[name]
        i_t, i_v = m.param_names.index("tau_y"), m.param_names.index(visc)
        for a in FACTORS:
            q = shift_parameters(m, p, a)
            assert q[i_t] == p[i_t], (name, a, "tau_y moved")
            assert q[i_v] != p[i_v], (name, a, f"{visc} did not move")


def test_the_exponent_on_k_is_n_not_one():
    """`k -> a^n k`, not `a k`. The two agree at `n = 1` and nowhere else, so this pins
    the exponent rather than the scaling."""
    m = MODELS["power_law"]
    k, n = 2.0, 0.45
    q = shift_parameters(m, np.array([k, n]), 4.0)
    assert q[0] == pytest.approx(k * 4.0 ** n)
    assert q[0] != pytest.approx(k * 4.0)
    assert q[1] == n


def test_negative_or_zero_shift_is_refused():
    for a in (0.0, -1.0):
        with pytest.raises(ValueError, match="must be positive"):
            shift_parameters(NEWTONIAN, np.array([1e-3]), a)


def test_structural_families_are_refused_for_the_right_reason():
    """Item 11's mistake was refusing something OpenFOAM could actually do, on a reason
    that was false. The reason here must name the STRUCTURE EQUATION and must not claim no
    counterpart exists -- the steady curve could be shifted; the relaxation time cannot."""
    for model in (STRUCTURAL_STEADY, STRUCTURAL_STEADY_SEPARATE):
        with pytest.raises(ThermalShiftNotRepresentable) as e:
            shift_parameters(model, np.ones(len(model.param_names)), 2.0)
        msg = str(e.value)
        assert "STRUCTURE EQUATION" in msg
        assert "no counterpart" not in msg.lower()
    assert set(SHIFTABLE) == set(PROBE)


# --- the shift factor is never invented ------------------------------------------------

def test_arrhenius_uses_the_scipy_gas_constant_and_needs_a_caller_activation_energy():
    from scipy.constants import R
    assert GAS_CONSTANT == R
    a = arrhenius_shift(2.0e4, 278.15, 298.15)
    expected = np.exp(2.0e4 / R * (1 / 278.15 - 1 / 298.15))
    assert a == pytest.approx(expected, rel=1e-15)
    assert a > 1.0, "colder must be more viscous for a positive activation energy"


def test_arrhenius_is_the_identity_at_the_reference_temperature():
    assert arrhenius_shift(5e4, 298.15, 298.15) == pytest.approx(1.0, abs=1e-15)


def test_arrhenius_refuses_non_kelvin_temperatures():
    """A caller passing celsius would silently get a wrong answer for anything above 0 C;
    below it the sign is impossible and this catches at least that."""
    with pytest.raises(ValueError, match="kelvin"):
        arrhenius_shift(2e4, -5.0, 25.0)


def test_shift_from_two_viscosities_is_the_ratio():
    assert shift_from_two_viscosities(43e-3, 24e-3) == pytest.approx(43.0 / 24.0)
    with pytest.raises(ValueError, match="positive"):
        shift_from_two_viscosities(0.0, 1.0)


def test_no_default_activation_energy_is_shipped_anywhere():
    """A plausible-looking default would be exactly the invented constant this repo
    forbids, and it would be invisible at the call site."""
    import inspect

    import biorheo.constitutive.shift as mod
    sig = inspect.signature(mod.arrhenius_shift)
    assert sig.parameters["E_a"].default is inspect.Parameter.empty


# --- the envelope, which is the trap ---------------------------------------------------

def test_the_envelope_is_divided_by_a():
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    for a in FACTORS:
        s = shift_record(rec, a)
        assert s.envelope.gammadot_min == pytest.approx(rec.envelope.gammadot_min / a, rel=1e-15)
        assert s.envelope.gammadot_max == pytest.approx(rec.envelope.gammadot_max / a, rel=1e-15)


def test_the_shifted_record_refuses_outside_its_own_shifted_envelope():
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    s = shift_record(rec, 4.0)
    with pytest.raises(EnvelopeError):
        s.eta_at(np.array([s.envelope.gammadot_max * 1.001]))
    # and it must still evaluate INSIDE, or the refusal is just a broken record
    inside = float(np.sqrt(s.envelope.gammadot_min * s.envelope.gammadot_max))
    assert np.isfinite(s.eta_at(np.array([inside]))[0])


def test_the_shifted_record_evaluates_to_the_shifted_curve():
    """End to end through the record layer, not just the parameter map: unit conversion,
    `params_si` and the envelope all sit between them."""
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    a = 2.0
    s = shift_record(rec, a)
    gd = np.logspace(np.log10(s.envelope.gammadot_min * 1.01),
                     np.log10(s.envelope.gammadot_max * 0.99), 50)
    assert np.allclose(s.eta_at(gd), a * rec.eta_at(a * gd), rtol=1e-12)


def test_the_identifiability_sidecar_is_dropped_not_carried():
    """An interval belongs to a fit against measured data. There is no data at the shifted
    temperature, so carrying the old interval across would be a confidence statement about
    measurements nobody made."""
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    assert rec.identifiability is not None, "the source record must HAVE one, or this is vacuous"
    s = shift_record(rec, 2.0)
    assert s.identifiability is None
    assert s.verdicts is None


def test_the_shifted_record_reports_as_a_point_estimate():
    """The consequence of dropping the sidecar, checked where a client would see it."""
    from biorheo.apps.report import DeviceSpec, build, render
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    spec = DeviceSpec(needle_radius=1.05e-4, needle_length=12.7e-3, barrel_radius=2.3e-3,
                      flow_rate=1e-6 / 60, label="t")
    d = build(shift_record(rec, 2.0), spec)
    assert "no_sidecar" in d.flags
    assert "POINT ESTIMATE" in render(d)


def test_the_shifted_record_says_so_in_its_notes_and_its_id():
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    s = shift_record(rec, 2.0)
    assert s.id.startswith(rec.id) and "shift" in s.id
    assert "TIME-TEMPERATURE SHIFTED" in s.notes
    assert "No data was measured at this temperature" in s.notes
    assert "; time-temperature shifted by a = 2" in s.envelope.geometry


def test_shift_envelope_can_record_the_new_temperature():
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    e = shift_envelope(rec.envelope, 2.0, temperature_C=5.0)
    assert e.temperature_C == 5.0
    assert shift_envelope(rec.envelope, 2.0).temperature_C == rec.envelope.temperature_C


# --- the collapse residual --------------------------------------------------------------

def test_collapse_recovers_a_shift_it_was_given():
    m, p = MODELS["cross"], PROBE["cross"]
    out = collapse_shift(GD, m(GD, *p), m(GD, *shift_parameters(m, p, 2.5)))
    assert out["a"] == pytest.approx(2.5, rel=1e-4)
    assert out["rms_log_residual"] < 1e-4
    assert out["n_used"] < out["n_total"], "shifting the axis must drop rates off the end"


def test_collapse_reports_a_large_residual_for_a_shape_mismatch():
    """**The half that makes the other half mean anything.** A best-fitting `a` always
    exists; without this, the routine would report a confident factor for two curves that
    do not collapse at all."""
    m, p = MODELS["cross"], PROBE["cross"]
    m2, p2 = MODELS["power_law"], PROBE["power_law"].copy()
    mid = np.array([float(np.sqrt(GD[0] * GD[-1]))])
    p2[0] *= float(m(mid, *p)[0] / m2(mid, *p2)[0])     # match magnitude, differ in shape
    out = collapse_shift(GD, m(GD, *p), m2(GD, *p2))
    assert out["rms_log_residual"] > 1e-3
    assert "meaningless rather than approximate" in out["note"]


def test_collapse_refuses_input_it_cannot_fit():
    with pytest.raises(ValueError, match="at least 4"):
        collapse_shift([1.0, 2.0], [1.0, 1.0], [2.0, 2.0])
    with pytest.raises(ValueError, match="shapes must match"):
        collapse_shift(GD, GD, GD[:-1])
    with pytest.raises(ValueError, match="positive"):
        collapse_shift([1.0, 2.0, 3.0, 4.0], [1.0, 1.0, 1.0, 0.0], [1.0, 1.0, 1.0, 1.0])


# --- end to end -------------------------------------------------------------------------

@pytest.mark.parametrize("a", FACTORS)
def test_the_newtonian_force_scales_exactly_as_a(a):
    """A closed form the force solver does not use: `dp = 8 mu L Q/(pi R^4)` is linear in
    `mu`, and the shift multiplies `mu` by `a`."""
    needle, syringe = Needle(radius=1.05e-4, length=12.7e-3), Syringe(barrel_radius=2.3e-3)
    Q, mu = 1e-6 / 60, 3.5e-3
    F0 = injection_force(Q, needle, syringe, NEWTONIAN, np.array([mu])).force
    F = injection_force(Q, needle, syringe, NEWTONIAN,
                        shift_parameters(NEWTONIAN, np.array([mu]), a)).force
    assert F == pytest.approx(a * F0, rel=1e-10)


# --- the shipped artifact -----------------------------------------------------------------

def test_every_gate_except_the_recorded_failure_passes(art):
    for name, g in art["gates"].items():
        if name == "G21_9_collapse_residual":
            continue
        assert g["pass"], (name, g)


def test_g219_is_recorded_failed_with_its_number_and_the_bar_is_not_moved(art):
    """A mis-specified bar is kept with its measurement, never edited out. The failing
    number must still be on the page."""
    g = art["gates"]["G21_9_collapse_residual"]
    assert g["pass"] is False
    assert g["clean"]["rms_log_residual"] > 1e-10
    assert "1e-10" in g["bar"] or "1e-10" in str(g["bar"])
    assert "interpolation error" in g["why_it_failed"]
    assert g["discriminating_half_passed"] is True
    text = MD.read_text(encoding="utf-8")
    assert "**FAILED**" in text and "has not been moved" in text


def test_g219_prime_is_stricter_in_kind_not_looser_in_size(art):
    """The replacement must be something a merely-small residual cannot pass: the residual
    has to CONVERGE at the interpolation order."""
    g = art["gates"]["G21_9p_residual_converges_at_the_interpolation_order"]
    assert g["pass"]
    assert 1.7 <= g["observed_order"] <= 2.3
    res = [r["rms_log_residual"] for r in g["ladder"]]
    assert res == sorted(res, reverse=True), "the residual must fall monotonically with h"
    assert res[0] / res[-1] > 10


def test_the_page_says_no_record_carries_a_temperature_axis():
    text = MD.read_text(encoding="utf-8")
    assert "No record in this registry carries a temperature axis" in text
    assert "no default activation energy is shipped" in text.lower()


def test_the_arrhenius_example_on_the_page_is_labelled_as_arithmetic_only(art):
    """A worked activation energy on a page about blood would be read as blood's. The
    artifact and the page both have to say it is not."""
    ex = art["arrhenius_example"]
    assert "not a property of any fluid in this registry" in ex["note"]
    assert "That is not a property of anything in this registry" in MD.read_text(encoding="utf-8")
