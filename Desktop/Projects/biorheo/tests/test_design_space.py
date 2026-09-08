"""Item 18: the needle-bore exponent and the derating factor.

Two layers, as everywhere else in this repo. The first tests
`biorheo.apps.design_space` against anchors that do not pass through the force solver --
the Newtonian `-4` and the power-law `-(3n+1)`, at parameter values chosen here rather
than read off a record. The second pins the gate verdicts in `data/needle_bore.json`
against the pre-registration `data/needle_bore_prereg.md`.

The discriminating half is asserted separately from the anchors: a routine that returned
`-4` for every fluid would pass G18.1 and G18.2 and mean nothing, so G18.3 is its own
test and is not folded into the gate loop.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from biorheo.apps.design_space import (
    bore_exponent, bore_exponent_closed_form, derating, max_flow_rate, newtonian_rule_error,
)
from biorheo.apps.injectability import Needle, Syringe, injection_force
from biorheo.constitutive.gnf import NEWTONIAN, POWER_LAW
from biorheo.registry import load_record

REPO = Path(__file__).resolve().parent.parent
JSON = REPO / "data" / "needle_bore.json"
MD = REPO / "data" / "needle_bore.md"
RECORDS = REPO / "src" / "biorheo" / "registry" / "records"
SCRIPT = REPO / "scripts" / "needle_bore.py"

NEEDLE = Needle(radius=1.05e-4, length=12.7e-3)
SYRINGE = Syringe(barrel_radius=2.3e-3)
Q = 1.0e-6 / 60.0


@pytest.fixture(scope="module")
def art() -> dict:
    assert JSON.exists(), f"{JSON} missing; run scripts/needle_bore.py"
    return json.loads(JSON.read_text(encoding="utf-8"))


# --- anchors outside the force solver ------------------------------------------------

@pytest.mark.parametrize("mu", [1e-3, 1.0, 37.0])
@pytest.mark.parametrize("Q_test", [1e-9, 1.0e-6 / 60.0, 1e-6])
def test_newtonian_bore_exponent_is_exactly_minus_four(mu, Q_test):
    """`dp = 8 mu L Q / (pi R^4)`, so `b = -4` for every viscosity, rate and length."""
    b = bore_exponent(Q_test, NEEDLE, SYRINGE, NEWTONIAN, np.array([mu]))
    assert b == pytest.approx(-4.0, abs=1e-6)


@pytest.mark.parametrize("n", [0.2, 0.45, 0.7, 1.0])
def test_power_law_bore_exponent_matches_the_closed_form(n):
    """`dp ~ Q^n R^-(3n+1)`. The closed form is evaluated from `n` alone and touches
    neither the quadrature nor the pressure inversion."""
    params = np.array([2.0, n])   # k, n
    b = bore_exponent(Q, NEEDLE, SYRINGE, POWER_LAW, params)
    assert bore_exponent_closed_form(POWER_LAW, params) == pytest.approx(-(3 * n + 1))
    assert b == pytest.approx(-(3 * n + 1), rel=1e-6)


def test_power_law_at_n_equal_one_returns_the_newtonian_value():
    """A discriminating negative case for the closed form: at `n = 1` the two anchors must
    agree, and at `n = 0.45` they must NOT."""
    assert bore_exponent_closed_form(POWER_LAW, np.array([2.0, 1.0])) == -4.0
    assert bore_exponent_closed_form(POWER_LAW, np.array([2.0, 0.45])) != -4.0


def test_closed_form_returns_none_where_there_is_none():
    """It must never fall back on the numerical route; its whole purpose is to be
    independent of it."""
    rec = load_record(RECORDS / "blood_horner2020_donor_A.yaml")
    assert bore_exponent_closed_form(rec.model, rec.params_si) is None


def test_shear_thinning_is_never_steeper_than_newtonian():
    """The one-sided physical claim, on a model built here rather than loaded."""
    for n in (0.3, 0.5, 0.8, 1.0):
        b = bore_exponent(Q, NEEDLE, SYRINGE, POWER_LAW, np.array([2.0, n]))
        assert b >= -4.0 - 1e-9, n


def test_newtonian_rule_error_is_zero_for_a_newtonian_fluid():
    """The rule IS the fluid's own scaling there, so the over-prediction must vanish."""
    out = newtonian_rule_error(Q, NEEDLE, SYRINGE, NEWTONIAN, np.array([3.5e-3]),
                               bore_ratio=0.7)
    assert abs(out["over_prediction"]) < 1e-9
    assert out["F_newtonian_rule"] == pytest.approx(out["F_true"], rel=1e-9)


def test_newtonian_rule_over_predicts_for_a_shear_thinning_fluid():
    """The discriminating half of the previous test: the routine must not return zero for
    everything."""
    out = newtonian_rule_error(Q, NEEDLE, SYRINGE, POWER_LAW, np.array([2.0, 0.45]),
                               bore_ratio=0.7)
    assert out["over_prediction"] > 0.1


def test_newtonian_rule_error_refuses_a_widening_bore():
    with pytest.raises(ValueError, match="bore_ratio"):
        newtonian_rule_error(Q, NEEDLE, SYRINGE, NEWTONIAN, np.array([1e-3]), bore_ratio=1.5)


# --- the budget crossing --------------------------------------------------------------

def test_max_flow_rate_inverts_the_newtonian_closed_form():
    """`F = 8 mu L Q A_s / (pi R^4)`, so `Q_max = budget pi R^4 / (8 mu L A_s)`, computed
    here without the root finder."""
    mu, budget = 5e-3, 25.0
    expected = (budget * np.pi * NEEDLE.radius ** 4
                / (8 * mu * NEEDLE.length * SYRINGE.area))
    got = max_flow_rate(budget, NEEDLE, SYRINGE, NEWTONIAN, np.array([mu]))
    assert got == pytest.approx(expected, rel=1e-8)


def test_max_flow_rate_raises_rather_than_clipping():
    """A budget the bracket cannot straddle is a finding, not a number to clip."""
    with pytest.raises(ValueError, match="never reaches"):
        max_flow_rate(1e12, NEEDLE, SYRINGE, NEWTONIAN, np.array([1e-3]))
    with pytest.raises(ValueError, match="unreachable from below"):
        max_flow_rate(1e-12, NEEDLE, SYRINGE, NEWTONIAN, np.array([1e-3]))


def test_derating_with_zero_width_intervals_is_exactly_one():
    """G18.6-prime's first half. Without this, a routine returning a constant 0.5 would
    pass the pre-registered `factor < 1` bar for every record."""
    p = np.array([5e-3])
    d = derating(25.0, NEEDLE, SYRINGE, NEWTONIAN, p, {"mu": (float(p[0]), float(p[0]))})
    assert d.factor == pytest.approx(1.0, abs=1e-12)
    assert d.q_upper_bound == pytest.approx(d.q_nominal, rel=1e-12)


def test_derating_with_a_real_interval_is_strictly_below_one():
    """G18.6-prime's second half."""
    p = np.array([5e-3])
    d = derating(25.0, NEEDLE, SYRINGE, NEWTONIAN, p, {"mu": (4e-3, 6e-3)})
    assert d.factor < 1.0
    assert d.binding_parameter == "mu" and d.binding_end == "hi"
    # Newtonian: Q_max ~ 1/mu exactly, so the factor is the viscosity ratio
    assert d.factor == pytest.approx(5e-3 / 6e-3, rel=1e-8)


def test_derating_refuses_when_no_interval_is_usable():
    with pytest.raises(ValueError, match="no usable profile end"):
        derating(25.0, NEEDLE, SYRINGE, NEWTONIAN, np.array([5e-3]), {})


def test_bore_exponent_excludes_friction_by_construction():
    """`b` must be a property of the fluid, so a syringe friction term cannot enter it.
    The signature has no `friction`; this pins that the force it differentiates is the
    frictionless one."""
    p = np.array([5e-3])
    f_no = injection_force(Q, NEEDLE, SYRINGE, NEWTONIAN, p).force
    f_yes = injection_force(Q, NEEDLE, SYRINGE, NEWTONIAN, p, friction=10.0).force
    assert f_yes == pytest.approx(f_no + 10.0)
    assert bore_exponent(Q, NEEDLE, SYRINGE, NEWTONIAN, p) == pytest.approx(-4.0, abs=1e-6)


# --- the shipped artifact -------------------------------------------------------------

def test_every_gate_passes(art):
    for name, g in art["gates"].items():
        if "pass" in g:
            assert g["pass"], (name, g)


def test_g181_newtonian_control(art):
    g = art["gates"]["G18_1_newtonian_control"]
    assert g["n"] == 20 and g["n_pass"] == 20
    assert g["worst"] < 1e-3


def test_g182_power_law_closed_form(art):
    g = art["gates"]["G18_2_power_law_closed_form"]
    assert g["n"] == 20 and g["n_pass"] == 20
    assert g["worst"] < 1e-3


def test_g183_is_the_discriminating_half_and_is_not_folded_into_the_anchors(art):
    """Asserted on its own: G18.1 and G18.2 are both passed by a routine that returns -4
    for every fluid, and this is the test that such a routine fails."""
    g = art["gates"]["G18_3_discriminating_half"]
    assert g["worst"] > 0.1
    assert g["n_above"] >= 1
    assert g["worst_record"].startswith("carbopol")


def test_g184_no_record_is_steeper_than_newtonian(art):
    g = art["gates"]["G18_4_never_steeper_than_newtonian"]
    assert g["n"] == 115 and not g["failures"]
    assert g["steepest"] >= -4.0 - 1e-6


def test_g185_split_prediction_came_out_as_pre_registered(art):
    """Both halves, and the uninformative flag, which is the outcome the repo has had to
    record before when a comparison decides nothing."""
    g = art["gates"]["G18_5_newtonian_rule_error_split"]
    assert g["blood_half_as_predicted"] and g["carbopol_half_as_predicted"]
    assert g["blood_n_below"] == g["blood_n"] == 20
    assert g["blood_worst"] < 0.10
    assert g["carbopol_all_above"] == g["carbopol_all_n"] == 15
    assert g["carbopol_least"] > 0.10
    assert not g["uninformative"]


def test_g186_prime_pins_both_ends_the_pre_registered_bar_cannot(art):
    g = art["gates"]["G18_6p_derating_negative_control"]
    assert g["factor_zero_width"] == pytest.approx(1.0, abs=1e-12)
    assert g["factor_real"] < 1.0
    assert "unfalsifiable" in g["note"]


def test_the_identifiability_split_is_a_tail_not_a_shift(art):
    """The medians must NOT separate and the worst cases must. Asserting the first half is
    the point: the section was first written claiming a shift, and the rows said otherwise."""
    g = art["gates"]["G18_6_derating"]
    med_a = g["median_factor_all_measurable"]
    med_b = g["median_factor_not_all_measurable"]
    assert abs(med_a - med_b) < 0.2, (med_a, med_b)
    assert g["worst_factor_not_all_measurable"] < g["worst_factor_all_measurable"] / 100


def test_no_gauge_or_standard_is_named_on_the_page():
    """Standing repo rule: the ISO 9626 comparison was retracted and the standard is still
    not in hand, so this item's page must not name a gauge."""
    text = MD.read_text(encoding="utf-8")
    for banned in ("27G", "30G", "27 G", "30 G", "ISO 9626:"):
        assert banned not in text, banned


def test_the_page_says_the_registry_has_no_mab():
    """The venture's first wedge fluid is not in this data and the page must say so, so a
    reader cannot mistake a blood demonstration for a formulation result."""
    text = MD.read_text(encoding="utf-8")
    assert "monoclonal antibody" in text and "registry contains none" in text.replace("**", "")


def test_every_record_in_the_registry_is_non_commercial(art):
    """**All 115 records IN THIS ARTIFACT'S FROZEN LIST are `commercial_use: false`** -- the 15 Carbopol ones under
    CC BY-NC-SA, the 100 blood ones because the source workbook states no licence at all.

    This is not a detail of this item; it is a standing constraint on the venture, because
    a paid deliverable cannot be built on any curve currently in the registry. The page has
    to say so, and this test fails the day someone adds a record and forgets to re-check.
    """
    n_nc = sum(1 for r in art["rows"] if not r["commercial_ok"])
    assert n_nc == len(art["rows"]) == 115
    text = MD.read_text(encoding="utf-8")
    assert "commercial_use: false" in text
    assert "cannot appear in a commercial deliverable" in text


@pytest.mark.slow
def test_rerunning_the_script_reproduces_every_gate(art, tmp_path):
    """RESTORES both artifacts before asserting, so a failure cannot destroy the reference
    it was checking against."""
    import shutil
    import subprocess
    import sys
    saved = {p: tmp_path / p.name for p in (JSON, MD)}
    for src, dst in saved.items():
        shutil.copy2(src, dst)
    try:
        r = subprocess.run([sys.executable, str(SCRIPT)], capture_output=True, text=True,
                           timeout=3600)
        assert r.returncode == 0, r.stdout[-3000:] + r.stderr[-2000:]
        new = json.loads(JSON.read_text(encoding="utf-8"))
    finally:
        for src, dst in saved.items():
            shutil.copy2(dst, src)
    for name, g in art["gates"].items():
        if "pass" in g:
            assert new["gates"][name]["pass"] == g["pass"], name
    for a, b in zip(art["rows"], new["rows"]):
        assert a["id"] == b["id"]
        assert b["b_Q1"] == pytest.approx(a["b_Q1"], rel=1e-9), a["id"]
