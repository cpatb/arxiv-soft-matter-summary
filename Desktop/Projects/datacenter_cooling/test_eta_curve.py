"""The containment CURVE (2026-09-08): gates from PREREG_eta_curve_2026-09-08.md.

Reads the committed artifact rather than re-running the CFD. What these pin is the answer
a buyer needs and the Ra ladder could not give: the SHAPE of phi(eta) between the two
endpoints it sampled.
"""
from __future__ import annotations

import json
import os

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, "cfd", "eta_curve_fit.json")
MD = os.path.join(HERE, "ETA_CURVE_2026-09-08.md")
PREREG = os.path.join(HERE, "PREREG_eta_curve_2026-09-08.md")


@pytest.fixture(scope="module")
def a() -> dict:
    with open(ART, encoding="utf-8") as f:
        return json.load(f)


@pytest.fixture(scope="module")
def md() -> str:
    with open(MD, encoding="utf-8") as f:
        return f.read()


def test_artifacts_exist():
    for p in (ART, MD, PREREG):
        assert os.path.exists(p), p


def test_the_sweep_is_a_sweep_not_two_points(a):
    assert a["n_cells"] == 36
    assert len(a["Ra_rungs"]) == 4
    for r in a["per_ra"]:
        assert r["n_points"] == 9, r


def test_gd1_control_reproduces_the_ra_ladder(a):
    """Without this, nothing below is comparable to the work it extends."""
    g = a["GD1"]
    assert g["verdict"] == "PASS"
    assert g["n_compared"] == 8
    assert g["worst_abs_diff"] <= g["bar"]
    # not vacuous: the two scripts are separate runs, so exact zero would be suspicious
    assert g["worst_abs_diff"] > 0.0


def test_gd2_is_recorded_FAILED_and_the_bar_was_not_moved(a):
    g = a["GD2"]
    assert g["verdict"] == "FAIL"
    assert g["n_rejected"] == 3 and g["n_rungs"] == 4
    assert g["bar_rms_ratio"] == 3.0


def test_gd2_failed_at_the_HIGHEST_Ra_which_is_the_whole_point(a):
    """GD2 demanded rejection on every rung while GD3 predicted the floor falls with Ra.
    The two gates were in tension and GD2 was the mis-specified one; the failure must be
    at the top rung or that reading is wrong."""
    rejected = [p["proportional_rejected"] for p in a["per_ra"]]
    assert rejected[-1] is False
    assert all(rejected[:-1])
    assert a["per_ra"][-1]["Ra"] == max(p["Ra"] for p in a["per_ra"])


def test_gd2_prime_is_stricter_in_kind_and_passes(a):
    g = a["GD2p"]
    assert g["verdict"] == "PASS"
    assert g["monotone_falling"] is True
    assert g["range_factor"] > 4.0
    foa = g["F_over_A_by_Ra"]
    assert foa == sorted(foa, reverse=True), foa


def test_gd3_both_registered_predictions_hold(a):
    g = a["GD3"]
    assert g["verdict"] == "PASS"
    assert g["F_falls_monotonically"] is True
    assert abs(g["A_relative_drop"]) <= g["bar_A_drop"]
    # the floor must actually move, or "it falls" is not a finding
    assert g["F_by_Ra"][0] / g["F_by_Ra"][-1] > 3.0


def test_the_leak_term_is_the_one_that_does_NOT_move(a):
    """The discriminating half: if A moved as much as F, the floor reading would not be
    the explanation. A broken fit that returned constants for both would pass the F test
    alone, so the SIZE of F's move is asserted above and A's stillness here."""
    g = a["GD3"]
    assert abs(g["A_relative_drop"]) < 0.05
    fmove = g["F_by_Ra"][0] / g["F_by_Ra"][-1] - 1.0
    assert fmove > 10 * abs(g["A_relative_drop"])


def test_the_tracer_floor_is_inadmissible_and_the_page_says_so(a, md):
    """F_tracer comes out negative. A recirculation fraction cannot be negative, so this
    is not a floor -- it means the linear+floor form does not fit the tracer at all."""
    g = a["GD4"]
    assert g["any_negative"] is True
    assert all(v < 0 for v in g["F_tracer_by_Ra"])
    assert "must not be quoted as a floor" in md
    assert "not a negative floor" in md
    # and the tracer really does go to ~0, which is why the fit has no floor to find
    assert max(g["tracer_at_max_eta_by_Ra"]) < 0.01


def test_eta_star_is_reported_per_Ra_and_rises_toward_a_hall(a):
    """GD5, the deliverable. No bar was set and none is invented; what is asserted is
    that it is reported per rung and that it moves in the direction the floor implies."""
    es = [p["eta_star"] for p in a["per_ra"]]
    assert len(es) == 4
    assert es == sorted(es), es
    assert 0.9 < min(es) < max(es) < 1.0


def test_the_page_keeps_its_own_scope_limits(md):
    for s in ("Nothing about a real hall", "flux = 0", "not a measured containment"):
        assert s in md, s
