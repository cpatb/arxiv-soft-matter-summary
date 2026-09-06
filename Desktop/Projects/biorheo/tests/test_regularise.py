"""Unit 0.3b: regularisation schemes carry their knob and report their own error.

Hand anchors (Bingham: tau_y = 1 Pa, k = 1 Pa s, n = 1; exact eta = 1/gd + 1):
  papanastasiou, m = 1 s, gd = 1:  (1 - e^-1)/1 + 1 = 0.6321206 + 1 = 1.6321206
                                   exact 2, induced error (1.6321206 - 2)/2 = -0.1839397
  bercovier_engelman, eps = 1, gd = 1: 1/sqrt(2) + 1 = 1.7071068; error -0.1464466
  clamp, eta_max = 10: exact at gd = 0.05 is 21 -> 10; at gd = 1 is 2 -> 2
                       bi-viscosity crossing 10 gd_c = 1 + gd_c -> gd_c = 1/9
"""

from __future__ import annotations

import numpy as np
import pytest

from biorheo.constitutive import (
    CASSON,
    CROSS,
    HERSCHEL_BULKLEY,
    SCHEMES,
    critical_rate,
    induced_error,
    regularise,
)

BINGHAM = (1.0, 1.0, 1.0)  # tau_y, k, n


# --- construction rules ---------------------------------------------------------------

def test_knob_is_mandatory_and_must_be_the_right_one():
    with pytest.raises(ValueError, match="exactly one knob"):
        regularise(HERSCHEL_BULKLEY, "papanastasiou")
    with pytest.raises(ValueError, match="exactly one knob"):
        regularise(HERSCHEL_BULKLEY, "papanastasiou", eps=1.0)
    with pytest.raises(ValueError, match="positive"):
        regularise(HERSCHEL_BULKLEY, "papanastasiou", m=0.0)
    with pytest.raises(ValueError, match="unknown scheme"):
        regularise(HERSCHEL_BULKLEY, "biviscosity", eta_max=1.0)


def test_regularised_model_names_its_scheme_and_knob():
    r = regularise(HERSCHEL_BULKLEY, "papanastasiou", m=100.0)
    assert r.name == "herschel_bulkley+papanastasiou"
    assert r.fixed == {"m": 100.0}
    assert r.regularisation == {"papanastasiou": 100.0}
    assert r.is_regularised and r.base is HERSCHEL_BULKLEY
    assert r.param_names == HERSCHEL_BULKLEY.param_names
    assert "m" in r.accepts_fixed()
    with pytest.raises(ValueError, match="already regularised"):
        regularise(r, "clamp", eta_max=1.0)


def test_yield_schemes_refuse_models_without_a_yield_term():
    with pytest.raises(ValueError, match="yield term"):
        regularise(CROSS, "papanastasiou", m=1.0)
    # but the clamp applies to anything
    r = regularise(CROSS, "clamp", eta_max=0.5)
    np.testing.assert_allclose(r(np.array([1e-6, 1e6]), 1.0, 0.1, 1.0, 1.0), [0.5, 0.1],
                               rtol=1e-5)


def test_every_scheme_has_a_named_knob_with_a_unit():
    assert SCHEMES == {"papanastasiou": ("m", "s"), "bercovier_engelman": ("eps", "1/s"),
                       "clamp": ("eta_max", "Pa*s")}


# --- hand values ----------------------------------------------------------------------

def test_papanastasiou_hand_value():
    r = regularise(HERSCHEL_BULKLEY, "papanastasiou", m=1.0)
    assert r(1.0, *BINGHAM) == pytest.approx(1.6321206, rel=1e-7)
    assert induced_error(r, 1.0, BINGHAM) == pytest.approx(-0.1839397, rel=1e-6)


def test_bercovier_engelman_hand_value():
    r = regularise(HERSCHEL_BULKLEY, "bercovier_engelman", eps=1.0)
    assert r(1.0, *BINGHAM) == pytest.approx(1.7071068, rel=1e-7)
    assert induced_error(r, 1.0, BINGHAM) == pytest.approx(-0.1464466, rel=1e-6)


def test_clamp_hand_values():
    r = regularise(HERSCHEL_BULKLEY, "clamp", eta_max=10.0)
    np.testing.assert_allclose(r(np.array([0.05, 1.0]), *BINGHAM), [10.0, 2.0])
    assert induced_error(r, 0.05, BINGHAM) == pytest.approx((10 - 21) / 21)
    assert induced_error(r, 1.0, BINGHAM) == pytest.approx(0.0)


# --- limits -----------------------------------------------------------------------------

def test_papanastasiou_is_finite_at_zero_rate_for_bingham():
    r = regularise(HERSCHEL_BULKLEY, "papanastasiou", m=50.0)
    # eta(0) = tau_y m + k = 50 + 1
    assert r(0.0, *BINGHAM) == pytest.approx(51.0)
    assert r(1e-9, *BINGHAM) == pytest.approx(51.0, rel=1e-6)
    assert induced_error(r, 0.0, BINGHAM) == -1.0


def test_papanastasiou_still_diverges_for_n_below_one():
    # The scheme regularises the yield term only; k gd^(n-1) is untouched.
    r = regularise(HERSCHEL_BULKLEY, "papanastasiou", m=50.0)
    p = (1.0, 1.0, 0.5)
    assert r(1e-8, *p) > 1e3
    assert np.isinf(r(0.0, *p))
    c = regularise(HERSCHEL_BULKLEY, "clamp", eta_max=1e3)
    assert c(0.0, *p) == 1e3  # only the clamp bounds the whole viscosity


def test_bercovier_engelman_zero_rate_limit():
    r = regularise(HERSCHEL_BULKLEY, "bercovier_engelman", eps=0.01)
    assert r(0.0, *BINGHAM) == pytest.approx(1.0 / 0.01 + 1.0)


def test_casson_papanastasiou_zero_rate_limit():
    r = regularise(CASSON, "papanastasiou", m=10.0)
    tau_y, eta_c = 0.004, 0.0036
    # (sqrt(tau_y m) + sqrt(eta_c))^2 = (sqrt(0.04) + 0.06)^2 = (0.2 + 0.06)^2 = 0.0676
    assert r(0.0, tau_y, eta_c) == pytest.approx(0.0676, rel=1e-9)


def test_all_schemes_converge_to_the_exact_law_as_the_knob_sharpens():
    gd = np.logspace(-2, 2, 9)
    p = (2.0, 0.5, 0.6)
    exact = HERSCHEL_BULKLEY(gd, *p)
    for scheme, knob in (("papanastasiou", "m"), ("bercovier_engelman", "eps"),
                         ("clamp", "eta_max")):
        # papanastasiou: m up; bercovier: eps down; clamp: eta_max up
        # exact eta at the grid's lowest rate, gd = 0.01, is 200 + 0.5 * 0.01^-0.4 = 203,
        # so the clamp is inert (error exactly 0) from eta_max = 1e3 upward
        values = {"m": (1.0, 1e2, 1e4), "eps": (1.0, 1e-2, 1e-4),
                  "eta_max": (1e1, 1e2, 1e3)}[knob]
        errs = [np.max(np.abs(induced_error(regularise(HERSCHEL_BULKLEY, scheme,
                                                       **{knob: v}), gd, p)))
                for v in values]
        assert errs[0] > errs[1] > errs[2], (scheme, errs)
        np.testing.assert_allclose(
            regularise(HERSCHEL_BULKLEY, scheme, **{knob: values[-1]})(gd, *p), exact,
            rtol=1e-3)


def test_regularisation_never_overpredicts():
    gd = np.logspace(-3, 3, 25)
    p = (2.0, 0.5, 0.6)
    for scheme, knob in (("papanastasiou", {"m": 5.0}), ("bercovier_engelman", {"eps": 0.2}),
                         ("clamp", {"eta_max": 50.0})):
        assert np.all(induced_error(regularise(HERSCHEL_BULKLEY, scheme, **knob), gd, p) <= 0)


# --- the bi-viscosity identity ---------------------------------------------------------

def test_clamp_equals_biviscosity_stress_form():
    # Bi-viscosity: sigma = eta_max gd below gd_c, tau_y + k gd^n above, with gd_c fixed
    # by stress continuity. Built here independently of the clamp and compared on a
    # grid straddling gd_c = 1/9 for the Bingham anchor.
    from scipy.optimize import brentq

    tau_y, k, n = 1.0, 1.0, 0.7
    eta_max = 10.0
    gd_c = brentq(lambda g: eta_max * g - (tau_y + k * g ** n), 1e-9, 1e3)
    gd = np.logspace(-3, 1, 41)
    sigma_bv = np.where(gd < gd_c, eta_max * gd, tau_y + k * gd ** n)
    r = regularise(HERSCHEL_BULKLEY, "clamp", eta_max=eta_max)
    np.testing.assert_allclose(r(gd, tau_y, k, n) * gd, sigma_bv, rtol=1e-12)
    # and the Bingham crossing is 1/9 by hand
    assert brentq(lambda g: 10 * g - (1 + g), 1e-9, 1e3) == pytest.approx(1 / 9)


# --- critical rate -----------------------------------------------------------------------

def test_critical_rate_bingham_papanastasiou():
    # |err| = e^{-gd}/(1 + gd) = 0.01 at gd = 3.1759 (hand bracket: 3.175 -> 0.010011,
    # 3.176 -> 0.009998).
    r = regularise(HERSCHEL_BULKLEY, "papanastasiou", m=1.0)
    gc = critical_rate(r, BINGHAM, rtol=0.01)
    assert 3.175 < gc < 3.176
    assert abs(induced_error(r, gc, BINGHAM)) == pytest.approx(0.01, rel=1e-6)
    assert abs(induced_error(r, 0.9 * gc, BINGHAM)) > 0.01


def test_critical_rate_refuses_when_uninformative():
    # m = 1e-9 s: |err| = e^{-1e-9 gd}/(1 + gd) is 1e-6 at the top of the grid, so a
    # 1e-7 tolerance is never met.
    r = regularise(HERSCHEL_BULKLEY, "papanastasiou", m=1e-9)
    with pytest.raises(ValueError, match="never drops"):
        critical_rate(r, BINGHAM, rtol=1e-7)
    # a tolerance below what a double can resolve in a relative difference is refused
    with pytest.raises(ValueError, match="resolve"):
        critical_rate(r, BINGHAM, rtol=1e-13)
    with pytest.raises(ValueError, match="below rtol"):
        critical_rate(regularise(HERSCHEL_BULKLEY, "clamp", eta_max=1e12), BINGHAM, rtol=0.5)


def test_yield_schemes_refuse_a_base_model_with_fixed_constants():
    """The Papanastasiou and Bercovier-Engelman schemes rebuild the law from the model's
    parameters and never see `base.fixed`, so a fixed constant would be dropped without
    a word. The clamp wraps `base.func` and does receive them, so it is still allowed.
    No shipped model has a fixed constant; this is the guard that keeps that a fact
    rather than an assumption (NEXT_STEPS item 5, referee item, closed 2026-09-05)."""
    from dataclasses import replace as _replace

    from biorheo.constitutive import get_model
    from biorheo.constitutive.regularise import regularise

    hb = get_model("herschel_bulkley")
    assert not hb.fixed, "a shipped model grew a fixed constant; re-read this test"
    with_fixed = _replace(hb, fixed={"rho": 1060.0})
    for scheme, knob in (("papanastasiou", {"m": 1e3}),
                         ("bercovier_engelman", {"eps": 1e-3})):
        with pytest.raises(ValueError, match="rho"):
            regularise(with_fixed, scheme, **knob)
    # the clamp is unaffected: it wraps base.func, which does get base.fixed
    regularise(hb, "clamp", eta_max=10.0)
