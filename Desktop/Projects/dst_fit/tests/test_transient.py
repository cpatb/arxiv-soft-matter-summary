"""Tests for the strain clock: the Wyart-Cates inversion and the recovery fit.

The load-bearing test here is
:func:`test_inversion_is_invariant_along_the_degenerate_direction`. The whole
reason ``gamma_0`` is worth measuring through the model rather than by
landmarking the viscosity is the claim that the inversion depends only on the
ratios ``phi_0/phi`` and ``phi_m/phi``, which are exactly what the package's
central degeneracy preserves. That claim is proved to machine precision below
rather than argued in prose.

Everything else is checked against synthetic data with a known answer, because
on the real records there is nothing to check against.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dst_fit.io_literature import LIN_TRANSIENTS, read_lin_transients
from dst_fit.models import PHI_0_RCP, WC_EXPONENT, WYART_CATES, frictional_fraction
from dst_fit.transient import (
    MIN_ADMISSIBLE_FRACTION,
    SMOOTH_MIN_SAMPLES,
    admissible_window,
    fit_recovery,
    fit_strain_clock,
    frictional_fraction_from_viscosity,
    frictionless_relative_viscosity,
    half_recovery_strain,
    smooth,
    smoothing_samples,
)

ROOT = Path(__file__).resolve().parents[1]
LIN_ROOT = ROOT / (
    "Hydrodynamic and Contact Contributions to Continuous Shear Thickening "
    "in Colloidal Suspensions"
)
needs_lin = pytest.mark.skipif(
    not (LIN_ROOT / "Fig2" / "Fig2").is_dir(), reason="Lin et al. deposit not present"
)


def wc_relative_viscosity(f, phi, phi_m, phi_0, n=WC_EXPONENT):
    """eta_r at a given frictional fraction. The forward map the inversion undoes."""
    f = np.asarray(f, dtype=float)
    phi_j = (1.0 - f) * phi_0 + f * phi_m
    return (1.0 - phi / phi_j) ** -float(n)


# --------------------------------------------------------------------------
# the inversion
# --------------------------------------------------------------------------
def test_inversion_round_trips_the_forward_model():
    """f -> eta_r -> f returns the original to machine precision."""
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    f = np.linspace(0.0, 1.0, 41)
    eta = wc_relative_viscosity(f, phi, phi_m, phi_0)
    back = frictional_fraction_from_viscosity(eta, phi, phi_m, phi_0)
    np.testing.assert_allclose(back, f, rtol=0, atol=1e-12)


def test_inversion_round_trips_at_a_non_default_exponent():
    """The exponent has to travel with the fit; n = 2.2 is Guy et al.'s value."""
    phi, phi_m, phi_0, n = 0.45, 0.558, 0.611, 2.2
    f = np.linspace(0.05, 0.95, 19)
    eta = wc_relative_viscosity(f, phi, phi_m, phi_0, n=n)
    back = frictional_fraction_from_viscosity(eta, phi, phi_m, phi_0, n=n)
    np.testing.assert_allclose(back, f, rtol=0, atol=1e-12)


def test_inversion_matches_the_package_forward_model():
    """Cross-check against models._wyart_cates rather than a local reimplementation."""
    phi, phi_m, phi_0, sigma_star = 0.40, 0.52, 0.64, 5.0
    model = WYART_CATES.with_fixed(eta_s=1.0, phi_0=phi_0)
    sigma = np.logspace(-1, 3, 60)
    eta_r = model(sigma, phi, phi_m, sigma_star)
    expected = frictional_fraction(sigma, sigma_star)
    got = frictional_fraction_from_viscosity(eta_r, phi, phi_m, phi_0)
    np.testing.assert_allclose(got, expected, rtol=1e-10, atol=1e-12)


def test_inversion_is_invariant_along_the_degenerate_direction():
    """The claim that makes gamma_0 quotable where phi_m is not.

    phi, phi_m and phi_0 enter Wyart-Cates only through phi/phi_J, so scaling
    all three by a common factor leaves every viscosity, and every inverted f,
    unchanged. That is the direction the fitted (phi, phi_m) valley runs along:
    a flow curve pins the ratios and neither number. The inversion is a function
    of the ratios alone, so f is constant along it -- to machine precision, not
    approximately.
    """
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    eta = np.logspace(np.log10(1.01), 2.0, 200)
    reference = frictional_fraction_from_viscosity(eta, phi, phi_m, phi_0)

    for scale in (0.90, 0.95, 1.0, 1.05, 1.10):
        moved = frictional_fraction_from_viscosity(
            eta, phi * scale, phi_m * scale, phi_0 * scale
        )
        np.testing.assert_allclose(moved, reference, rtol=0, atol=1e-12)


def test_inversion_is_not_invariant_to_a_change_of_shape():
    """The invariance is along one direction only, not a claim that nothing matters.

    Moving phi_m while holding phi and phi_0 changes the ratio phi_m/phi and must
    change f. Without this, the test above would also pass on a function that
    ignored its arguments.
    """
    eta = np.logspace(np.log10(1.5), 1.5, 50)
    a = frictional_fraction_from_viscosity(eta, 0.40, 0.52, 0.64)
    b = frictional_fraction_from_viscosity(eta, 0.40, 0.50, 0.64)
    assert np.max(np.abs(a - b)) > 0.05


def test_inversion_goes_negative_below_the_frictionless_branch_rather_than_clipping():
    """A dip below the model's own floor is reported, not hidden."""
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    floor = frictionless_relative_viscosity(phi, phi_0)
    below = frictional_fraction_from_viscosity(np.array([0.5 * floor]), phi, phi_m, phi_0)
    assert below[0] < 0.0
    at = frictional_fraction_from_viscosity(np.array([floor]), phi, phi_m, phi_0)
    np.testing.assert_allclose(at, [0.0], atol=1e-12)


def test_inversion_exceeds_one_above_the_frictional_branch():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    ceiling = (1.0 - phi / phi_m) ** -2.0
    above = frictional_fraction_from_viscosity(np.array([1.5 * ceiling]), phi, phi_m, phi_0)
    assert above[0] > 1.0


def test_inversion_returns_nan_at_or_below_the_solvent_viscosity():
    got = frictional_fraction_from_viscosity(
        np.array([0.5, 1.0, 1.0 - 1e-15, 2.0]), 0.40, 0.52, 0.64
    )
    assert np.isnan(got[:3]).all()
    assert np.isfinite(got[3])


def test_inversion_propagates_nan_and_handles_empty():
    got = frictional_fraction_from_viscosity(np.array([np.nan, np.inf, 5.0]), 0.4, 0.52, 0.64)
    assert np.isnan(got[0])
    assert np.isfinite(got[2])
    assert frictional_fraction_from_viscosity(np.array([]), 0.4, 0.52, 0.64).size == 0


@pytest.mark.parametrize(
    "phi, phi_m, phi_0",
    [(0.52, 0.40, 0.64), (0.40, 0.70, 0.64), (-0.1, 0.52, 0.64), (0.40, 0.52, 0.52)],
)
def test_inversion_refuses_a_bad_packing_ordering(phi, phi_m, phi_0):
    with pytest.raises(ValueError, match="phi"):
        frictional_fraction_from_viscosity(np.array([5.0]), phi, phi_m, phi_0)


def test_frictionless_viscosity_agrees_with_the_model_at_f_zero():
    phi, phi_0 = 0.40, PHI_0_RCP
    model = WYART_CATES.with_fixed(eta_s=1.0, phi_0=phi_0)
    # sigma far below sigma* drives f to zero, so the model tends to the floor.
    at_zero = float(model(np.array([1e-30]), phi, 0.52, 5.0)[0])
    np.testing.assert_allclose(
        frictionless_relative_viscosity(phi, phi_0), at_zero, rtol=1e-12
    )


# --------------------------------------------------------------------------
# smoothing
# --------------------------------------------------------------------------
def test_smoothing_width_tracks_the_sample_spacing():
    """Two records sampled differently must be smoothed over the same strain."""
    fine = np.arange(0.0, 10.0, 5e-4)
    coarse = np.arange(0.0, 10.0, 1e-3)
    assert smoothing_samples(fine, window=0.05) == 100
    assert smoothing_samples(coarse, window=0.05) == 50


def test_smoothing_width_has_a_floor():
    assert smoothing_samples(np.array([0.0, 1.0]), window=0.05) == SMOOTH_MIN_SAMPLES


def test_smoothing_attenuates_an_isolated_spike():
    """The ordinary job: a few bad samples must not survive as a landmark."""
    y = np.full(2000, 10.0)
    y[1000:1005] = 2.0
    out = smooth(y, size=101)
    assert out.min() > 9.5  # 5 samples of 8 deficit spread over 101
    assert np.isclose(out[0], 10.0) and np.isclose(out[-1], 10.0)


def test_the_reversal_minimum_is_the_interior_dip_not_an_end_artefact():
    """The defect this module was built around, guarded end to end.

    Several Lin records fall over their last handful of samples. Taken as a
    global argmin that end artefact becomes "the reversal minimum", which puts
    the fit window past the end of the recovery. A boxcar alone does not fix it:
    smoothing with ``mode="nearest"`` extends the final value, so a terminal
    spike is attenuated but not erased. The interior restriction is what fixes
    it, and this checks the two working together.
    """
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.05, 0.80, phi, phi_m, phi_0)
    floor = frictionless_relative_viscosity(phi, phi_0)
    eta[-6:] = 0.10 * floor  # deeper than any interior point, at the very end

    win = admissible_window(g, eta, phi, phi_m, phi_0)
    assert win.gamma_min < 1.0, "the minimum was taken from the end-of-record artefact"
    rec = fit_recovery(g, eta, 2.0, phi, phi_m, phi_0)
    np.testing.assert_allclose(rec.gamma_eff, 0.35, rtol=0.05)


def test_smooth_passes_short_records_through():
    y = np.array([1.0, 2.0, 3.0])
    np.testing.assert_array_equal(smooth(y, size=101), y)


# --------------------------------------------------------------------------
# the recovery fit
# --------------------------------------------------------------------------
def synthetic_transient(gamma_eff, f_0, f_ss, phi, phi_m, phi_0, seed=0, noise=0.0):
    """A record built by running the model forwards, for tests with a known answer."""
    gamma = np.arange(5e-4, 10.0, 5e-4)
    f = f_ss + (f_0 - f_ss) * np.exp(-gamma / gamma_eff)
    eta = wc_relative_viscosity(f, phi, phi_m, phi_0)
    if noise > 0:
        rng = np.random.default_rng(seed)
        eta = eta * (1.0 + noise * rng.standard_normal(eta.size))
    return gamma, eta


def test_recovery_recovers_a_known_strain_constant():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0)
    rec = fit_recovery(g, eta, 2.0, phi, phi_m, phi_0)
    assert rec.success
    np.testing.assert_allclose(rec.gamma_eff, 0.35, rtol=0.02)
    np.testing.assert_allclose(rec.f_plateau, 0.80, rtol=0.02)


def test_recovery_survives_multiplicative_noise():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0, seed=7, noise=0.02)
    rec = fit_recovery(g, eta, 2.0, phi, phi_m, phi_0)
    np.testing.assert_allclose(rec.gamma_eff, 0.35, rtol=0.10)


def test_recovery_is_invariant_along_the_degenerate_direction():
    """The end-to-end version of the inversion invariance, on a fitted constant."""
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0, seed=3, noise=0.01)
    base = fit_recovery(g, eta, 2.0, phi, phi_m, phi_0).gamma_eff
    for scale in (0.92, 1.0, 1.08):
        moved = fit_recovery(g, eta, 2.0, phi * scale, phi_m * scale, phi_0 * scale).gamma_eff
        np.testing.assert_allclose(moved, base, rtol=1e-9)


def test_recovery_reports_the_time_constant_consistently():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0)
    rec = fit_recovery(g, eta, 4.0, phi, phi_m, phi_0)
    np.testing.assert_allclose(rec.t_eff, rec.gamma_eff / 4.0, rtol=1e-12)


def test_admissible_window_refuses_a_record_that_never_reaches_the_model():
    """A record entirely below the frictionless branch has no invertible part.

    This is where the "the model cannot describe this" contract lives. It is
    :func:`admissible_window` that raises, not :func:`fit_recovery`: see the next
    test for why the two were deliberately separated.
    """
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    floor = frictionless_relative_viscosity(phi, phi_0)
    g = np.arange(5e-4, 10.0, 5e-4)
    eta = np.full(g.size, 0.6 * floor)
    eta[-50:] = 0.9 * floor  # still below the floor everywhere
    with pytest.raises(ValueError, match="frictionless branch"):
        admissible_window(g, eta, phi, phi_m, phi_0)


def test_recovery_still_measures_a_clock_below_the_frictionless_floor():
    """The clock survives where the f-inversion does not, and that is the point.

    ``fit_recovery`` fits ``u = phi_J/phi``, which is defined and monotone for
    every ``eta_r > 1``, so it returns a clock even when the whole record sits
    under the model's frictionless floor and ``f`` is negative throughout. On
    Lin's PMMA at ``gamma_dot = 0.05`` that is the actual situation, and
    refusing there would discard the only record at the low end of the rate
    ladder -- which is the end the rate trend is measured against.

    What is *not* licensed below the floor is reading the fitted amplitude as a
    physical frictional fraction. The clock is invariant to the affine map
    between ``f`` and ``u``; the amplitudes are not.
    """
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    floor = frictionless_relative_viscosity(phi, phi_0)
    g = np.arange(5e-4, 10.0, 5e-4)
    lam = 0.8
    # The recovery is built exponential in u, not in eta: u is the variable the
    # clock is defined on, and eta is a nonlinear function of it. Building the
    # synthetic in eta instead and demanding the same lambda back would be
    # testing the fitter against a curve that is not the model (it recovers
    # 0.713 for a nominal 0.800, an 11% bias that is real, not numerical).
    u = 1.690 + (2.366 - 1.690) * np.exp(-g / lam)
    eta = (1.0 - 1.0 / u) ** -2.0
    assert eta.max() < floor, "this synthetic is meant to sit under the floor"
    rec = fit_recovery(g, eta, 1.0, phi, phi_m, phi_0)
    assert rec.window.depth_below_floor < 1.0
    assert rec.window.admissible_fraction == 0.0
    np.testing.assert_allclose(rec.gamma_eff, lam, rtol=0.02)


def test_recovery_rejects_a_non_positive_rate():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0)
    with pytest.raises(ValueError, match="shear rate"):
        fit_recovery(g, eta, 0.0, phi, phi_m, phi_0)


def test_window_reports_the_depth_below_the_floor():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0)
    floor = frictionless_relative_viscosity(phi, phi_0)
    eta[:200] = 0.5 * floor  # a dip the model cannot reach
    win = admissible_window(g, eta, phi, phi_m, phi_0)
    assert win.depth_below_floor < 1.0
    assert win.admissible_fraction >= MIN_ADMISSIBLE_FRACTION
    assert win.usable


def test_window_refuses_unsorted_strain():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0)
    g = g[::-1].copy()
    with pytest.raises(ValueError, match="increasing"):
        admissible_window(g, eta, phi, phi_m, phi_0)


# --------------------------------------------------------------------------
# separating the two clocks
# --------------------------------------------------------------------------
def test_strain_clock_recovers_known_constants():
    """Build gamma_eff from Richards Eq. 4 and read the two clocks back off it."""
    gamma_0, t_r = 0.20, 0.24
    rates = np.array([0.5, 1.0, 2.0, 5.0])
    gamma_eff = 1.0 / (1.0 / gamma_0 + 1.0 / (rates * t_r))
    clock = fit_strain_clock(rates, gamma_eff)
    assert clock.form_holds
    np.testing.assert_allclose(clock.gamma_0, gamma_0, rtol=1e-9)
    np.testing.assert_allclose(clock.t_r, t_r, rtol=1e-9)
    np.testing.assert_allclose(clock.r_squared, 1.0, atol=1e-12)


def test_strain_clock_flags_a_falling_gamma_eff_as_outside_the_form():
    """gamma_eff decreasing with rate cannot come from an additive time clock."""
    rates = np.array([0.5, 1.0, 2.0, 5.0])
    gamma_eff = np.array([1.77, 1.02, 0.48, 0.41])  # the measured trend's shape
    clock = fit_strain_clock(rates, gamma_eff)
    assert not clock.form_holds
    assert clock.t_r < 0


def test_strain_clock_prefers_strain_when_the_data_collapses_on_strain():
    rates = np.array([0.5, 1.0, 2.0, 5.0])
    clock = fit_strain_clock(rates, np.full(4, 1.2))
    assert clock.strain_spread < clock.time_spread
    np.testing.assert_allclose(clock.time_spread, 10.0, rtol=1e-12)


def test_strain_clock_needs_three_rates():
    with pytest.raises(ValueError, match="three rates"):
        fit_strain_clock(np.array([1.0, 2.0]), np.array([0.5, 0.6]))


def test_strain_clock_refuses_non_positive_inputs():
    with pytest.raises(ValueError, match="positive"):
        fit_strain_clock(np.array([1.0, 2.0, 3.0]), np.array([0.5, -0.6, 0.7]))


# --------------------------------------------------------------------------
# the reader
# --------------------------------------------------------------------------
def test_transient_panels_are_registered():
    assert set(LIN_TRANSIENTS) == {"Fig_2b", "Fig_2d", "Fig_3a"}
    assert LIN_TRANSIENTS["Fig_2b"][1] == "Fig_2a"
    assert LIN_TRANSIENTS["Fig_2d"][1] == "Fig_2c"
    assert LIN_TRANSIENTS["Fig_3a"][1] is None


def test_reader_refuses_an_unknown_panel():
    with pytest.raises(KeyError, match="unknown panel"):
        read_lin_transients(".", "Fig_9z")


@needs_lin
@pytest.mark.parametrize(
    "panel, expected",
    [
        ("Fig_2b", (0.5, 1.0, 2.0, 5.0)),
        ("Fig_2d", (0.05, 0.1, 0.2, 1.0)),
        ("Fig_3a", (0.02, 0.13, 0.2, 0.74)),
    ],
)
def test_reader_parses_every_rate_from_inconsistent_file_names(panel, expected):
    """0.5, 1.0, 0.05 and 0.20 are all formatted differently in the deposit."""
    t = read_lin_transients(LIN_ROOT, panel)
    assert t.rates == expected
    assert len(t.strain) == len(t.eta_r) == len(expected)


@needs_lin
def test_reader_returns_strictly_increasing_strain():
    t = read_lin_transients(LIN_ROOT, "Fig_2b")
    for g in t.strain:
        assert np.all(np.diff(g) > 0)


@needs_lin
def test_records_are_uniformly_sampled_in_strain():
    """The fit weights every sample equally, which is only fair if this holds."""
    t = read_lin_transients(LIN_ROOT, "Fig_2b")
    for g in t.strain:
        d = np.diff(g)
        np.testing.assert_allclose(d, np.median(d), rtol=1e-6)


# --------------------------------------------------------------------------
# the fit-free landmark
# --------------------------------------------------------------------------
def test_half_recovery_strain_matches_the_exponential_it_was_built_from():
    """For a clean exponential starting at f = 0, gamma_50 = ln(2) * gamma_eff."""
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.0, 0.80, phi, phi_m, phi_0)
    g50, plateau = half_recovery_strain(g, eta, phi, phi_m, phi_0)
    np.testing.assert_allclose(g50, np.log(2.0) * 0.35, rtol=0.05)
    np.testing.assert_allclose(plateau, 0.80, rtol=0.02)


def test_half_recovery_strain_needs_no_fit_window():
    """It is a landmark, so it survives a record the exponential fit would refuse."""
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.0, 0.80, phi, phi_m, phi_0, seed=11, noise=0.03)
    g50, _ = half_recovery_strain(g, eta, phi, phi_m, phi_0)
    assert np.isfinite(g50) and g50 > 0


def test_half_recovery_strain_is_invariant_along_the_degenerate_direction():
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0)
    base = half_recovery_strain(g, eta, phi, phi_m, phi_0)[0]
    for scale in (0.92, 1.08):
        moved = half_recovery_strain(g, eta, phi * scale, phi_m * scale, phi_0 * scale)[0]
        np.testing.assert_allclose(moved, base, rtol=1e-12)


def test_half_recovery_strain_is_measured_from_the_reversal_minimum():
    """Not from the start of the record, which is where the dip has not happened yet.

    Prepending a pre-reversal decay moves the record's origin but not its
    minimum, so the landmark must not move. This is the property that lets
    gamma_50 be compared across rates whose dips sit at different strains.
    """
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.0, 0.80, phi, phi_m, phi_0)
    base = half_recovery_strain(g, eta, phi, phi_m, phi_0)[0]

    step = float(np.median(np.diff(g)))
    lead = np.arange(1, 401) * step
    g_shift = np.concatenate([lead, lead[-1] + g])
    # A decaying lead-in that stays above the dip, so the minimum is unchanged.
    eta_shift = np.concatenate([np.linspace(9.0, eta[0], lead.size), eta])

    shifted = half_recovery_strain(g_shift, eta_shift, phi, phi_m, phi_0)[0]
    np.testing.assert_allclose(shifted, base, rtol=0.02)


@pytest.mark.parametrize("bad", [0.0, 1.0, -0.5, 1.5])
def test_half_recovery_strain_refuses_a_fraction_outside_the_unit_interval(bad):
    phi, phi_m, phi_0 = 0.40, 0.52, 0.64
    g, eta = synthetic_transient(0.35, 0.02, 0.80, phi, phi_m, phi_0)
    with pytest.raises(ValueError, match="fraction"):
        half_recovery_strain(g, eta, phi, phi_m, phi_0, fraction=bad)
