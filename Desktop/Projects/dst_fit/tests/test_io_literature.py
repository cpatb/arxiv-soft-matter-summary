"""Tests for the literature readers and the constant sets they carry.

The deposited data is not redistributed with the package, so the tests that need
it skip when it is absent. Everything that can be tested without it -- the
weight-to-volume conversion, the constant sets, the refusals -- is tested
unconditionally.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from dst_fit.io_literature import (
    CONSTANT_SETS,
    GUY_FIG1A_PHI,
    GUY_FIG1B_PHI,
    GUY_PMMA_404,
    GUY_PMMA_3770,
    RICHARDS_CORNSTARCH,
    log_crossing,
    mass_to_volume_fraction,
    peak_log_slope,
    read_guy_branch,
    read_guy_flow_curves,
    read_guy_phi_m_column,
    read_lin_reversal,
    read_richards_series,
    verify_guy_size_assignment,
)
from dst_fit.models import WYART_CATES

ROOT = Path(__file__).resolve().parents[1]
RICHARDS_DIR = (
    ROOT
    / "Competing Timescales Lead to Oscillations in Shear-Thickening Suspensions"
    / "Figure1"
    / "Plotting"
)
FIG1 = ROOT / "Towards a unified description of the rheology of hard-particle suspensions" / "Fig 1" / "Fig 1"
FIG1A = FIG1 / "Fig 1a" / "Fig1a.xls"
FIG1B = FIG1 / "Fig 1b" / "Fig1b.xlsx"
FIG1C = FIG1 / "Fig 1c" / "Fig1c.xlsx"

needs_richards = pytest.mark.skipif(
    not RICHARDS_DIR.is_dir(), reason="Richards et al. deposit not present"
)
needs_guy = pytest.mark.skipif(
    not FIG1C.is_file(), reason="Guy et al. deposit not present"
)


# --------------------------------------------------------------------------
# constant sets
# --------------------------------------------------------------------------
def test_constant_sets_are_the_published_numbers():
    """Pinned so a later edit to a constant has to be deliberate."""
    r = RICHARDS_CORNSTARCH
    assert (r.eta_s, r.phi_0, r.n, r.beta) == (6.0e-3, 0.546, 2.0, 0.94)
    assert r.phi_m_published == 0.457
    assert r.sigma_star_published == 5.1
    assert r.phi_units == "weight"

    assert GUY_PMMA_404.eta_s == 2.4e-3
    assert GUY_PMMA_3770.eta_s == 2.83e-3, "the two sizes have different solvents"
    for g in (GUY_PMMA_404, GUY_PMMA_3770):
        assert g.n == 2.2
        assert g.phi_m_published == 0.558 and g.phi_m_published_err == 0.005
        assert g.phi_units == "volume"
    assert GUY_PMMA_404.sigma_star_published == 200.0
    assert GUY_PMMA_3770.sigma_star_published == 2.0


def test_constant_sets_bind_straight_onto_the_model():
    m = WYART_CATES.with_fixed(**RICHARDS_CORNSTARCH.as_fixed())
    assert m.fixed == {"eta_s": 6.0e-3, "phi_0": 0.546, "n": 2.0, "beta": 0.94}
    assert m.bounds[1][1] == pytest.approx(0.545), "phi_m must be capped below phi_0"
    assert set(CONSTANT_SETS) == {c.name for c in CONSTANT_SETS.values()}


# --------------------------------------------------------------------------
# the conversion that is deliberately not applied
# --------------------------------------------------------------------------
def test_mass_to_volume_endpoints_and_monotonicity():
    w = np.linspace(0.0, 1.0, 21)
    phi = mass_to_volume_fraction(w, 1.55, 1.126)
    assert phi[0] == pytest.approx(0.0)
    assert phi[-1] == pytest.approx(1.0)
    assert np.all(np.diff(phi) > 0)


def test_mass_to_volume_is_the_identity_when_densities_match():
    w = np.array([0.1, 0.38, 0.5, 0.9])
    np.testing.assert_allclose(mass_to_volume_fraction(w, 1.3, 1.3), w, rtol=1e-12)


def test_mass_to_volume_is_not_a_constant_rescaling():
    """Which is the reason it cannot be folded into phi_0 and forgotten."""
    w = np.array([0.38, 0.44, 0.50])
    ratio = mass_to_volume_fraction(w, 1.55, 1.126) / w
    assert ratio.max() - ratio.min() > 0.02


def test_mass_to_volume_refuses_nonsense():
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        mass_to_volume_fraction(np.array([1.2]), 1.55, 1.126)
    with pytest.raises(ValueError, match="positive"):
        mass_to_volume_fraction(np.array([0.4]), -1.0, 1.126)


# --------------------------------------------------------------------------
# Richards
# --------------------------------------------------------------------------
@needs_richards
def test_richards_series_shape_and_units():
    sweeps = read_richards_series(RICHARDS_DIR)
    assert [s.phi for s in sweeps] == [0.38, 0.40, 0.42, 0.44, 0.46, 0.48, 0.50]
    for s in sweeps:
        assert s.metadata["phi_units"] == "weight", "these are not volume fractions"
        assert np.all(np.diff(s.stress) > 0), "must be ordered on the stress axis"
        assert np.all(s.viscosity > 0) and np.all(np.isfinite(s.viscosity))
        # viscosity is the derived quantity; the identity must hold exactly
        np.testing.assert_allclose(s.viscosity * s.shear_rate, s.stress, rtol=1e-12)
        assert s.metadata["shear_rate_error"].shape == s.shear_rate.shape


@needs_richards
def test_richards_stress_consistency_is_reported_as_untestable():
    """stress_reported is None on purpose: a self-comparison would fake a pass."""
    from dst_fit import check_stress_consistency

    s = read_richards_series(RICHARDS_DIR)[0]
    out = check_stress_consistency(s)
    assert out["ok"] is False and out["reason"] == "no stress column"


@needs_richards
def test_richards_curves_thicken():
    """If they did not, fitting Wyart-Cates to them would be meaningless."""
    for s in read_richards_series(RICHARDS_DIR):
        i = int(np.argmin(s.viscosity))
        assert s.viscosity[-1] / s.viscosity[i] > 2.0


# --------------------------------------------------------------------------
# Guy
# --------------------------------------------------------------------------
@needs_guy
def test_guy_phi_m_readable_without_fitting():
    """phi + (phi_m - phi) is constant across the deposited rows, and equals 0.560."""
    assert read_guy_phi_m_column(FIG1C) == pytest.approx(0.560, abs=1e-3)


@needs_guy
def test_guy_branch_sheets_are_ordered_and_positive():
    for sheet, n in [("d=3770 nm eta_1", 15), ("d=3770 nm eta_2", 6),
                     ("d=404 nm eta_1", 6)]:
        phi, eta_r = read_guy_branch(FIG1C, sheet)
        assert phi.size == n
        assert np.all(np.diff(phi) > 0)
        assert np.all(eta_r >= 1.0), "a suspension cannot be thinner than its solvent"


@needs_guy
def test_guy_branch_rejects_an_unknown_sheet():
    with pytest.raises(KeyError, match="no sheet"):
        read_guy_branch(FIG1C, "d=1 nm eta_9")


@needs_guy
def test_guy_flow_curves_are_stress_ordered_and_deduplicated():
    for path, phi_map in ((FIG1A, GUY_FIG1A_PHI), (FIG1B, GUY_FIG1B_PHI)):
        for s in read_guy_flow_curves(path, phi_map):
            assert np.all(np.diff(s.stress) >= 0), f"{s.sheet}: not ordered by stress"
            assert np.all(s.viscosity > 0)
            assert s.metadata["phi_label"] == s.sheet
            # exact duplicate rows must be gone, but repeated stresses carrying
            # different viscosities are real DST data and must survive
            rows = np.column_stack([s.stress, s.viscosity])
            assert np.unique(rows, axis=0).shape[0] == rows.shape[0]


@needs_guy
def test_repeated_stresses_are_counted_not_silently_dropped():
    """One imposed stress, several shear rates: the S-curve, which no eta(sigma) fits."""
    by_label = {s.sheet: s for s in read_guy_flow_curves(FIG1A, GUY_FIG1A_PHI)}
    assert by_label["59%"].metadata["n_repeated_stresses"] > 0
    assert by_label["40%"].metadata["n_repeated_stresses"] == 0


@needs_guy
def test_guy_unmapped_sheet_gets_no_phi_rather_than_a_guess():
    """59% and 61% have no entry in the authors' phi table, so phi must be None."""
    by_label = {s.sheet: s for s in read_guy_flow_curves(FIG1A, GUY_FIG1A_PHI)}
    assert by_label["59%"].phi is None
    assert by_label["61%"].phi is None
    assert by_label["40%"].phi == pytest.approx(0.39296)


@needs_guy
def test_size_assignment_recovers_the_stated_solvent_viscosities():
    """The check that fixes which workbook is which size, from the numbers.

    Its output is a measurement, not a boolean: the ratio of the measured
    lowest-stress viscosity to the tabulated relative viscosity is the solvent
    viscosity, and it has to come back at the value the paper states.
    """
    a = verify_guy_size_assignment(FIG1A, FIG1C, "d=404 nm eta_1", GUY_FIG1A_PHI,
                                   eta_s=GUY_PMMA_404.eta_s)
    b = verify_guy_size_assignment(FIG1B, FIG1C, "d=3770 nm eta_1", GUY_FIG1B_PHI,
                                   eta_s=GUY_PMMA_3770.eta_s)
    assert a.matched and b.matched
    assert np.median(a.implied_eta_s) == pytest.approx(2.4e-3, rel=0.03)
    assert np.median(b.implied_eta_s) == pytest.approx(2.83e-3, rel=0.05)


@needs_guy
def test_the_wrong_solvent_viscosity_is_caught_not_absorbed():
    """2.4 mPa.s on the large particles is the error DATA_INVENTORY originally had."""
    wrong = verify_guy_size_assignment(FIG1B, FIG1C, "d=3770 nm eta_1", GUY_FIG1B_PHI,
                                       eta_s=GUY_PMMA_404.eta_s)
    assert not wrong.matched
    assert min(wrong.ratio) > 1.0, "and it fails high, by roughly 2.83/2.4"


@needs_guy
def test_size_assignment_refuses_a_phi_not_in_the_branch_sheet():
    with pytest.raises(ValueError, match="no phi within"):
        verify_guy_size_assignment(FIG1A, FIG1C, "d=3770 nm eta_1", GUY_FIG1A_PHI,
                                   eta_s=GUY_PMMA_404.eta_s)


# --------------------------------------------------------------------------
# Lin et al.: the reversal split, and the landmarks compared against it
# --------------------------------------------------------------------------
LIN_DIR = (
    ROOT
    / "Hydrodynamic and Contact Contributions to Continuous Shear Thickening in Colloidal Suspensions"
    / "Fig2"
    / "Fig2"
)
needs_lin = pytest.mark.skipif(
    not LIN_DIR.is_dir(), reason="Lin et al. deposit not present"
)


@needs_lin
@pytest.mark.parametrize("system", ["Fig_2a", "Fig_2c"])
def test_lin_decomposition_closes_exactly(system):
    """contact + hydro = total. If this drifts, the columns are misread."""
    split = read_lin_reversal(LIN_DIR, system)
    assert split.closure_error().max() < 1e-6


@needs_lin
@pytest.mark.parametrize("system", ["Fig_2a", "Fig_2c"])
def test_lin_reduced_stress_needs_no_solvent_viscosity(system):
    """Sigma = eta_r * gamma_dot is sigma / eta_s, and both sides are ordered."""
    split = read_lin_reversal(LIN_DIR, system)
    np.testing.assert_allclose(split.sweep_sigma, split.sweep_eta_r * split.sweep_rate)
    np.testing.assert_allclose(split.rev_sigma, split.rev_total * split.rev_rate)
    assert np.all(np.diff(split.sweep_rate) > 0)
    assert np.all(np.diff(split.rev_rate) > 0)


@needs_lin
@pytest.mark.parametrize("system", ["Fig_2a", "Fig_2c"])
def test_lin_contact_fraction_is_a_fraction_and_rises(system):
    split = read_lin_reversal(LIN_DIR, system)
    chi = split.contact_fraction
    assert np.all((chi > 0) & (chi < 1))
    assert chi[-1] > chi[0], "contacts must take over as stress rises"
    assert chi.max() < 0.8, "hydrodynamic stress never switches off entirely"


@needs_lin
def test_lin_refuses_an_unknown_system():
    with pytest.raises(KeyError, match="unknown system"):
        read_lin_reversal(LIN_DIR, "Fig_9z")


def test_log_crossing_interpolates_in_log_not_linear():
    """The abscissa spans decades, so a linear interpolation lands in the wrong place."""
    x = np.array([1.0, 100.0])
    y = np.array([0.0, 1.0])
    got = log_crossing(x, y, 0.5)
    assert got == pytest.approx(10.0), "the log-midpoint of 1 and 100 is 10, not 50.5"


def test_log_crossing_refuses_to_extrapolate():
    x = np.array([1.0, 2.0, 4.0])
    y = np.array([0.1, 0.2, 0.3])
    assert np.isnan(log_crossing(x, y, 0.9))
    assert np.isnan(log_crossing(x, y, 0.05))


def test_log_crossing_finds_an_interior_target():
    x = np.array([1.0, 10.0, 100.0])
    y = np.array([0.2, 0.4, 0.9])
    got = log_crossing(x, y, 0.4)
    assert got == pytest.approx(10.0)


def test_peak_log_slope_finds_a_planted_power_law_step():
    """A curve with one steep decade must report that decade, and its slope."""
    x = np.logspace(0, 4, 200)
    # flat, then slope 1 for a decade around x = 100, then flat again
    ly = np.clip(np.log10(x) - 1.5, 0.0, 1.0)
    where, slope = peak_log_slope(x, 10.0**ly)
    assert 30 < where < 320
    assert slope == pytest.approx(1.0, rel=0.05)


def test_peak_log_slope_needs_three_points():
    with pytest.raises(ValueError, match="three points"):
        peak_log_slope(np.array([1.0, 2.0]), np.array([1.0, 2.0]))
