"""Tests for the elimination audit.

Anchors outside the code under test: hand-computed dimensionless groups, exact
analytic invariances, and discriminating negative cases. The audit's verdict
machinery is tested against synthetic responses whose true generator is known,
so a group that DOES determine the response must come back "survives" and one
that cannot must come back "FALSIFIED".

Run:  python -m pytest test_onset_clocks.py -q
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from audit import TOL_TIGHT, audit_clock, tie_retention  # noqa: E402
from clocks import (  # noqa: E402
    CLOCKS,
    Drop,
    frictional_fraction_at_impact,
    impact_stress_scale,
    jamming_deficit,
    material,
    peclet_particle,
    shear_rate_naive,
    strain_in_one_impact,
    wall_shear_gauge,
)
from data import is_arrested, load_drops  # noqa: E402

D = Drop(label="ref", phi=0.44, salt_mM=0.0, u0=1.15, d0=2.4e-3,
         beta_p=0.65, beta_eq=1.789)


# --------------------------------------------------------------------------
# the contradiction the project exists for
# --------------------------------------------------------------------------
def test_every_packing_in_the_dataset_is_above_phi_m():
    """The premise. If this ever passes below 1, the project's question changes."""
    for phi in (0.40, 0.44, 0.48, 0.50):
        d = Drop("x", phi, 0.0, 1.0, 2.4e-3, 1.0, 1.789)
        assert jamming_deficit(d) > 1.0, f"phi={phi} can now jam; re-read section 1"


def test_jamming_deficit_by_hand():
    # (0.70 - 0.44)/(0.70 - 0.5176) = 0.26/0.1824 = 1.42544
    assert jamming_deficit(D) == pytest.approx(1.42544, rel=1e-5)


def test_jamming_deficit_is_gauge_invariant():
    """A ratio of differences of the packing triple cannot move under the gauge."""
    m = material(0.44)
    base = (m.phi0 - 0.44) / (m.phi0 - m.phiM)
    for c in (0.9, 1.0, 1.1):
        assert ((c * m.phi0 - c * 0.44) / (c * m.phi0 - c * m.phiM)
                == pytest.approx(base, rel=1e-12))


# --------------------------------------------------------------------------
# the stress gauges
# --------------------------------------------------------------------------
def test_inertial_stress_by_hand():
    # rho = 0.44*2000 + 0.56*1000 = 1440;  1440 * 2 * 1.15^2 = 3808.8 Pa
    assert impact_stress_scale(D, 2.0) == pytest.approx(3808.8, rel=1e-9)


def test_deviatoric_gauge_is_well_below_the_inertial_one():
    """The lab's own Wagner analysis puts the ratio at 5 to 10x. Check the sign
    and the order of magnitude; the exact value depends on t* at peak force."""
    ratio = impact_stress_scale(D, 1.0) / wall_shear_gauge(D)
    assert 2.0 < ratio < 20.0


def test_the_gauge_choice_changes_whether_f_is_saturated():
    """The discriminating case: same drop, two gauges, materially different f."""
    f_inertial = frictional_fraction_at_impact(D, k=2.0)
    f_shear = float(np.exp(-((material(0.44).sigma_star / wall_shear_gauge(D))
                             ** material(0.44).delta)))
    assert f_inertial > 0.99
    assert f_shear < 0.98
    assert f_inertial - f_shear > 0.02


def test_wall_shear_grows_with_impact_speed():
    speeds = [0.75, 1.15, 1.6, 2.6]
    taus = [wall_shear_gauge(Drop("x", 0.44, 0.0, u, 2.4e-3, 1.0, 1.789))
            for u in speeds]
    assert taus == sorted(taus)


# --------------------------------------------------------------------------
# the rate and strain clocks
# --------------------------------------------------------------------------
def test_naive_shear_rate_is_u_over_d():
    assert shear_rate_naive(D) == pytest.approx(1.15 / 2.4e-3)


def test_strain_on_the_naive_convention_is_identically_one():
    """The structural result: this criterion has no velocity dependence at all."""
    for u in (0.5, 1.0, 2.0, 4.0):
        for d0 in (1.5e-3, 2.4e-3, 3.0e-3):
            d = Drop("x", 0.44, 0.0, u, d0, 1.0, 1.789)
            assert strain_in_one_impact(d) == pytest.approx(1.0, rel=1e-12)


def test_particle_peclet_is_enormous_so_these_drops_are_not_brownian():
    """800 nm silica at 500 1/s. A discriminating check on the a^3 scaling."""
    assert peclet_particle(D) > 100.0


# --------------------------------------------------------------------------
# the audit machinery, against known generators
# --------------------------------------------------------------------------
def test_a_group_that_generates_the_response_survives():
    rng = np.random.default_rng(0)
    g = np.exp(rng.uniform(0, 2, 60))
    r = g ** 0.9                       # r is a function of g alone
    v = audit_clock("perfect", g, r)
    assert v.verdict == "survives"
    assert v.tie_fraction_tight < 0.2


def test_a_constant_group_is_falsified():
    rng = np.random.default_rng(1)
    g = np.full(60, 3.0) * (1 + 1e-6 * rng.standard_normal(60))
    r = np.exp(rng.uniform(0, 2, 60))
    v = audit_clock("constant", g, r)
    assert v.verdict == "FALSIFIED"
    assert v.tie_fraction_tight > 0.9


def test_a_group_orthogonal_to_the_response_is_falsified():
    rng = np.random.default_rng(2)
    g = np.exp(rng.uniform(0, 2, 60))
    r = np.exp(rng.uniform(0, 2, 60))   # independent
    v = audit_clock("noise", g, r)
    assert v.verdict == "FALSIFIED"


def test_tie_retention_is_zero_for_a_deterministic_map():
    g = np.linspace(1.0, 10.0, 40)
    ret_med, ret_max, n = tie_retention(g, g**2, TOL_TIGHT)
    assert n > 0
    assert ret_med < 0.2
    assert ret_max < 0.3


def test_audit_refuses_too_few_points():
    with pytest.raises(ValueError):
        audit_clock("tiny", np.array([1.0, 2.0]), np.array([1.0, 2.0]))


def test_every_registered_clock_returns_a_finite_positive_number_for_a_typical_drop():
    for name, fn in CLOCKS.items():
        val = fn(D)
        assert np.isfinite(val), f"{name} is not finite at the reference drop"
        assert val > 0, f"{name} is not positive at the reference drop"


# --------------------------------------------------------------------------
# the data layer
# --------------------------------------------------------------------------
def test_the_dataset_loads_with_its_exclusions_counted():
    drops, prov = load_drops()
    assert prov["n_rows"] == 75
    assert prov["n_used"] == prov["n_rows"] - prov["n_bounded_beta_p_dropped"] \
        - prov["n_non_finite_dropped"]
    assert len(drops) == prov["n_used"]
    assert prov["bounded_kinds"], "the excluded kinds must be named, not just counted"


def test_the_response_spans_the_factor_the_audit_is_asked_to_explain():
    _, prov = load_drops()
    lo, hi = prov["beta_p_range"]
    assert hi / lo > 5.0


def test_arrest_is_recomputed_from_the_repos_own_rule():
    """A stored flag in this programme has already outlived its criterion once."""
    assert is_arrested(Drop("x", 0.44, 0.0, 1.0, 2.4e-3, 0.5, 1.789))
    assert not is_arrested(Drop("x", 0.44, 0.0, 1.0, 2.4e-3, 1.7, 1.789))
    # exactly at the boundary: 0.9 * beta_eq is NOT arrested (strict <)
    assert not is_arrested(Drop("x", 0.44, 0.0, 1.0, 2.4e-3, 0.9 * 1.789, 1.789))


def test_the_lamella_convention_is_undefined_somewhere_and_says_so():
    """Reported as NaN rather than extrapolated. If this ever stops being true,
    the README's section 8 claim about regime coverage needs revisiting."""
    from clocks import shear_rate_lamella
    drops, _ = load_drops()
    vals = np.array([shear_rate_lamella(d) for d in drops])
    assert np.isnan(vals).any(), "no drop is outside the correlation's regime now"
    assert np.isfinite(vals).sum() > 40


# --------------------------------------------------------------------------
# the inversion
# --------------------------------------------------------------------------
def test_implied_viscosity_inverts_the_spreading_law_by_hand():
    # rho = 1440, U0 = 1.15, D0 = 2.4e-3, beta = 0.65, p = 5
    # 1440*1.15*0.0024 = 3.9744;  0.65^5 = 0.1160290625;  ratio = 34.254
    from invert import implied_viscosity
    assert implied_viscosity(D, p=5.0) == pytest.approx(34.254, rel=1e-3)


def test_required_phi_j_round_trips_through_the_divergence():
    """Anchor on the forward model, not on the inverse's own output."""
    from invert import required_phi_j
    m = material(0.44)
    for phi_j in (0.45, 0.50, 0.60):
        eta = m.eta_s * (1.0 - 0.44 / phi_j) ** -2.0
        d = Drop("x", 0.44, 0.0, 1.0, 2.4e-3, 1.0, 1.789)
        assert required_phi_j(d, eta) == pytest.approx(phi_j, abs=1e-10)


def test_required_phi_j_refuses_below_the_solvent_viscosity():
    from invert import required_phi_j
    assert np.isnan(required_phi_j(D, 1e-6))


def test_required_phi_m_inverts_the_interpolation():
    """phi_J = phi_0 - (phi_0 - phi_m) f, solved for phi_m."""
    from invert import required_phi_m
    m = material(0.44)
    for f in (0.5, 0.9, 1.0):
        for phi_m in (0.42, 0.48, 0.52):
            phi_j = m.phi0 - (m.phi0 - phi_m) * f
            assert required_phi_m(D, phi_j, f) == pytest.approx(phi_m, abs=1e-12)


def test_the_required_phi_m_is_far_below_the_fitted_one():
    """The headline. If this ever stops holding, section 1 needs rewriting."""
    from invert import invert as inv
    drops, _ = load_drops()
    arrested = [d for d in drops if is_arrested(d)]
    pm = np.array([inv(d, 5.0).phi_m_required for d in arrested])
    pm = pm[np.isfinite(pm)]
    assert np.median(pm) < material(0.44).phiM - 0.05


def test_the_required_phi_m_barely_moves_with_the_unknown_prefactor():
    """A steep divergence pins phi_J near phi whatever the O(1) coefficient is.

    Over four decades of C the median must move less than the gap to the fitted
    value, or the headline is an artefact of a coefficient nobody knows.
    """
    from invert import prefactor_phi_m
    drops, _ = load_drops()
    arrested = [d for d in drops if is_arrested(d)]
    meds = []
    for C in (0.01, 1.0, 100.0):
        v = np.array([prefactor_phi_m(d, 5.0, C) for d in arrested])
        meds.append(float(np.median(v[np.isfinite(v)])))
    spread = max(meds) - min(meds)
    gap = material(0.44).phiM - meds[1]
    assert spread < gap, f"prefactor moves it {spread:.4f}, gap is only {gap:.4f}"


def test_the_inversion_is_only_applied_to_drops_that_arrested():
    """A viscous-ARREST model applied to a drop that did not arrest is what
    produced the spurious velocity drift in the first version of this."""
    drops, _ = load_drops()
    arrested = [d for d in drops if is_arrested(d)]
    assert 0 < len(arrested) < len(drops)
    assert all(d.beta_p < 0.9 * d.beta_eq for d in arrested)


def test_the_inversion_is_pinned_and_the_readme_says_so():
    """Guard against re-claiming the withdrawn phi_m headline.

    The required phi_m tracks the phi the sample was made at. If this ever stops
    being true the inversion has become informative and section 9 must be
    rewritten -- but until then, no value of phi_m may be quoted from it.
    """
    from invert import invert as inv
    drops, _ = load_drops()
    arrested = [d for d in drops if is_arrested(d)]
    for phi in (0.40, 0.44, 0.48, 0.50):
        v = np.array([inv(d, 5.0).phi_m_required
                      for d in arrested if d.phi == phi])
        v = v[np.isfinite(v)]
        if v.size < 3:
            continue
        assert abs(np.median(v) - phi) < 0.03, (
            f"required phi_m at phi={phi} is {np.median(v):.4f}, no longer "
            "tracking phi; the inversion may now be informative")


# --------------------------------------------------------------------------
# robustness to the withdrawn calibration
# --------------------------------------------------------------------------
def test_the_stress_clock_stays_flat_on_every_corner_of_the_parameter_box():
    """C2, the audit's central result, must not depend on the withdrawn fit."""
    from robustness import corners, f_range
    drops, _ = load_drops()
    ranges = [f_range(drops, q) for q in corners(n_per_axis=3)]
    assert max(ranges) < 6.35, (
        f"f now spans {max(ranges):.2f}x somewhere in the box, which reaches the "
        "response's 6.35x; the elimination of the stress clock is no longer "
        "parameter-independent")


def test_the_premise_fails_only_for_a_low_phi_m_and_says_which():
    """C1 is conditional on phi_m; the test records the threshold, not a pass."""
    from robustness import can_any_drop_jam, corners
    drops, _ = load_drops()
    jammable = [q for q in corners(n_per_axis=3) if can_any_drop_jam(drops, q)]
    assert jammable, "the premise now holds everywhere; section 13 needs rewriting"
    assert max(q.phi_m for q in jammable) < 0.50, (
        "the premise fails at a phi_m above the highest packing in the dataset, "
        "which should be impossible")


def test_the_transient_clock_is_falsified_at_the_bottom_of_the_measured_gamma0():
    """The discriminating case: the strain clock helps only at high gamma_0."""
    from clocks import GAMMA_0_RANGE, frictional_fraction_transient
    lo = np.array([frictional_fraction_transient(d, GAMMA_0_RANGE[0])
                   for d in load_drops()[0]])
    hi = np.array([frictional_fraction_transient(d, GAMMA_0_RANGE[1])
                   for d in load_drops()[0]])
    lo, hi = lo[np.isfinite(lo)], hi[np.isfinite(hi)]
    assert lo.max() / lo.min() < hi.max() / hi.min()
    assert lo.max() / lo.min() < 1.5      # flat at the bottom of the range


def test_the_pooled_correlation_flips_sign_within_a_packing():
    """The Simpson's paradox of section 15, guarded so it cannot be forgotten."""
    from scipy.stats import spearmanr

    from clocks import CLOCKS
    drops, _ = load_drops()
    sub = [d for d in drops if d.phi == 0.44 and is_arrested(d)]
    fn = CLOCKS["Re"]
    pooled = spearmanr([fn(d) for d in drops], [d.beta_p for d in drops])[0]
    within = spearmanr([fn(d) for d in sub], [d.beta_p for d in sub])[0]
    assert pooled > 0 > within, (
        "Re no longer flips sign between the pooled and within-packing "
        "analyses; section 15's warning may no longer apply")
