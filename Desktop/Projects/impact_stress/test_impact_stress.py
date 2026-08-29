"""
test_impact_stress.py — validation suite.

Philosophy (lab CLAUDE.md): treat every function as guilty until verified.
Four kinds of test here, in increasing order of how much they prove:

  1. LIMITS      — does the formula reduce to the known answer in a limit?
  2. IDENTITY    — do two algebraically-equivalent forms agree numerically?
                   (catches transcription errors in the rewriting)
  3. LITERATURE  — does an independently-coded printed equation agree with ours?
  4. CROSS-CHECK — does a quantity computed two physically different ways agree?
                   (the only kind that can catch a wrong derivation)

Run:  python -m pytest test_impact_stress.py -q
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

import dilatancy as dl
import materials as mat
import stress as st

# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------
SUSP = mat.lab_silica()
DROP = mat.suspension_as_fluid(SUSP)
R_LAB = 1.25e-3      # 2.5 mm drop diameter, mid-range of the lab's 2.3-2.9 mm
U_LAB = 2.3          # band B centre


def make(U=U_LAB, R=R_LAB, target=None, drop=None):
    return st.Impact(drop=drop or DROP, U=U, R=R,
                     target=target or st.rigid_target())


# ==========================================================================
# 1. LIMITS
# ==========================================================================
def test_water_hammer_rigid_limit_is_rho_c_U():
    imp = make(target=st.rigid_target())
    assert_allclose(st.water_hammer_pressure(imp), DROP.rho * DROP.c * U_LAB, rtol=1e-12)


def test_water_hammer_matched_target_is_exactly_half():
    """A drop hitting a deep bath of ITSELF sees exactly half the rigid-wall
    water-hammer pressure. This is the sharpest single statement about what a
    fluid substrate does in the compressible stage."""
    imp = make(target=st.liquid_target(DROP))
    assert_allclose(st.water_hammer_pressure(imp),
                    0.5 * DROP.rho * DROP.c * U_LAB, rtol=1e-12)


def test_water_hammer_free_surface_is_zero():
    free = st.Target(name="vacuum", rho=0.0, Z=0.0, E=0.0)
    assert st.water_hammer_pressure(make(target=free)) == 0.0


def test_compliance_factors_rigid_and_matched():
    assert_allclose(st.compliance_factors(make())["acoustic_factor"], 1.0)
    assert_allclose(st.compliance_factors(make())["inertial_factor"], 1.0)
    m = st.compliance_factors(make(target=st.liquid_target(DROP)))
    assert_allclose(m["acoustic_factor"], 0.5, rtol=1e-12)
    assert_allclose(m["inertial_factor"], 0.5, rtol=1e-12)


def test_wyart_cates_limits():
    s = SUSP
    assert_allclose(s.f_ss(np.array([1e-9])), 0.0, atol=1e-12)          # no stress, no friction
    assert_allclose(s.f_ss(np.array([1e12])), 1.0, rtol=1e-6)           # huge stress, all frictional
    assert s.eta_frictional > s.eta_0                                    # thickening, not thinning
    assert_allclose(s.phi_J(0.0), s.phi0)
    assert_allclose(s.phi_J(1.0), s.phiM)


def test_dilute_limit_recovers_solvent_viscosity():
    dilute = mat.lab_silica(phi=1e-9)
    assert_allclose(dilute.eta_0, dilute.eta_s, rtol=1e-6)


# ==========================================================================
# 2. IDENTITY — algebraically equivalent forms must agree numerically
# ==========================================================================
def test_wagner_pressure_dimensional_and_selfsimilar_agree():
    imp = make()
    t_star = np.array([0.01, 0.1, 0.5, 1.0])
    xi = np.array([0.0, 0.3, 0.7, 0.95])
    for ts in t_star:
        t = ts * imp.tau
        d = st.contact_radius(imp, t)
        for x in xi:
            assert_allclose(st.wagner_pressure(imp, t, x * d),
                            st.wagner_pressure_selfsimilar(imp, ts, x), rtol=1e-10)


def test_wagner_axis_matches_selfsimilar_at_xi_zero():
    imp = make()
    ts = np.array([0.01, 0.3, 1.0, 3.0])
    assert_allclose(st.wagner_pressure_axis(imp, ts),
                    st.wagner_pressure_selfsimilar(imp, ts, 0.0), rtol=1e-12)


def test_shear_dimensional_and_selfsimilar_agree():
    imp = make()
    for ts in [0.01, 0.1, 1.0]:
        t = ts * imp.tau
        d = st.contact_radius(imp, t)
        for x in [0.1, 0.5, 0.9, 0.98]:
            assert_allclose(st.wall_shear_stress(imp, t, x * d),
                            st.wall_shear_selfsimilar(imp, ts, x), rtol=1e-10)


def test_onset_velocity_closed_form_matches_rootfind():
    U_root = dl.onset_velocity(SUSP, R_LAB, criterion="mean", t_star=1.0)
    U_closed = dl.onset_velocity_closed_form(SUSP.sigma_star, R_LAB, DROP.rho,
                                             DROP.eta, t_star=1.0)
    assert_allclose(U_root, U_closed, rtol=1e-6)


# ==========================================================================
# 3. LITERATURE — independently coded printed equations
# ==========================================================================
def test_contact_radius_matches_CSG22_printed_form():
    """CSG22 Fig. 2a:  r0 = sqrt(6 D V t)/2.  Ours: d = sqrt(3 R U t)."""
    imp = make()
    t = np.logspace(-7, -3, 25)
    r0_paper = np.sqrt(6.0 * imp.D * imp.U * t) / 2.0
    assert_allclose(st.contact_radius(imp, t), r0_paper, rtol=1e-12)


def test_wagner_pressure_matches_CSG22_eq10_printed_form():
    """CSG22 Eq. [10] as printed:
        p = (3 sqrt2 / 2 pi) rho V^2 D / sqrt(3 V D t - 2 r^2)
    coded here directly from the paper, with no algebra in common with ours."""
    imp = make()
    rho, V, D = imp.drop.rho, imp.U, imp.D
    for ts in [0.005, 0.05, 0.5]:
        t = ts * imp.tau
        d = st.contact_radius(imp, t)
        r = np.array([0.0, 0.25, 0.5, 0.75, 0.9]) * d
        p_paper = (3.0 * np.sqrt(2.0) / (2.0 * np.pi)) * rho * V**2 * D / np.sqrt(
            3.0 * V * D * t - 2.0 * r**2)
        assert_allclose(st.wagner_pressure(imp, t, r), p_paper, rtol=1e-11)


def test_wall_shear_matches_CSG22_eq13_printed_form():
    """CSG22 Eq. [13] as printed:
        tau = [sqrt6 rho^1/2 V^3/2 D^1/2 eta^1/2 / pi^3/2] * 2r/(3 V D t - 2 r^2)
    """
    imp = make()
    rho, V, D, eta = imp.drop.rho, imp.U, imp.D, imp.drop.eta
    for ts in [0.005, 0.05, 0.5]:
        t = ts * imp.tau
        d = st.contact_radius(imp, t)
        r = np.array([0.1, 0.3, 0.6, 0.9]) * d
        pref = np.sqrt(6.0) * np.sqrt(rho) * V**1.5 * np.sqrt(D) * np.sqrt(eta) / np.pi**1.5
        tau_paper = pref * 2.0 * r / (3.0 * V * D * t - 2.0 * r**2)
        assert_allclose(st.wall_shear_stress(imp, t, r), tau_paper, rtol=1e-11)


def test_selfsimilar_amplitude_prefactors():
    """The two amplitudes quoted in the report must be what the code produces."""
    imp = make()
    # p_axis = (sqrt3/pi) rho U^2 t*^(-1/2) = 0.5513 rho U^2 t*^(-1/2)
    assert_allclose(float(st.wagner_pressure_axis(imp, 1.0)) / imp.p_dyn,
                    np.sqrt(3.0) / np.pi, rtol=1e-12)
    # <tau> = 1.1353 rho U^2 Re_R^(-1/2) t*^(-1/2)
    assert_allclose(float(st.wall_shear_mean(imp, 1.0)) / (imp.p_dyn * imp.Re_R**-0.5),
                    10.7 / (3.0 * np.pi), rtol=1e-12)
    # pressure rim amplification = 1/sqrt(2 * DELTA_OVER_D) (to O(Delta))
    ratio = float(st.wagner_pressure_peak(imp, 1.0)) / float(st.wagner_pressure_axis(imp, 1.0))
    assert_allclose(ratio, 1.0 / np.sqrt(1.0 - (1 - st.DELTA_OVER_D) ** 2), rtol=1e-12)
    assert 5.0 < ratio < 6.0


def _repo_closure_module():
    import importlib.util
    import sys
    from pathlib import Path

    p = Path(r"C:\Users\compl\Github\Drop Impact\closure\closure.py")
    if not p.exists():
        pytest.skip("drop-impact repo not present")
    spec = importlib.util.spec_from_file_location("repo_closure", p)
    m = importlib.util.module_from_spec(spec)
    sys.modules["repo_closure"] = m
    spec.loader.exec_module(m)
    return m


def test_wyart_cates_matches_repo_closure_implementation():
    """CROSS-IMPLEMENTATION: our WyartCates must agree with the drop-impact
    repo's own SuspensionClosure — written independently, and the code the
    calibration was fitted with.

    Compared at phi = 0.55, NOT at the paper's phi = 0.44, because the repo's
    class refuses to instantiate below phiM (see the next test)."""
    m = _repo_closure_module()
    s = mat.lab_silica(phi=0.55)
    repo = m.SuspensionClosure(
        eta_s=s.eta_s, phi0=s.phi0, phiM=s.phiM, sigma_star=s.sigma_star,
        delta=s.delta, phi=s.phi, tau_relax=2.0, name="silica", placeholder=False)

    sig = np.logspace(-1, 4, 40)
    assert_allclose(s.f_ss(sig), repo.f_ss(sig), rtol=1e-12)
    f = np.linspace(0.0, 1.0, 21)
    assert_allclose(s.eta_of_f(f), repo.eta(f), rtol=1e-12)
    assert_allclose(s.eta_0, repo.eta0, rtol=1e-12)


def test_repo_closure_rejects_its_own_fitted_silica_at_the_papers_phi():
    """FINDING, pinned as a test so it cannot be lost.

    closure.py asserts phiM < phi < phi0. The fitted silica set has
    phiM = 0.5176 and the paper's operating point is phi = 0.44, so the repo's
    own closure CANNOT be constructed at the volume fraction the experiments
    are run at. Physically the same statement: at phi = 0.44 < phiM this
    material cannot shear-jam at any stress in the Wyart-Cates picture, only
    thicken continuously. The drops nevertheless arrest.

    That is not a bug in either code — it is a real, load-bearing tension, and
    it is discussed in the report."""
    m = _repo_closure_module()
    with pytest.raises(ValueError):
        m.SuspensionClosure(eta_s=0.004736, phi0=0.70, phiM=0.5176,
                            sigma_star=44.2327, delta=1.3852, phi=0.44,
                            tau_relax=2.0, name="silica")


# ==========================================================================
# 4. CROSS-CHECK — two physically different routes to the same number
# ==========================================================================
def test_drag_force_integral_reproduces_CSG22_eq14_prefactor():
    """Integrate Eq. [13] over the contact area ourselves and compare with the
    printed drag force, CSG22 Eq. [14]:  F_d = 5.35 eta^1/2 rho^1/2 V^2 D sqrt(t).

    This is the strongest test in the file: it checks the shear-stress field,
    the screening cut-off, and the mean-stress normalisation all at once.

    Expected residual ~11%: the published 5.35 corresponds to an effective
    cut-off Delta = 0.0144 sqrt(VDt) rather than the 0.02 quoted in the text.
    Because the integral depends on Delta only logarithmically, a 40% change in
    Delta moves F_d by 11%. That is documented, not swept away.
    """
    imp = make()
    t = 3.0e-5
    d = float(st.contact_radius(imp, t))
    r_out = d * (1.0 - st.DELTA_OVER_D)
    r = np.linspace(0.0, r_out, 400_001)
    tau = st.wall_shear_stress(imp, t, r)
    F_num = 2.0 * np.pi * np.trapezoid(np.nan_to_num(tau) * r, r)

    F_paper = 5.35 * np.sqrt(imp.drop.eta) * np.sqrt(imp.drop.rho) * imp.U**2 * imp.D * np.sqrt(t)
    assert_allclose(F_num, F_paper, rtol=0.15)
    assert 0.85 < F_num / F_paper < 0.95   # pin the known 11% residual


def test_mean_shear_equals_drag_over_area():
    """<tau> must equal (published drag force)/(contact area) by construction."""
    imp = make()
    for ts in [0.05, 0.3, 1.0]:
        t = ts * imp.tau
        d = float(st.contact_radius(imp, t))
        F = 5.35 * np.sqrt(imp.drop.eta * imp.drop.rho) * imp.U**2 * imp.D * np.sqrt(t)
        assert_allclose(float(st.wall_shear_mean(imp, ts)), F / (np.pi * d**2), rtol=1e-10)


def test_compressible_crossover_two_independent_routes_agree():
    """The end of the compressible stage, estimated two ways:
       (a) geometric: contact edge slows to the sound speed
       (b) matching:  Wagner axis pressure falls to the water-hammer pressure
    They share no algebra. They must agree to O(1) - here, a factor 1.64."""
    x = st.wagner_hammer_crossover(make())
    ratio = x["t_star_shock_detach"] / x["t_star_pressure_match"]
    assert_allclose(ratio, 0.5 / (3.0 / np.pi**2), rtol=1e-12)
    assert 1.0 < ratio < 2.0


def test_frictional_contour_inverts_the_shear_field():
    """xi_c from Eq. (17) must be exactly where tau equals sigma*."""
    imp = make()
    for ts in [0.05, 0.5, 2.0]:
        for s in [5.0, 44.2, 500.0]:
            xi_c = float(dl.frictional_contour(imp, s, ts))
            if 0.0 < xi_c < 1.0:
                assert_allclose(float(st.wall_shear_selfsimilar(imp, ts, xi_c)), s, rtol=1e-9)


def test_area_fraction_monotone_in_velocity_and_bounded():
    frac = [float(dl.frictional_area_fraction(make(U=U), SUSP.sigma_star, 1.0))
            for U in [0.2, 0.5, 1.0, 2.0, 4.0]]
    assert all(0.0 <= f <= 1.0 for f in frac)
    assert all(b >= a for a, b in zip(frac, frac[1:]))


def test_onset_velocity_scales_as_sigma_star_to_the_two_thirds():
    """U_c ~ sigma*^(2/3) is the falsifiable prediction; check the code obeys it."""
    U1 = dl.onset_velocity(SUSP, R_LAB, sigma_star=20.0)
    U2 = dl.onset_velocity(SUSP, R_LAB, sigma_star=160.0)
    assert_allclose(U2 / U1, 8.0 ** (2.0 / 3.0), rtol=1e-6)


def test_shear_gauge_is_smaller_than_inertial_gauge_by_sqrt_Re():
    imp = make()
    ratio = dl.lambda_shear(imp, 44.2, 1.0) / dl.lambda_inertial(imp, 44.2)
    assert_allclose(ratio, (10.7 / (3.0 * np.pi)) * imp.Re_R**-0.5, rtol=1e-12)


# ==========================================================================
# 5. EDGE CASES AND GUARDS
# ==========================================================================
def test_pressure_is_nan_outside_the_contact_patch():
    imp = make()
    t = 1e-5
    d = float(st.contact_radius(imp, t))
    assert np.isnan(st.wagner_pressure(imp, t, 1.5 * d))
    assert np.isnan(st.wall_shear_stress(imp, t, 1.5 * d))


def test_shear_stress_vanishes_on_the_axis():
    """Required by axisymmetry — CSG22 make the point explicitly."""
    imp = make()
    assert_allclose(float(st.wall_shear_stress(imp, 1e-5, 0.0)), 0.0, atol=1e-12)


def test_wyart_cates_rejects_impossible_parameters():
    with pytest.raises(ValueError):
        mat.WyartCates(name="bad", eta_s=1e-3, phi0=0.6, phiM=0.5, sigma_star=10.0,
                       delta=1.0, phi=0.9, a=4e-7, rho_p=2000, rho_f=1000, gamma=0.072)
    with pytest.raises(ValueError):
        mat.WyartCates(name="bad", eta_s=1e-3, phi0=0.6, phiM=0.5, sigma_star=-1.0,
                       delta=1.0, phi=0.44, a=4e-7, rho_p=2000, rho_f=1000, gamma=0.072)


def test_suspension_density_matches_repo_pinned_value():
    """The repo pins rho_susp = 1440 kg/m^3 at phi = 0.44
    (tests/test_collapse_hunt.py:172)."""
    assert_allclose(mat.lab_silica(phi=0.44).rho, 1440.0, rtol=1e-12)


def test_wood_sound_speed_between_phase_extremes_and_below_both():
    """Wood's law for solid-in-liquid gives a speed between the two phases'
    -- and for silica in water it must exceed water's, since silica is stiffer
    per unit density here. Guard against a sign/reciprocal slip."""
    c = mat.lab_silica().c_wood
    assert 1000.0 < c < 5900.0
    assert c > mat.WATER.c


def test_no_nans_in_the_headline_report():
    rep = dl.dilatancy_report(make(), SUSP)
    assert all(np.isfinite(v) for k, v in rep.items() if isinstance(v, float))


# ==========================================================================
# 6. DIMENSIONAL ANALYSIS (pint) — catches an algebra slip that scaling misses
# ==========================================================================
def test_key_formulas_are_dimensionally_consistent():
    pint = pytest.importorskip("pint")
    u = pint.UnitRegistry()
    rho = 1440 * u.kg / u.m**3
    U = 2.3 * u.m / u.s
    R = 1.25 * u.mm
    eta = 0.0343 * u.Pa * u.s
    c = 1600 * u.m / u.s
    t = 30 * u.us

    assert (rho * c * U).to(u.Pa).dimensionality == u.Pa.dimensionality          # water hammer
    assert (rho * U**2).to(u.Pa).dimensionality == u.Pa.dimensionality           # dynamic
    d = np.sqrt(3 * R * U * t)
    assert d.to(u.m).dimensionality == u.m.dimensionality                        # contact radius
    p = rho * U**2 * R / d
    assert p.to(u.Pa).dimensionality == u.Pa.dimensionality                      # Wagner
    tau = np.sqrt(rho * eta * U**3 / R)
    assert tau.to(u.Pa).dimensionality == u.Pa.dimensionality                    # shear gauge
    Re_R = (rho * U * R / eta).to_base_units()
    assert Re_R.dimensionless


# ==========================================================================
# 7. THE INERTIAL-DILATANT NUMBER (Mobaseri et al., arXiv:2601.12642)
# ==========================================================================
import inertial_dilatant as idn  # noqa: E402


def test_ID_two_routes_agree():
    """ID = F_i/F_d (Eq. 6) must equal t_d/t_i, the paper's own restatement."""
    for sysx in [idn.cornstarch_M26(0.42), idn.lab_silica_dilatant(0.44, 0.42)]:
        for U0 in [0.5, 2.3, 7.2]:
            assert_allclose(idn.ID(sysx, U0, 4e-3),
                            idn.ID_from_times(sysx, U0, 4e-3), rtol=1e-12)


def test_ID_scales_linearly_with_U0_and_as_a_squared():
    """The two dependences that carry the whole argument."""
    s1 = idn.lab_silica_dilatant(0.44, 0.42)
    assert_allclose(idn.ID(s1, 4.6, 2.5e-3) / idn.ID(s1, 2.3, 2.5e-3), 2.0, rtol=1e-12)
    # ID ~ kappa ~ a^2, at fixed phi, phi_c, rho
    big = idn.DilatantSystem(name="x", phi=0.44, phi_c=0.42, a=8e-6, rho=1440.0)
    small = idn.DilatantSystem(name="y", phi=0.44, phi_c=0.42, a=0.8e-6, rho=1440.0)
    assert_allclose(idn.ID(big, 2.3, 2.5e-3) / idn.ID(small, 2.3, 2.5e-3), 100.0, rtol=1e-12)


def test_carman_kozeny_permeability_matches_hand_value():
    """kappa = (1-phi)^3 a^2/(180 phi^2), checked by hand for cornstarch."""
    s = idn.cornstarch_M26(0.40)
    hand = (1 - 0.40) ** 3 * (10e-6) ** 2 / (180 * 0.40**2)
    assert_allclose(s.kappa, hand, rtol=1e-12)
    assert_allclose(hand, 7.5e-13, rtol=1e-3)


def test_friction_coefficient_A_at_the_papers_beta0():
    """A = [(1-sinB)/(1+sinB)] exp(pi tanB); B0 = 44 deg gives A ~ 3.7-4.0."""
    s = idn.cornstarch_M26(0.42)
    assert 3.5 < s.A < 4.5


def test_no_dilatant_branch_below_phi_c():
    """phi < phi_c: there is no Reynolds dilatancy, so no dilatant force, so ID
    is infinite (inertia is unopposed). Must not silently return a number."""
    s = idn.lab_silica_dilatant(phi=0.44, phi_c=0.47)
    assert s.dphi < 0
    assert np.isinf(idn.ID(s, 2.3, 2.5e-3))
    assert idn.dilatant_pressure(s, 2.3, 2.5e-3) == 0.0


def test_darcy_pressure_for_cornstarch_matches_measured_order_of_magnitude():
    """Jerome, Vandenberghe & Forterre (PRL 117, 098003) MEASURED -10 to -100 kPa
    of pore suction under impact into a dense glass-bead bath. Our Eq.-(3)
    estimate for cornstarch must land in that decade — an absolute check
    against a direct measurement, not against another formula."""
    s = idn.cornstarch_M26(0.40)
    P = idn.dilatant_pressure(s, 2.0, 4e-3)
    assert 1e4 < P < 1e6, f"got {P:.3g} Pa"


def test_silica_darcy_pressure_massively_exceeds_its_own_ceiling():
    """The finding that motivates the capped extension: Eq. (3) at a = 800 nm
    returns pore suctions two to three decades above what the pore fluid can
    sustain (cavitation ~1 atm; interface pull-in 2 gamma/a)."""
    s = idn.lab_silica_dilatant(0.44, 0.42)
    r = idn.report(s, 2.3, 2.5e-3)
    assert r["darcy_overshoot"] > 50.0
    assert r["P_ceiling_Pa"] <= idn.P_ATM


def test_ID_capped_is_size_independent_once_the_cap_binds():
    """Once the ceiling binds, ID_cap = rho U^2/(pi A P_ceiling) — no kappa, no a."""
    a1 = idn.DilatantSystem(name="a1", phi=0.44, phi_c=0.42, a=0.8e-6, rho=1440.0)
    a2 = idn.DilatantSystem(name="a2", phi=0.44, phi_c=0.42, a=0.4e-6, rho=1440.0)
    # both must be cap-bound at these conditions
    assert idn.report(a1, 2.3, 2.5e-3)["darcy_overshoot"] > 1
    assert idn.report(a2, 2.3, 2.5e-3)["darcy_overshoot"] > 1
    # a2's ceiling is capillary-limited only if 2g/a < 1 atm; both here are 1 atm
    assert_allclose(idn.ID_capped(a1, 2.3, 2.5e-3), idn.ID_capped(a2, 2.3, 2.5e-3), rtol=1e-12)


def test_silica_is_far_deeper_in_the_dilatant_regime_than_cornstarch():
    """THE HEADLINE RATIO. Prefactor-free, since the unknown O(1) constants in
    F_d cancel between two systems evaluated the same way."""
    c = idn.compare(idn.lab_silica_dilatant(0.44, 0.42), idn.cornstarch_M26(0.42),
                    U0=2.3, D0=3e-3)
    assert c["ratio_a_over_b"] < 0.05, c
    assert 100 < 1.0 / c["kappa_ratio"] < 400


def test_reentrant_velocity_capped_is_reachable_uncapped_is_not():
    """The capped model's falsifiable prediction: a re-entrant liquid regime for
    silica at tens of m/s. The uncapped Eq. (6) puts it out of reach entirely."""
    s = idn.lab_silica_dilatant(0.44, 0.42)
    U_cap = idn.reentrant_velocity(s, 2.5e-3, capped=True)
    U_unc = idn.reentrant_velocity(s, 2.5e-3, capped=False)
    assert 5.0 < U_cap < 200.0, U_cap
    assert U_unc > 10 * U_cap


# ==========================================================================
# 8. TESTS ADDED AFTER THE ADVERSARIAL REVIEW
#    Every one of these would have caught a bug that shipped.
# ==========================================================================
def test_peak_force_recovers_the_universal_high_Re_values():
    """CAUGHT A REAL BUG: CSG22 Eq. [4] was transcribed inverted.
    At Re -> infinity the quartic's root must be ttil_max = 0.17 and the peak
    force must be the pinned 0.85."""
    fast = st.Impact(drop=mat.Fluid(name="inviscid", rho=1000.0, c=1500.0,
                                    eta=1e-9, gamma=0.072),
                     U=2.0, R=1e-3, target=st.rigid_target())
    pf = st.peak_force(fast)
    assert_allclose(pf["t_tilde_max"], 0.17, rtol=0.02)
    assert_allclose(pf["F_tilde_max"], 0.85, rtol=0.02)


def test_peak_force_rises_as_Re_falls():
    """CSG22 Section 2.1: 'When Re is below 200, Fmax increases and tmax
    decreases with decreasing Re.' The two must move OPPOSITELY."""
    Fs, ts = [], []
    for eta in [1e-4, 1e-3, 1e-2, 5e-2]:
        i = make(drop=mat.Fluid(name="x", rho=1440.0, c=1600.0, eta=eta, gamma=0.072))
        pf = st.peak_force(i)
        Fs.append(pf["F_tilde_max"])
        ts.append(pf["t_tilde_max"])
    assert all(b > a for a, b in zip(Fs, Fs[1:])), Fs      # F rises as eta rises
    assert all(b < a for a, b in zip(ts, ts[1:])), ts      # t falls as eta rises
    assert_allclose([f * t for f, t in zip(Fs, ts)],
                    [st.F_MAX_TILDE * st.T_MAX_TILDE] * len(Fs), rtol=1e-12)


def test_pressure_integrates_to_the_published_impact_force():
    """CROSS-CHECK the reviewers flagged as missing and free: integrating the
    Wagner pressure field over the wetted disc must give CSG22 Eq. [1],
    Ftilde = (3 sqrt6 / 2) ttil^(1/2), prefactor included. This checks Eq. [10]
    and Eq. [1] against each other with no shared algebra."""
    imp = make()
    for ts in [1e-4, 1e-3, 1e-2]:
        t = ts * imp.tau
        d = float(st.contact_radius(imp, t))
        # integrate in theta with r = d sin(theta): the 1/sqrt(d^2-r^2)
        # singularity is integrable but a uniform mesh converges only as
        # sqrt(N), so change variables. The FIELD is still the coded one.
        th = np.linspace(0.0, np.pi / 2 * (1 - 1e-12), 200_001)
        r = d * np.sin(th)
        integ = np.nan_to_num(st.wagner_pressure(imp, t, r)) * r * d * np.cos(th)
        F = 2 * np.pi * np.trapezoid(integ, th)
        F_paper = (3 * np.sqrt(6) / 2) * np.sqrt(ts / 2) * imp.drop.rho * imp.D**2 * imp.U**2
        assert_allclose(F, F_paper, rtol=2e-3)


def test_shear_field_from_an_independent_stokes_layer_argument():
    """The strongest kind of test: rebuild CSG22 Eq. [13] from different physics.

    Wall stress of a Stokes first problem is tau = eta u / sqrt(pi nu t_loc).
    Take u = the Wagner disc slip velocity (2U/pi) r/sqrt(d^2-r^2), and
    t_loc = (d^2-r^2)/(3RU), the age of the contact line at radius r. No part of
    this shares algebra with the transcribed Eq. [13]."""
    imp = make()
    for ts in [1e-3, 1e-2, 1e-1]:
        t = ts * imp.tau
        d = float(st.contact_radius(imp, t))
        r = np.array([0.2, 0.5, 0.8, 0.95]) * d
        s = np.sqrt(d**2 - r**2)
        u_slip = (2.0 * imp.U / np.pi) * r / s
        t_loc = s**2 / (3.0 * imp.R * imp.U)
        tau_indep = imp.drop.eta * u_slip / np.sqrt(np.pi * imp.drop.nu * t_loc)
        assert_allclose(st.wall_shear_stress(imp, t, r), tau_indep, rtol=1e-12)


def test_shear_vanishes_on_a_bath_of_the_same_liquid():
    """The sharpest prediction in the module, and it must not be an accident of
    floating point: matched densities give an exact symmetry plane."""
    same = st.liquid_target(DROP)
    assert st.shear_sharing_factor(make(target=same))["slip_ratio"] == 0.0
    assert st.shear_sharing_factor(make(target=st.rigid_target()))["slip_ratio"] == 1.0
    denser = st.liquid_target(mat.Fluid(name="2x", rho=2 * DROP.rho, c=1500.0,
                                        eta=1e-3, gamma=0.072))
    assert_allclose(st.shear_sharing_factor(make(target=denser))["slip_ratio"],
                    1.0 / 3.0, rtol=1e-12)


def test_inertial_factor_is_withheld_for_a_solid_target():
    """CAUGHT A REAL BUG: the inertia-sharing law was being reported for elastic
    solids, where it does not apply and is wrong by ~2.5x."""
    solid = st.compliance_factors(make(target=st.elastic_target(1030.0, 1e6)))
    assert np.isnan(solid["inertial_factor"])
    assert "NOT applicable" in solid["applicable_law"]
    liquid = st.compliance_factors(make(target=st.liquid_target(mat.WATER)))
    assert np.isfinite(liquid["inertial_factor"])


def test_glass_is_not_acoustically_rigid():
    """'Glass is within 0.1% of rigid' was wrong: it takes ~15% off."""
    f = st.compliance_factors(make(target=st.glass_target()))["acoustic_factor"]
    assert 0.80 < f < 0.90


def test_compliance_threshold_reproduces_howlands_worked_example():
    """Howland et al. PRL 117, 184502 (2016): E_c = 93 kPa for R = 1.5 mm,
    U = 2.36 m/s, t_ej - t_i = 30 us, rho = 789 kg/m^3. Their criterion
    p_peak >~ E with p_peak = rho R U/(t_ej-t_i) is algebraically C >~ t*_ej."""
    rho, R, U, dt = 789.0, 1.5e-3, 2.36, 30e-6
    assert_allclose(rho * R * U / dt, 93.1e3, rtol=0.01)
    C = rho * U**2 / (rho * R * U / dt)      # = Cauchy number at E = p_peak
    assert_allclose(C, U * dt / R, rtol=1e-12)   # = t*_ej


def test_ID_capped_reduces_to_ID_when_the_cap_does_not_bind():
    """CAUGHT A REAL BUG: ID_capped carried a pi that Eq. 6 does not, so the two
    disagreed by exactly pi in the limit where they must agree."""
    big = idn.DilatantSystem(name="coarse", phi=0.44, phi_c=0.42, a=200e-6, rho=1440.0)
    assert idn.report(big, 2.3, 2.5e-3)["darcy_overshoot"] < 1.0   # cap inactive
    assert_allclose(idn.ID_capped(big, 2.3, 2.5e-3), idn.ID(big, 2.3, 2.5e-3), rtol=1e-12)


def test_pore_throat_ceiling_never_binds_for_either_system():
    """The corrected capillary criterion uses the pore THROAT (0.155 R), not the
    particle size. Both systems then sit above 1 atm, so the ceiling is always
    cavitation - which is what removes the size dependence from ID_capped."""
    for s in [idn.cornstarch_M26(0.41), idn.lab_silica_dilatant(0.44, 0.42)]:
        assert s.p_capillary > idn.P_ATM
        assert_allclose(s.p_dilatant_ceiling, idn.P_ATM, rtol=1e-12)


def test_poroelastic_ratio_separates_cornstarch_from_colloids():
    """Jerome et al.'s own validity criterion for their Darcy-Reynolds equation.
    Cornstarch must pass comfortably; 800 nm silica must not."""
    cs = idn.poroelastic_time_ratio(idn.cornstarch_M26(0.41), 5.1, 4e-3)
    si = idn.poroelastic_time_ratio(idn.lab_silica_dilatant(0.44, 0.42), 2.3, 2.5e-3)
    assert cs["darcy_incompressible_valid"] and cs["ratio"] < 0.05
    assert not si["darcy_incompressible_valid"] and si["ratio"] > 0.5
    assert si["ratio"] / cs["ratio"] > 30


def test_onset_velocity_raises_rather_than_returning_a_bare_nan():
    with pytest.raises(ValueError):
        dl.onset_velocity(SUSP, R_LAB, sigma_star=1e-6, U_bracket=(1.0, 2.0))


def test_stagnation_term_is_not_silently_dropped_when_asked_for():
    imp = make()
    bare = float(st.wagner_pressure_axis(imp, 0.315))
    full = float(st.wagner_pressure_axis(imp, 0.315, include_stagnation=True))
    assert_allclose(full - bare, 0.5 * imp.p_dyn, rtol=1e-12)
    # the omitted term is ~50% OF THE SINGULAR PART at t* ~ 0.3 (34% of the total)
    assert 0.45 < (full - bare) / bare < 0.55
