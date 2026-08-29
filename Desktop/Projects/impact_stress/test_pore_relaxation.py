"""
test_pore_relaxation.py — validation suite for pore_relaxation.py.

Same four kinds of test as test_impact_stress.py, in increasing order of what
they prove:

  1. LIMITS      — does the formula reduce to the known answer in a limit?
  2. IDENTITY    — do two algebraically-equivalent forms agree numerically?
  3. LITERATURE  — does an independently-known printed value come back out?
  4. CROSS-CHECK — does a quantity computed two physically different ways agree,
                   or agree with an independent MEASUREMENT?

The tests that carry the most weight here are the two in
`TestAgainstMeasurement`: the stiffness ceiling on the dilatant pressure,
Eq. (D3), is a new claim in this project and it is checked against Jerome,
Vandenberghe & Forterre's measured pore suction in glass beads and against
their own extrapolation to cornstarch.

Run:  python -m pytest test_pore_relaxation.py -q
"""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

import inertial_dilatant as idl
import pore_relaxation as pr

D0_LAB = 2.5e-3     # m, lab drop diameter
U0_LAB = 2.3        # m/s, band B centre
D0_M26 = 4.0e-3     # m, [M26] cornstarch drop
U0_M26 = 5.1        # m/s, [M26] Regime-2/3 boundary


# ==========================================================================
# 1. LIMITS
# ==========================================================================
class TestLimits:

    def test_rigid_skeleton_recovers_fluid_only_storage(self):
        """K_p -> inf must leave only the pore-fluid compressibility."""
        phi = 0.50
        S_rigid = pr.storage_coefficient(phi, K_p=1e30)
        assert_allclose(S_rigid, (1 - phi) * pr.CHI_WATER, rtol=1e-10)

    def test_soft_skeleton_recovers_collective_diffusion(self):
        """K_p small => D_p -> kappa K_p / mu, the collective diffusivity."""
        phi, kappa, K_p = 0.50, 1.8e-15, 1.0e3
        D = pr.pore_diffusivity(kappa, K_p, phi)
        assert_allclose(D, kappa * K_p / pr.MU_W, rtol=1e-3)

    def test_storage_terms_are_additive_and_series_dominated_by_softer(self):
        """The larger compliance wins; for every system here that is 1/K_p."""
        p = pr.silica_poro(phi=0.50)
        S = pr.storage_coefficient(0.50, p.K_p)
        assert (1.0 / p.K_p) / S > 0.99, "network should dominate storage"

    def test_tau_p_at_particle_length_is_size_independent(self):
        """L = a makes tau_p independent of a, because kappa ~ a^2.

        This is the counterexample to 'smaller particles always drain slower'
        and it is why the drainage LENGTH has to be stated with every number.
        """
        phi, K_p = 0.50, 1.0e3
        taus = []
        for a in (800e-9, 10e-6, 170e-6):
            kappa = (1 - phi) ** 3 * a**2 / (180 * phi**2)
            taus.append(pr.pore_relaxation_time(a, kappa, K_p, phi))
        assert_allclose(taus, taus[0], rtol=1e-12)

    def test_no_dilatant_branch_below_onset(self):
        s = idl.lab_silica_dilatant(phi=0.40, phi_c=0.42)
        assert s.dphi < 0
        assert pr.dilatant_pressure_darcy(s, 1e-3, 1e3) == 0.0
        assert pr.dilatant_pressure_ceiling(s, 1e3) == 0.0


# ==========================================================================
# 2. IDENTITIES
# ==========================================================================
class TestIdentities:

    def test_D1_equals_D2_at_the_same_drainage_length(self):
        """p_dil = (mu/kappa) dphi L^2 gd  ==  K_p dphi gd tau_p(L).

        This is the check that the derivation linking the dilatancy closure to
        the pore-relaxation time was transcribed correctly. It holds only in the
        skeleton-dominated limit, so the tolerance carries the fluid-storage
        correction (which is ~1e-6 relative here).
        """
        p = pr.silica_poro(phi=0.50)
        s, K = p.sys, p.K_p
        L, gd = 0.49 * D0_LAB, U0_LAB / D0_LAB
        tau_p = p.tau_p(L)
        assert_allclose(pr.dilatant_pressure_darcy(s, L, gd),
                        pr.dilatant_pressure_from_identity(K, s, gd, tau_p),
                        rtol=1e-4)

    def test_capped_is_min_of_the_two(self):
        p = pr.cornstarch_poro(phi=0.42)
        s, K = p.sys, p.K_p
        L, gd = 0.49 * D0_M26, U0_M26 / D0_M26
        assert_allclose(pr.dilatant_pressure_capped(s, K, L, gd),
                        min(pr.dilatant_pressure_darcy(s, L, gd),
                            pr.dilatant_pressure_ceiling(s, K)), rtol=0)

    def test_ceiling_strain_is_referenced_to_phi_not_porosity(self):
        """p_max = K_p dphi/phi. Regression guard on a real bug.

        The first version of (D3) used dphi/(1-phi), which mixes conventions:
        with incompressible grains the particle volume phi*V is fixed, so
        eps_v = -dphi/phi, and the same phi appears in K_p = phi dP/dphi.
        The two forms coincide exactly at phi = 0.5, which is why the original
        error survived the first pass at the silica operating point.
        """
        K = 1.0e3
        s = idl.lab_silica_dilatant(phi=0.44, phi_c=0.42)
        assert_allclose(pr.dilatant_pressure_ceiling(s, K), K * 0.02 / 0.44,
                        rtol=1e-12)
        wrong = K * s.dphi / (1 - s.phi)
        assert_allclose(pr.dilatant_pressure_ceiling(s, K) / wrong,
                        (1 - 0.44) / 0.44, rtol=1e-12)
        # the two conventions are identical at phi = 0.5 and only there
        half = idl.lab_silica_dilatant(phi=0.50, phi_c=0.42)
        assert_allclose(pr.dilatant_pressure_ceiling(half, K),
                        K * half.dphi / (1 - half.phi), rtol=1e-12)

    def test_slowest_mode_is_the_documented_fraction_of_the_scaling(self):
        """tau_1 = (4/pi^2) L^2/D_p, about 2.5x faster than the scaling."""
        p = pr.silica_poro(phi=0.50)
        L = 0.49 * D0_LAB
        exact = pr.slowest_mode_time(L, p.kappa, p.K_p, p.sys.phi)
        assert_allclose(exact / p.tau_p(L), 4.0 / np.pi**2, rtol=1e-12)
        assert 0.40 < exact / p.tau_p(L) < 0.41

    def test_De_is_tau_ratio(self):
        assert_allclose(pr.drainage_number(2.0e-3, D0_LAB, U0_LAB),
                        2.0e-3 / (D0_LAB / U0_LAB), rtol=1e-12)

    def test_tau_p_scales_as_L_squared(self):
        p = pr.silica_poro(phi=0.50)
        assert_allclose(p.tau_p(2e-3) / p.tau_p(1e-3), 4.0, rtol=1e-12)

    def test_tau_p_inverse_in_Kp_when_network_dominates(self):
        p = pr.silica_poro(phi=0.50)
        assert_allclose(p.tau_p(1e-3, K_p=1e3) / p.tau_p(1e-3, K_p=1e4),
                        10.0, rtol=1e-3)


# ==========================================================================
# 3. LITERATURE VALUES
# ==========================================================================
class TestLiterature:

    def test_carnahan_starling_at_hard_sphere_freezing(self):
        """CS must reproduce the known hard-sphere freezing pressure.

        Molecular dynamics puts the freezing transition at phi = 0.494 with
        compressibility factor Z = betaP/n = 12.3. Carnahan-Starling is known to
        overestimate slightly there, by a few percent -- so agreement to 5% is
        the pass, and agreement to 0.1% would mean the reference is wrong.
        """
        phi = 0.494
        Z = (1 + phi + phi**2 - phi**3) / (1 - phi) ** 3
        assert_allclose(Z, 12.3, rtol=0.05)
        assert Z > 12.3, "CS is expected to sit slightly above the MD value"

    def test_osmotic_modulus_is_tiny_for_800nm_colloid(self):
        """A Brownian 800 nm hard-sphere colloid at phi = 0.5 has K ~ 1 Pa.

        Order-of-magnitude anchor: n kT = 7.6e-3 Pa at phi = 0.5 for 800 nm
        spheres, and CS multiplies that by O(10-100) for K = phi dPi/dphi.
        The point of the test is that it is NOT kPa, which is what makes
        colloidal consolidation take days.
        """
        K = pr.osmotic_modulus_hard_sphere(0.50, 800e-9)
        assert 0.05 < K < 50.0, f"K_osm = {K:g} Pa is off the expected scale"

    def test_debye_length_matches_the_standard_formula(self):
        """lambda_D[nm] = 0.304 / sqrt(I[M]) for a 1:1 electrolyte at 25 C.

        Computed here at 25 C to match the constant in the printed formula;
        the module's default is 20 C, a 0.8% difference.
        """
        lam = pr.debye_length(1.0, T=298.15)     # 1 mM
        assert_allclose(lam * 1e9, 0.304 / np.sqrt(1e-3), rtol=0.02)

    def test_carman_kozeny_approaches_1e_3_d2_at_random_close_packing(self):
        """kappa ~ 1e-3 d^2 at phi = 0.64 is the standard RCP value."""
        phi, d = 0.64, 1.0
        kappa = (1 - phi) ** 3 * d**2 / (180 * phi**2)
        assert_allclose(kappa, 1.0e-3, rtol=0.4)

    def test_permeability_ratio_silica_to_cornstarch_is_the_a_squared_ratio(self):
        """At matched phi the whole kappa contrast is (a_si/a_cs)^2 = 1/156."""
        phi = 0.42
        si = idl.lab_silica_dilatant(phi=phi, phi_c=0.40)
        cs = idl.cornstarch_M26(phi=phi)
        assert_allclose(si.kappa / cs.kappa, (800e-9 / 10e-6) ** 2, rtol=1e-12)


# ==========================================================================
# 4. CROSS-CHECKS AGAINST MEASUREMENT  (the load-bearing tests)
# ==========================================================================
class TestAgainstMeasurement:

    def test_ceiling_reproduces_JVF_measured_suction_in_glass_beads(self):
        """THE LOAD-BEARING TEST. [JVF16] measure -1e4 to -1e5 Pa of pore
        suction in 170 um glass beads. Eq. (D3), p_max = K_p dphi/(1-phi), must
        land within a decade of that.

        History, because it is the reason the module looks the way it does: with
        the pressure-free stiffness k_n ~ E_glass a this test returned 1.5e8 Pa,
        three decades high, and that failure is what forced the Hertzian
        pressure-dependent stiffness. The affine modulus is still an upper bound
        and dphi = 0.005 is still an estimate, so landing at or somewhat above
        the measured range is expected; landing below it, or two decades above,
        would falsify (D3).
        """
        p = pr.glass_beads_JVF()
        ceiling = pr.dilatant_pressure_ceiling(p.sys, p.K_p)
        assert 1.0e4 <= ceiling <= 1.0e6, (
            f"(D3) gives {ceiling:.3g} Pa for JVF's glass beads against a "
            "measured 1e4-1e5 Pa.")

    def test_pressure_free_stiffness_is_the_three_decade_error_it_is_labelled(self):
        """Guard the deprecation: k_n ~ E a must still be shown to be ~1e3 too
        stiff at JVF's confinement, so nobody quietly reinstates it."""
        a, phi, nu, E = 170e-6, 0.58, 0.25, 70.0e9
        k_free = pr.contact_stiffness_hertz_scale(E, a)
        k_hertz = pr.contact_stiffness_hertz(E, nu, a, phi, pr.P_CONF_JVF)
        assert 1.0e2 < k_free / k_hertz < 1.0e4

    def test_cornstarch_ceiling_sits_at_the_measured_impact_stress_scale(self):
        """CONSISTENCY, NOT VALIDATION -- read the docstring before citing it.

        [JVF16]'s closing paragraph extrapolates 1-10 MPa for cornstarch, but
        that extrapolation is itself uncapped Darcy, i.e. the very quantity (D3)
        exists to bound, so agreeing with it would prove nothing. The
        independent number is the jammed-column stress measured directly under
        impact by Waitukaitis & Jaeger, Nature 487, 205 (2012), which is
        1e5-1e6 Pa. (D3) must land in or near that band.
        """
        p = pr.cornstarch_poro(phi=0.42)
        ceiling = pr.dilatant_pressure_ceiling(p.sys, p.K_p)
        assert 1.0e4 <= ceiling <= 1.0e7, (
            f"(D3) gives {ceiling:.3g} Pa for cornstarch against a measured "
            "impact stress scale of 1e5-1e6 Pa.")

    def test_ceiling_is_far_below_uncapped_darcy_for_silica(self):
        """The whole reason (D3) exists: uncapped Darcy gives 1e7-1e8 Pa for
        800 nm silica, which no pore fluid can reach."""
        p = pr.silica_poro(phi=0.50)
        gd = U0_LAB / D0_LAB
        darcy = pr.dilatant_pressure_darcy(p.sys, 0.49 * D0_LAB, gd)
        ceiling = pr.dilatant_pressure_ceiling(p.sys, p.K_p)
        assert darcy > 1.0e6, f"expected the known 1e7-1e8 Pa problem, got {darcy:.3g}"
        assert ceiling < darcy / 1.0e3, "the ceiling must actually bind for silica"

    def test_silica_is_undrained_and_cornstarch_is_marginal(self):
        """The headline comparison, as an assertion.

        At the contact-patch drainage length, 800 nm silica must be deep in the
        undrained limit (De >> 1) while 10 um cornstarch sits near De ~ 1.
        """
        si = pr.silica_poro(phi=0.50)
        cs = pr.cornstarch_poro(phi=0.42)
        De_si = pr.drainage_number(si.tau_p(0.49 * D0_LAB), D0_LAB, U0_LAB)
        De_cs = pr.drainage_number(cs.tau_p(0.49 * D0_M26), D0_M26, U0_M26)
        assert De_si > 1.0e3, f"silica De = {De_si:.3g}, expected >> 1"
        assert De_cs < 1.0e2, f"cornstarch De = {De_cs:.3g}, expected O(1)"
        assert De_si / De_cs > 1.0e3

    def test_module_agrees_with_existing_rigid_skeleton_result(self):
        """Cross-check against inertial_dilatant.poroelastic_time_ratio.

        Setting K_p = inf and dropping the porosity factor must reproduce the
        2026-08-06 brief's tau_d/tau_i = 0.81 for silica. The two differ ONLY by
        the (1-phi) that this module carries and the older one does not, so the
        ratio of the two answers must equal 1/(1-phi).
        """
        phi = 0.44
        s = idl.lab_silica_dilatant(phi=phi)
        old = idl.poroelastic_time_ratio(s, U0_LAB, D0_LAB)
        ratio_old = old["ratio"] if isinstance(old, dict) else old
        tau_new = pr.pore_relaxation_time(D0_LAB, s.kappa, K_p=1e30, phi=phi)
        ratio_new = pr.drainage_number(tau_new, D0_LAB, U0_LAB)
        assert_allclose(ratio_new / ratio_old, 1.0 - phi, rtol=1e-6)


# ==========================================================================
# 5. EDGES AND GUARDS
# ==========================================================================
class TestEdges:

    @pytest.mark.parametrize("phi", [0.0, 1.0, -0.1, 1.5])
    def test_phi_out_of_range_raises(self, phi):
        with pytest.raises(ValueError):
            pr.storage_coefficient(phi, K_p=1e3)

    @pytest.mark.parametrize("bad", [0.0, -1.0])
    def test_nonpositive_moduli_and_lengths_raise(self, bad):
        with pytest.raises(ValueError):
            pr.storage_coefficient(0.5, K_p=bad)
        with pytest.raises(ValueError):
            pr.pore_relaxation_time(bad, 1e-15, 1e3, 0.5)
        with pytest.raises(ValueError):
            pr.pore_diffusivity(bad, 1e3, 0.5)

    def test_zero_salt_raises_rather_than_returning_infinity(self):
        with pytest.raises(ValueError):
            pr.debye_length(0.0)

    def test_report_is_finite_everywhere_it_should_be(self):
        p = pr.silica_poro(phi=0.50)
        rep = pr.report(p, D0_LAB, U0_LAB)
        for k, v in rep.items():
            if isinstance(v, float):
                assert np.isfinite(v), f"{k} is not finite"
            elif isinstance(v, dict):
                for kk, vv in v.items():
                    if isinstance(vv, float):
                        assert np.isfinite(vv), f"{k}.{kk} is not finite"

    def test_salt_stiffens_the_network(self):
        """k_n ~ 1/lambda_D, so more salt = shorter screening = stiffer contact.

        A prediction of the salt-screening series, asserted so it cannot be
        silently inverted by a later edit, and checked on both the standalone
        stiffness and the packaged silica system.
        """
        soft = pr.contact_stiffness_repulsive(44.2327, 800e-9, lam=30e-9)
        stiff = pr.contact_stiffness_repulsive(44.2327, 800e-9, lam=1e-9)
        assert_allclose(stiff / soft, 30.0, rtol=1e-12)
        p = pr.silica_poro(phi=0.50)
        assert p.K_p_high > p.K_p > p.K_p_low
        assert_allclose(p.K_p_high / p.K_p_low, 30.0, rtol=1e-12)

    def test_edl_modulus_is_independent_of_Z_and_phi(self):
        """K_p = P a/(6 lambda_D): Z and phi cancel exactly. Checked against the
        two-step affine route, which still has both in it."""
        P, a, phi, lam = 44.2327, 800e-9, 0.52, 10e-9
        direct = pr.network_modulus_edl(P, a, lam)
        two_step = pr.network_modulus_affine(
            pr.contact_stiffness_edl(P, a, phi, lam), a, phi)
        assert_allclose(direct, two_step, rtol=1e-12)
        assert_allclose(direct, pr.network_modulus_edl(P, a, lam), rtol=0)
