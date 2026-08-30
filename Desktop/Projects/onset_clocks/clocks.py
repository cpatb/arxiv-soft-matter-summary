"""Candidate criteria for what arrests an impacting dense-suspension drop.

The problem, in the repo's own words
------------------------------------
``impact_stress/materials.py`` carries this as a caveat on the lab's own
calibration:

    At phi = 0.44 < phiM = 0.5176 the suspension can NEVER jam by this model,
    at any stress; it can only thicken continuously. The drops in the
    experiment demonstrably do arrest. That mismatch is real.

Two escapes have been offered. The first is that the volume-fraction axis of
the rheometer dataset is offset from the drops' by "an unquantified amount".
That offset is now quantified -- see
``tier0_wyart_cates_closure/CALIBRATION_2026-08-30.md`` -- and it moves ``phi_m``
the *wrong way*, so it cannot rescue the model. The second is that the impact
stress is so large that the material thickens into arrest anyway. That is what
this module exists to test, along with every other candidate on the table.

The method: eliminate by degeneracy, not by correlation
-------------------------------------------------------
Each candidate is a dimensionless group ``G`` computed per drop from
independently-measured inputs. If ``G`` is what sets the arrested footprint
``beta_p``, then two drops with the same ``G`` must have the same ``beta_p``.
So the test is not "does ``G`` correlate with ``beta_p``" -- with 75 drops and a
strong velocity axis nearly everything correlates -- but:

    among drops whose G agree to within its own resolution, how much of the
    total spread in beta_p survives?

A group that leaves the full spread among its own ties is falsified regardless
of how well it correlates. This is the "compare a detector against its own
quantum" rule, applied to a theory rather than to an instrument.

Units
-----
SI throughout: m, s, kg, Pa, Pa.s, kg/m^3. Every group returned is
dimensionless and is asserted so by the tests.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

IMPACT_STRESS = Path(r"C:\Users\compl\Desktop\Projects\impact_stress")
sys.path.insert(0, str(IMPACT_STRESS))

from conventions import solve_Lc  # noqa: E402
from materials import lab_silica, suspension_as_fluid  # noqa: E402
from pore_relaxation import (  # noqa: E402
    debye_length,
    drainage_lengths,
)
from stress import Impact, glass_target, peak_force, wall_shear_mean  # noqa: E402

# Added-mass coefficient of an impact-driven jamming front, measured:
# k = 2 to 12 across Waitukaitis & Jaeger, Nature 487, 205 (2012); Peters,
# Majumdar & Jaeger, Nature 532, 214 (2016); Han, Peters & Jaeger.
K_ADDED_MASS_LO, K_ADDED_MASS_HI = 2.0, 12.0

# Boltzmann constant, CODATA 2018 (exact).
K_B = 1.380649e-23
T_LAB = 293.15  # K; no temperature is recorded with these recordings

# Largest sigma*(0 mM)/sigma*(100 mM) ratio the repo's own collapse hunt was
# able to invert out of the salt series, in
# Github/Drop Impact/analysis_2026-07-27/collapse_scoreboard.csv. Used only as
# a generous upper bound on how far salt could possibly move sigma*, so that
# the stress clock is eliminated against the most favourable case for it rather
# than against the nominal one.
SIGMA_STAR_SALT_SPAN = 7.53


@dataclass(frozen=True)
class Drop:
    """One impact, from the per-drop table. All SI except phi and salt."""

    label: str
    phi: float
    salt_mM: float
    u0: float        # m/s
    d0: float        # m
    beta_p: float    # arrested footprint / d0, dimensionless
    beta_eq: float   # relaxed footprint / d0


def material(phi: float):
    """The lab's own 800 nm silica at this packing, with its own caveats.

    Reused from ``impact_stress.materials`` rather than re-fitted, so that
    everything here is consistent with the calibration the drop-impact repo
    already uses, including its documented failure of its own fit gate.
    """
    return lab_silica(phi=phi)


# --------------------------------------------------------------------------
# stress clocks
# --------------------------------------------------------------------------
def impact_stress_scale(d: Drop, k: float) -> float:
    """Sigma = rho k U0^2 [Pa]. k = 1/2 recovers the bare dynamic pressure."""
    return material(d.phi).rho * k * d.u0**2


def wall_shear_gauge(d: Drop) -> float:
    """Area-averaged Wagner wall shear stress at peak force [Pa].

    The **deviatoric** gauge, and the right one. ``rho k U0^2`` is a total
    stress; under undrained loading -- which this silica is, by five to six
    decades (``impact_stress/PORE_PHYSICS.md``) -- the pore fluid carries the
    isotropic part and only the deviatoric part generates particle pressure.
    ``impact_stress/STRESS_AND_DILATANCY_BRIEF.md`` establishes this for exactly
    this geometry and states the consequence: rho U0^2 runs 40 to 800x above
    sigma* across the whole accessible velocity range and therefore has no
    reachable threshold, while the wall shear is 5 to 10x smaller and does cross
    sigma* near the bottom of the window.

    Evaluated at the time of peak force, which is where the impact is doing the
    most, and which the same module computes from CSG22 Eqs. [3]-[4] rather than
    assuming the high-Re value.

    <tau> = C rho U^2 Re_R^(-1/2) t*^(-1/2), C = 10.7/(3 pi), band 0.76 to 1.14
    because the underlying integral is logarithmically divergent in the
    screening length. The band is carried by :data:`MEAN_SHEAR_COEFF_BAND`.
    """
    m = material(d.phi)
    imp = Impact(drop=suspension_as_fluid(m), U=d.u0, R=0.5 * d.d0,
                 target=glass_target())
    return float(wall_shear_mean(imp, peak_force(imp)["t_star_max"]))


def frictional_fraction_shear(d: Drop, sigma_star_factor: float = 1.0) -> float:
    """f at the deviatoric gauge rather than at rho k U0^2."""
    m = material(d.phi)
    ratio = (m.sigma_star * sigma_star_factor) / wall_shear_gauge(d)
    return float(np.exp(-(ratio ** m.delta)))


def shear_over_sigma_star(d: Drop) -> float:
    """The deviatoric gauge in units of the onset stress. Dimensionless."""
    return wall_shear_gauge(d) / material(d.phi).sigma_star


def sigma_over_sigma_star(d: Drop, k: float = K_ADDED_MASS_LO,
                          sigma_star_factor: float = 1.0) -> float:
    """Impact stress in units of the onset stress. Dimensionless."""
    m = material(d.phi)
    return impact_stress_scale(d, k) / (m.sigma_star * sigma_star_factor)


def frictional_fraction_at_impact(d: Drop, k: float = K_ADDED_MASS_LO,
                                  sigma_star_factor: float = 1.0) -> float:
    """f = exp(-(sigma*/Sigma)^delta) at the impact stress. Dimensionless, in (0, 1).

    ``sigma_star_factor`` scales sigma* so the clock can be tested against the
    most generous salt-driven shift the data could support, rather than only at
    its nominal value.
    """
    m = material(d.phi)
    ratio = (m.sigma_star * sigma_star_factor) / impact_stress_scale(d, k)
    return float(np.exp(-(ratio ** m.delta)))


def jamming_deficit(d: Drop) -> float:
    """f required for phi_J(f) to reach phi: (phi0 - phi)/(phi0 - phiM).

    Below 1 the material can jam at a finite stress; at or above 1 it cannot
    jam at any stress, because even a fully frictional network jams looser than
    the packing it is at. Dimensionless, and invariant under the model's
    ``(phi, phi_m, phi_0) -> c(...)`` gauge.
    """
    m = material(d.phi)
    return (m.phi0 - d.phi) / (m.phi0 - m.phiM)


# --------------------------------------------------------------------------
# rate and strain clocks
# --------------------------------------------------------------------------
def shear_rate_naive(d: Drop) -> float:
    """gammadot = U0/D0 [1/s], convention (A) of impact_stress/conventions.py.

    This is the convention the lab's own gammadot_c = 632 1/s is quoted on.
    """
    return d.u0 / d.d0


def eta_of_rate(m):
    """Callable gammadot [1/s] -> viscosity [Pa.s] on the model's steady branch.

    Wyart-Cates is written stress-controlled, so the rate-controlled viscosity
    is obtained by generating the steady flow curve and interpolating in log
    space. Legitimate here only because these packings are all below phiM, so
    the curve is single-valued in rate; :func:`_check_monotone` refuses if it
    ever is not.
    """
    # A grid wide enough that the L_c bisection, which brackets L_c/D0 over
    # [1e-4, 1], never asks for a rate outside it. np.interp clamps at the ends,
    # and a clamped viscosity silently breaks the root bracket rather than
    # raising, so the grid is widened instead of the clamp being tolerated.
    fc = m.flow_curve(sigma_grid=np.logspace(-8, 10, 4000))
    gd, eta = np.asarray(fc["gammadot"], float), np.asarray(fc["eta"], float)
    ok = np.isfinite(gd) & np.isfinite(eta) & (gd > 0) & (eta > 0)
    gd, eta = gd[ok], eta[ok]
    order = np.argsort(gd)
    gd, eta = gd[order], eta[order]
    if np.any(np.diff(gd) <= 0):
        raise ValueError(
            f"{m.name}: the flow curve is not single-valued in shear rate, so "
            "an eta(gammadot) does not exist. This packing back-bends.")

    def eta_fn(rate):
        return float(np.exp(np.interp(np.log(float(rate)), np.log(gd), np.log(eta))))

    return eta_fn


def lamella_solution(d: Drop) -> dict | None:
    """The Mobaseri PNAS-2025 self-consistent lamella solution, or None.

    The only characteristic length among the field's four competing shear-rate
    conventions that was validated against a measured dissipation integral
    rather than asserted. Solved self-consistently because Re depends on the
    viscosity at the rate L_c itself sets.

    Returns ``None`` where the root is not bracketed inside ``L_c < D0``. That
    is not a numerical failure: it happens when ``Re`` falls to order 1, where
    the correlation returns a "lamella" thicker than the drop, i.e. outside the
    regime it was fitted in. It is reported as undefined rather than
    extrapolated, and the count of such drops is part of the audit's result --
    a convention that does not cover the experiment cannot be what sets it.
    """
    m = material(d.phi)
    try:
        return solve_Lc(m.rho, d.u0, d.d0, m.gamma, eta_of_rate(m), m.eta_0)
    except ValueError:
        return None


def lamella_thickness(d: Drop) -> float:
    """L_c [m], or NaN where the correlation is out of its regime."""
    sol = lamella_solution(d)
    return float(sol["L_c_m"]) if sol else float("nan")


def shear_rate_lamella(d: Drop) -> float:
    """gammadot = U0/L_c [1/s], convention (D). Larger than (A) by D0/L_c."""
    lc = lamella_thickness(d)
    return d.u0 / lc if np.isfinite(lc) else float("nan")


def strain_in_one_impact(d: Drop, rate=shear_rate_naive) -> float:
    """Strain accumulated over the inertial impact time t = D0/U0. Dimensionless.

    For the naive convention this is identically 1 by construction, which is
    itself the point: a strain criterion built on U0/D0 has no velocity
    dependence at all and cannot produce a velocity onset. The lamella
    convention does vary, because L_c/D0 depends on Re and We.
    """
    return rate(d) * (d.d0 / d.u0)


def reynolds(d: Drop) -> float:
    """Re = rho U0 D0 / eta_0, on the unthickened (f = 0) branch."""
    m = material(d.phi)
    return m.rho * d.u0 * d.d0 / m.eta_0


def weber(d: Drop) -> float:
    """We = rho U0^2 D0 / gamma."""
    m = material(d.phi)
    return m.rho * d.u0**2 * d.d0 / m.gamma


def stokes_particle(d: Drop, rate=shear_rate_naive) -> float:
    """Particle Stokes number St = rho_p a^2 gammadot / eta_f.

    Whether a particle follows the suspending fluid over one unit of strain.
    """
    m = material(d.phi)
    return m.rho_p * m.a**2 * rate(d) / m.eta_s


def peclet_particle(d: Drop, rate=shear_rate_naive) -> float:
    """Pe = 6 pi eta_f a^3 gammadot / (k_B T). Brownian versus advective."""
    m = material(d.phi)
    return 6.0 * np.pi * m.eta_s * m.a**3 * rate(d) / (K_B * T_LAB)



# --------------------------------------------------------------------------
# candidates added 2026-08-30 after the first audit, because nothing on the
# original list could carry BOTH the velocity and the salt axis at fixed phi,
# and those are the only two knobs that move inside one material.
# --------------------------------------------------------------------------
# At nominal 0 mM the real screening is set by dissolved CO2 and the silica's
# own counterions. debye_length() refuses c -> 0 rather than returning infinity,
# so a residual is supplied explicitly and flagged. 0.01 mM is the value that
# module's own docstring names.
RESIDUAL_IONIC_MM = 0.01


def debye_over_radius(d: Drop) -> float:
    """kappa^-1 / a, the DLVO form of the salt axis. Dimensionless.

    The salt column is a label until it is expressed as a screening length. This
    is the cheapest way to do that without a zeta measurement, and it is exact
    given the ionic strength -- but it assumes the nominal LiCl concentration IS
    the ionic strength, which at 0 mM it is not.
    """
    m = material(d.phi)
    c = d.salt_mM if d.salt_mM > 0 else RESIDUAL_IONIC_MM
    return float(debye_length(c) / m.a)


def pore_deborah(d: Drop, length: str = "contact") -> float:
    """Pore-pressure Deborah number: drainage time over the impact time.

    De_p = tau_p / (D0/U0), with tau_p = L^2/D_p the time for pore pressure to
    relax over the drainage length L. De_p >> 1 is undrained -- the pore fluid
    carries the load and the particle network does not feel it.

    It is computed at the contact-patch length, which ``drainage_lengths``
    documents as the one governing whether the spreading lamella can drain.

    **It does NOT carry the salt axis, and an earlier version of this docstring
    claimed it did.** The permeability used is Carman-Kozeny, which is purely
    geometric: ``kappa ~ a^2 (1-phi)^3 / phi^2``. Salt would enter through an
    electroviscous correction or through the drained modulus ``K_p``, and
    neither is included. What this group actually carries is ``U0`` and ``phi``.
    """
    m = material(d.phi)
    L = drainage_lengths(d.d0, m.a, d.phi)[length]
    # Carman-Kozeny permeability, the standard estimate for a packed bed.
    kappa = m.a**2 * (1.0 - d.phi) ** 3 / (45.0 * d.phi**2)
    # Skeleton-dominated pore diffusivity, D_p = kappa K_p / mu_f. K_p is the
    # network's drained bulk modulus; the closure does not carry one, so the
    # osmotic scale kT/a^3 is used and the ABSOLUTE value of De_p is therefore
    # order-of-magnitude only. Its VARIATION across the dataset, which is all
    # the audit uses, is unaffected because the same K_p is used throughout.
    K_p = K_B * T_LAB / m.a**3
    D_p = kappa * K_p / 1.0e-3
    tau_p = L**2 / D_p
    return float(tau_p / (d.d0 / d.u0))


# The frictional-network formation strain, measured from Lin et al.'s shear
# reversal transients by dst_fit's gate G11: 0.54 to 1.59 on 2 um silica and
# 0.71 to 1.64 on 1.4 um PMMA, over three estimators. Every one of 21 estimates
# exceeds the 0.2 the Basilisk solver borrows from Richards' corn starch.
GAMMA_0_RANGE = (0.54, 1.64)
GAMMA_0_DEFAULT = 1.0


def frictional_fraction_transient(d: Drop, gamma_0: float = GAMMA_0_DEFAULT) -> float:
    """f with the strain clock, not at steady state.

    The steady closure sets ``f = f_ss(sigma)`` instantly. Richards et al.'s
    Eq. 4 says it does not: ``df/dt = -(gammadot/gamma_0)(f - f_ss)``, so f
    relaxes toward its target over an accumulated strain ``gamma_0``. Starting
    from an unsheared drop, after accumulating strain ``gamma_acc``:

        f = f_ss(sigma) * (1 - exp(-gamma_acc / gamma_0))

    This is the one candidate that puts a velocity dependence into the
    activation half through a route the steady audit could not see: the strain
    accumulated during an impact depends on U0 through the lamella length, so a
    lagging f does too.

    ``gamma_acc`` uses the Mobaseri lamella convention, the only one of the four
    that varies with velocity at all -- on ``U0/D0`` the accumulated strain is
    identically 1 and this reduces to a constant multiple of ``f_ss``.
    """
    f_ss = frictional_fraction_shear(d)
    gamma_acc = strain_in_one_impact(d, rate=shear_rate_lamella)
    if not np.isfinite(gamma_acc):
        return float("nan")
    return float(f_ss * (1.0 - np.exp(-gamma_acc / gamma_0)))


CLOCKS = {
    "f_shear (deviatoric gauge)": frictional_fraction_shear,
    "f_transient (strain clock, gamma_0=0.54)":
        lambda d: frictional_fraction_transient(d, gamma_0=GAMMA_0_RANGE[0]),
    "f_transient (strain clock, gamma_0=1.64)":
        lambda d: frictional_fraction_transient(d, gamma_0=GAMMA_0_RANGE[1]),
    "f_shear (deviatoric, sigma* x7.53 for salt)":
        lambda d: frictional_fraction_shear(d, sigma_star_factor=SIGMA_STAR_SALT_SPAN),
    "kappa^-1 / a (Debye)": debye_over_radius,
    # Only the contact length is registered. The drop length differs from it by
    # a constant factor of (1/0.49)^2, so it is the same clock rescaled and
    # returns identical rank statistics; listing both would double-count it.
    "pore Deborah (contact length)": pore_deborah,
    "<tau>/sigma* (deviatoric gauge)": shear_over_sigma_star,
    "f_impact (k=2)": lambda d: frictional_fraction_at_impact(d, k=K_ADDED_MASS_LO),
    "f_impact (k=12)": lambda d: frictional_fraction_at_impact(d, k=K_ADDED_MASS_HI),
    "f_impact (k=1/2, no added mass)": lambda d: frictional_fraction_at_impact(d, k=0.5),
    "f_impact (k=2, sigma* x7.53 for salt)":
        lambda d: frictional_fraction_at_impact(d, k=K_ADDED_MASS_LO,
                                                sigma_star_factor=SIGMA_STAR_SALT_SPAN),
    "Sigma/sigma* (k=2)": lambda d: sigma_over_sigma_star(d, k=K_ADDED_MASS_LO),
    "jamming deficit (phi0-phi)/(phi0-phiM)": jamming_deficit,
    "gammadot D0 (convention A) [1/s]": shear_rate_naive,
    "gammadot Lc (Mobaseri) [1/s]": shear_rate_lamella,
    "strain in one impact, convention A": strain_in_one_impact,
    "strain in one impact, Lc": lambda d: strain_in_one_impact(d, rate=shear_rate_lamella),
    "Re": reynolds,
    "We": weber,
    "St (particle)": stokes_particle,
    "Pe (particle)": peclet_particle,
}
