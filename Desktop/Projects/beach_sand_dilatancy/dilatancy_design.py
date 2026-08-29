"""
dilatancy_design.py -- the calculations that decide whether the beach-sand protocol is
worth running, done before any sand is collected.

PROTOCOL.md poses four questions and estimates two of them in a table.  This module turns
those estimates into computable, testable functions and answers the design question each
one implies:

  Q1  Which grain diameter enters the Darcy-Reynolds law for polydisperse sand?
      DESIGN QUESTION: how well-sorted may the sand be before the candidate diameters
      (d50, d10, Sauter) stop making distinguishable predictions?
      ANSWER: they are distinguishable everywhere realistic.  Even at C_u = 2, which is
      a well-sorted beach sand, the d10 and d50 forms of a_max differ by ~60%, far above
      any plausible measurement scatter.  Q1 is well posed.

  Q4  Where does the V^1.6 law have to saturate, because suction cannot exceed about one
      atmosphere?
      DESIGN QUESTION: is that ceiling inside a reachable drop-tower speed range?
      ANSWER: for 0.2 mm sand it lands at 1.96 m/s, i.e. inside Strader's own
      1-2.5 m/s window, which is why their theory may already be brushing it.  The
      ceiling scales as d^2, so a coarse-sand control is unreachable and a fine-sand
      case is reachable with a 0.19 m drop.  That contrast is the cleanest possible test.
      PROTOCOL.md §7.4 now treats this as a MEASUREMENT of psi rather than a yes/no check,
      since psi is the unknown and the departure speed determines it.

Everything is derived from two textbook relations and one dimensional argument:

  Kozeny-Carman   k = n^3 d_32^2 / (180 (1-n)^2),  n = porosity
  Darcy-Reynolds  dp = eta psi V D / k
  cavitation      dp cannot exceed p_atm - p_vapour

The weakest number in the chain is the dilation-angle factor psi ~ 0.1, and every
function that uses it takes it as an argument so its influence can be swept rather than
buried.  PROTOCOL.md says the same thing; this module makes it operational.

SI units throughout.  Run `python dilatancy_design.py` for the report and figures.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

# --------------------------------------------------------------------------- #
# Constants, with provenance.
# --------------------------------------------------------------------------- #
RHO_W = 1000.0          # kg/m^3   pore water
RHO_S = 2650.0          # kg/m^3   quartz sand grains
ETA_W = 1.0e-3          # Pa.s     water at 20 C
P_ATM = 101325.0        # Pa
P_VAP = 2340.0          # Pa       water vapour pressure at 20 C
G = 9.81

# NOTATION WARNING, and it cost a factor of 7.6 while writing this module.
# PROTOCOL.md uses "phi" for two different quantities: the POROSITY in the
# Kozeny-Carman line of §1 ("k ~ 1e-3 d^2, phi ~ 0.4") and the PACKING (solid) FRACTION
# everywhere else ("phi_0 = 0.560 vs 0.604").  A natural sand pack has porosity ~0.40
# and solid fraction ~0.60, so reading one for the other changes k by (0.6/0.4)^3 *
# (0.6/0.4)^2 = 7.6x and moves every derived speed with it.
# This module therefore never uses the bare name "phi".  It takes POROSITY explicitly.
POROSITY_DEFAULT = 0.40   # void fraction of a natural sand pack; solid fraction ~0.60
PSI_DEFAULT = 0.10        # dilation-angle factor; the weakest number in the chain
KC_CONST = 180.0          # Kozeny-Carman constant

# Strader, Causley, Dijksman & Clark, Phys. Rev. Fluids 7, 054304 (2022):
# a_max ~ d^-0.8 eta^0.4 V^1.6 for monodisperse beads at phi0 ~ 0.59.
STRADER_EXP = dict(d=-0.8, eta=0.4, V=1.6)


def _check_positive(**kw):
    for name, value in kw.items():
        arr = np.asarray(value, dtype=float)
        if not np.all(np.isfinite(arr)) or np.any(arr <= 0.0):
            raise ValueError(f"{name} must be finite and positive, got {value!r}")


# --------------------------------------------------------------------------- #
# 1. Grain-size distribution descriptors.
# --------------------------------------------------------------------------- #
def lognormal_quantile(d50, sigma_ln, p):
    """Grain diameter at mass-fraction-finer p, for a lognormal mass distribution."""
    _check_positive(d50=d50)
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie strictly inside (0, 1)")
    return float(d50 * np.exp(norm.ppf(p) * sigma_ln))


def sauter_diameter(d50, sigma_ln):
    """Sauter (surface-area) mean diameter d_32 of a lognormal MASS distribution.

    The specific surface of a volume distribution f(d) is S = int 6 f(d)/d dd, and
    d_32 = 6/S = 1 / E_vol[1/d].  For a lognormal volume distribution with median d50 and
    log-standard-deviation s, E_vol[1/d] = exp(s^2/2)/d50, so

        d_32 = d50 * exp(-s^2 / 2).

    This is the diameter Kozeny-Carman actually wants, because permeability is set by
    surface area per unit volume rather than by a length percentile.
    """
    _check_positive(d50=d50)
    return float(d50 * np.exp(-0.5 * sigma_ln**2))


def uniformity_coefficient(sigma_ln):
    """C_u = d60/d10 for a lognormal distribution.

    d60/d10 = exp((z60 - z10) s) with z60 = 0.2533, z10 = -1.2816, so
    C_u = exp(1.5349 s).  Inverting gives sigma_ln from a measured C_u.
    """
    return float(np.exp((norm.ppf(0.60) - norm.ppf(0.10)) * sigma_ln))


def sigma_ln_from_Cu(Cu):
    """Inverse of uniformity_coefficient.  C_u = 1 is monodisperse."""
    if Cu < 1.0:
        raise ValueError("C_u cannot be below 1")
    return float(np.log(Cu) / (norm.ppf(0.60) - norm.ppf(0.10)))


def size_descriptors(d50, Cu):
    """All three candidate diameters for a sand of median d50 and sorting C_u."""
    s = sigma_ln_from_Cu(Cu)
    return dict(sigma_ln=s, Cu=float(Cu), d50=float(d50),
                d10=lognormal_quantile(d50, s, 0.10),
                d60=lognormal_quantile(d50, s, 0.60),
                d_sauter=sauter_diameter(d50, s))


# --------------------------------------------------------------------------- #
# 2. Permeability and the Darcy-Reynolds suction.
# --------------------------------------------------------------------------- #
def permeability_kozeny_carman(d_eff, porosity=POROSITY_DEFAULT, const=KC_CONST):
    """Intrinsic permeability k (m^2) from Kozeny-Carman.

        k = n^3 d^2 / (const * (1 - n)^2),   n = porosity (void fraction)

    At n = 0.40 this returns k = 9.88e-4 d^2, reproducing the k ~ 1e-3 d^2 shorthand and
    the k = 4e-11 m^2 entry for 0.2 mm sand in PROTOCOL.md §2.  Passing the solid
    fraction here instead of the porosity inflates k by ~7.6x; see the notation warning
    at the top of this file.
    """
    _check_positive(d_eff=d_eff)
    if not 0.0 < porosity < 1.0:
        raise ValueError("porosity (void fraction) must lie in (0, 1)")
    return float(porosity**3 * d_eff**2 / (const * (1.0 - porosity) ** 2))


def darcy_reynolds_suction(V, d_eff, D=0.020, eta=ETA_W, psi=PSI_DEFAULT,
                           porosity=POROSITY_DEFAULT, cap_at_cavitation=True):
    """Pore suction required to feed dilation at rate ~psi V over intruder scale D.

    dp = eta psi V D / k.  Capped at p_atm - p_vapour by default, because suction beyond
    that cavitates or draws air in and the pack cannot sustain it.  The cap is the whole
    content of Q4, so it is on by default and can be switched off to see the uncapped law.
    """
    _check_positive(V=V, d_eff=d_eff, D=D, eta=eta, psi=psi)
    k = permeability_kozeny_carman(d_eff, porosity)
    dp = eta * psi * np.asarray(V, float) * D / k
    if cap_at_cavitation:
        dp = np.minimum(dp, P_ATM - P_VAP)
    return dp


def hardening_ratio(V, d_eff, D=0.020, **kw):
    """dp / (rho_s V^2): suction against the inertial grain stress.

    Above ~1 the pack hardens under the intruder; below ~1 it drains and the intruder
    sinks.  This is the number PROTOCOL.md §2 tabulates.
    """
    dp = darcy_reynolds_suction(V, d_eff, D=D, **kw)
    return dp / (RHO_S * np.asarray(V, float) ** 2)


def cavitation_speed(d_eff, D=0.020, eta=ETA_W, psi=PSI_DEFAULT,
                     porosity=POROSITY_DEFAULT):
    """Impact speed at which the required suction reaches the cavitation ceiling.

    V_cav = k (p_atm - p_vap) / (eta psi D), and k ~ d^2, so V_cav ~ d^2.  Above this
    speed the suction cannot grow, so a_max must depart from Strader's V^1.6.
    """
    k = permeability_kozeny_carman(d_eff, porosity)
    return float(k * (P_ATM - P_VAP) / (eta * psi * D))


# --------------------------------------------------------------------------- #
# 3. Q1: do the candidate diameters make distinguishable predictions?
# --------------------------------------------------------------------------- #
def amax_ratio_between_diameters(d_a, d_b, exponent=STRADER_EXP["d"]):
    """Ratio of predicted peak deceleration when the law is evaluated on d_a vs d_b.

    Strader's a_max ~ d^-0.8, so using a smaller characteristic diameter predicts a
    LARGER a_max.  The ratio is what an experiment has to resolve.
    """
    _check_positive(d_a=d_a, d_b=d_b)
    return float((d_a / d_b) ** exponent)


def q1_discrimination(d50, Cu, exponent=STRADER_EXP["d"]):
    """How far apart are the d50, d10 and Sauter predictions for this sand?"""
    s = size_descriptors(d50, Cu)
    return dict(
        **s,
        ratio_d10_over_d50=amax_ratio_between_diameters(s["d10"], s["d50"], exponent),
        ratio_sauter_over_d50=amax_ratio_between_diameters(s["d_sauter"], s["d50"],
                                                           exponent),
        ratio_d10_over_sauter=amax_ratio_between_diameters(s["d10"], s["d_sauter"],
                                                           exponent))


def psiD_upper_bound(d_eff, V_max, eta=ETA_W, porosity=POROSITY_DEFAULT):
    """Upper bound on the product psi*D implied by a clean, unsaturated V^n scaling.

    If an experiment reports an unbroken power law up to V_max on grains of size d_eff,
    the required suction must have stayed below the cavitation ceiling over that whole
    range.  Setting dp = eta psi V_max D / k <= p_atm - p_vap gives

        psi * D  <=  k (p_atm - p_vap) / (eta V_max).

    The bound is on the PRODUCT because the intruder scale D and the dilation factor psi
    are degenerate in the Darcy-Reynolds estimate; quoting a psi alone requires stating
    the D it assumes.
    """
    _check_positive(d_eff=d_eff, V_max=V_max, eta=eta)
    k = permeability_kozeny_carman(d_eff, porosity)
    return float(k * (P_ATM - P_VAP) / (eta * V_max))


def psi_bound_from_strader(d_fine=53e-6, V_max=2.5, D_assumed=0.020, **kw):
    """What does Strader et al.'s own clean scaling imply about psi?

    They report a_max ~ d^-0.8 eta^0.4 V^1.6 holding up to V ~ 2.5 m/s with beads down to
    53 um.  If psi were 0.1, the suction demanded at their finest size and highest speed
    would exceed one atmosphere by more than an order of magnitude, and the law would have
    had to break.  It did not, so psi must be smaller.

    This is the cheapest available constraint on the loosest number in PROTOCOL.md, and it
    comes from the literature rather than from new sand.
    """
    psiD = psiD_upper_bound(d_fine, V_max, **kw)
    return dict(d_fine=d_fine, V_max=V_max, D_assumed=D_assumed,
                psiD_max=psiD, psi_max=float(psiD / D_assumed),
                factor_below_protocol=float(PSI_DEFAULT / (psiD / D_assumed)))


def psiD_lower_bound(d_eff, V, porosity=POROSITY_DEFAULT, eta=ETA_W, rho_s=RHO_S):
    """Smallest psi*D for which the suction still beats the inertial grain stress.

    Requiring dp >= rho_s V^2, i.e. hardening rather than draining, gives

        psi * D  >=  rho_s V k / eta.

    If a published law is dilatancy-controlled at a given grain size, psi*D must lie
    above this.
    """
    _check_positive(d_eff=d_eff, V=V, eta=eta)
    k = permeability_kozeny_carman(d_eff, porosity)
    return float(rho_s * V * k / eta)


def consistency_window(d_fine=53e-6, d_coarse=850e-6, V_max=2.5, **kw):
    """Is there ANY psi*D consistent with both ends of a grain-size sweep?

    Two requirements pull in opposite directions:
      (a) the fine end must not cavitate, or the power law would break there:
          psi*D <= psiD_upper_bound(d_fine, V_max)
      (b) the coarse end must still be pore-pressure controlled, or the same law could
          not be dilatancy physics there:  psi*D >= psiD_lower_bound(d_coarse, V_max)

    Since k ~ d^2, (a) scales as d_fine^2 and (b) as d_coarse^2, so a wide grain-size
    sweep squeezes the window as (d_fine/d_coarse)^2.  For Strader's 53-850 um range that
    is a factor of 257, and the window turns out to be EMPTY by a factor of 43.0.

    Read this as a statement about the simple estimate, not about their measurements: a
    single constant psi with the intruder diameter as the only length cannot span their
    range.  Something else is operating, and finding out what is a better question than
    the one PROTOCOL.md §2 currently asks.
    """
    hi = psiD_upper_bound(d_fine, V_max, **{k: v for k, v in kw.items()
                                            if k in ("eta", "porosity")})
    lo = psiD_lower_bound(d_coarse, V_max, **{k: v for k, v in kw.items()
                                              if k in ("eta", "porosity", "rho_s")})
    return dict(psiD_max_no_cavitation=hi, psiD_min_for_hardening=lo,
                consistent=bool(lo <= hi),
                shortfall_factor=float(lo / hi) if hi > 0 else np.inf,
                d_fine=d_fine, d_coarse=d_coarse, V_max=V_max)


def match_field_regime(d_eff, D_lab, V_lab, D_field=0.10, V_field=1.0):
    """The lab and the field sit in the same drainage regime when D/V matches.

    The hardening ratio is dp/(rho_s V^2) = eta psi D / (k rho_s V), so at fixed sand it
    depends on the intruder only through D/V.  Matching a footstep therefore requires

        D_lab / V_lab = D_field / V_field,

    which for a 0.1 m foot at 1 m/s means either a 0.2 m intruder at 2 m/s or, far more
    cheaply, the existing 20 mm sphere at 0.20 m/s -- a 2 mm drop.  PROTOCOL.md §5 takes the
    second option; do not buy a 200 mm sphere.  Returns both, as D_lab_required and
    V_lab_required, plus the ratio achieved.

    Caveat: this matches the DRAINAGE group only.  The cavitation demand goes as V*D, which
    the slow-drop fix reduces by 25x relative to a footstep, so a real footstep on fine sand
    is cavitation-limited while the matched lab drop is not.  Matching both groups at once
    forces D_lab = D_field and V_lab = V_field, i.e. full scale.
    """
    _check_positive(D_lab=D_lab, V_lab=V_lab, D_field=D_field, V_field=V_field)
    target = D_field / V_field
    achieved = D_lab / V_lab
    return dict(target_D_over_V=float(target), achieved_D_over_V=float(achieved),
                mismatch_factor=float(target / achieved),
                D_lab_required=float(target * V_lab),
                V_lab_required=float(D_lab / target))


def minimum_Cu_to_discriminate(d50, resolvable_fraction=0.15,
                               exponent=STRADER_EXP["d"], Cu_max=10.0):
    """Smallest sorting coefficient at which the d10 and d50 predictions differ by more
    than `resolvable_fraction` (default 15%, a realistic scatter for a_max).

    This is the design number: below it the experiment cannot tell the hypotheses apart,
    above it it can.
    """
    Cus = np.linspace(1.0 + 1e-6, Cu_max, 4000)
    for Cu in Cus:
        r = q1_discrimination(d50, Cu, exponent)["ratio_d10_over_d50"]
        if abs(r - 1.0) >= resolvable_fraction:
            return float(Cu)
    return float("nan")
