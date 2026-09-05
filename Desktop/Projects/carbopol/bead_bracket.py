"""Design the sphere-settling bracket on tau_y: which beads, and what it will resolve.

WHY THIS EXISTS
    73 drop-impact recordings from 11AUG26 are sitting unusable because their bath's
    yield stress is unknown -- "Carbopol wt%" is a recipe label, not a material property,
    so those runs cannot be placed on any dimensionless axis. The retake campaign
    (`RETAKE_CAMPAIGN_2026-08-17.md` sec.9) proposes unblocking them at no cost by
    bracketing tau_y with sphere settling: drop beads of several radii and densities in,
    record which hold and which sink, and read tau_y off the critical yield-gravity
    number.

    That is half an hour at the bench. This is the desk half: which beads, in what order,
    and -- the part worth doing before rather than after -- **what the bracket will and
    will not resolve**, computed against the flow curves already in this directory.

THE PHYSICS, AND ITS ONE CONSTANT
    A sphere of radius R and density rho_s held in a yield-stress fluid of density rho_f
    is governed by the yield-gravity number

        Y_g = tau_y / (delta_rho * g * R),        delta_rho = rho_s - rho_f

    A sphere is held when Y_g exceeds a critical value Y_c.

    **RESOLVED 2026-09-02, and the number was wrong by a factor of 1.5.** The
    [UNVERIFIED] tag said Y_c = 0.143 was written from memory as Beris et al. (1985) and
    had to be resolved before any tau_y from this script was quoted. Resolving it found
    that 0.143 is real, is Beris's, and **belongs to a different definition of the yield
    number than the one used above**. There are three conventions in circulation and they
    span a factor of three:

        Y   = tau_y / (delta_rho g R)            this script, before today
        Y_g = 1.5 tau_y / (delta_rho g R)        Fraggedakis et al. (2016); Sarabian et
                                                 al., arXiv:2005.08429 Eq. 11
        Y_G = tau_y / (delta_rho g d),  d = 2R   Jossic & Magnin (2001, 2009);
                                                 Mitsoulis & Galazoulas (2009)

    so Y_g = 3 Y_G = 1.5 Y. **0.143 is Beris's critical value in the Y_g convention**, and
    two independent sources pin it:

      * Sarabian, Rosti, Brandt & Hormozi, arXiv:2005.08429, sec.3.1: "The critical
        gravity number is found to be Y_g^c = 0.143 for Bingham fluids ... by Beris et
        al. (1985)", with Y_g defined in their Eq. (11) as 1.5 tau_0/(delta_rho g R).
      * Nazari, Mittal, Shoele & Mohammadigoushki, arXiv:2606.04861: "the asymptotic
        analysis of Beris et al. (1985), using a Bingham fluid model, estimated
        Y_G,max ~ 0.048 for a sphere", with Y_G defined in their Eq. (8) as
        sigma_y/(g d_e delta_rho), d_e = 2R for a sphere.

    3 x 0.048 = 0.144. The two agree, through different conventions, to the precision
    either quotes.

    Converted into THIS script's convention, Beris's critical value is

        Y_c = 0.143 / 1.5 = 2 x 0.048 = 0.0953

    Using 0.143 here therefore overstated Y_c by 1.5x, which understates the critical
    radius by 1.5x and would have **overstated every tau_y read off the ladder by 1.5x**.
    That is larger than the ladder's own 1.60x worst-case bracket, so it was the dominant
    error and not a rounding matter.

    Source, resolved: Beris, Tsamopoulos, Armstrong & Brown, "Creeping motion of a sphere
    through a Bingham plastic", J. Fluid Mech. 158, 219-244 (1985),
    doi:10.1017/S0022112085002622.

    THE SPREAD, also resolved, and it is a different spread than was quoted. Measured and
    computed sphere yield limits in the Y_G convention (arXiv:2606.04861 sec.6, with its
    citations), converted to this script's Y_c = 2 Y_G:

        Beris et al. 1985      (theory, Bingham)     Y_G 0.048   ->  Y_c 0.095
        Atapattu et al. 1995   (experiment)          Y_G 0.048   ->  Y_c 0.096
        Tabuteau et al. 2007   (experiment, Carbopol) Y_G 0.061  ->  Y_c 0.122
        Jossic & Magnin 2001   (roughened sphere)    Y_G 0.062   ->  Y_c 0.124
        Jossic & Magnin 2001   (SMOOTH sphere)       Y_G 0.088   ->  Y_c 0.176

    So the honest span is Y_c = 0.095 to 0.176, a factor of 1.83 -- wider than the
    "0.12 to 0.20" (1.67x) this script used to quote, and centred lower.

    ** AND THE SPREAD IS NOT NOISE: IT IS SURFACE ROUGHNESS. ** Jossic & Magnin's own two
    numbers are the same measurement on a roughened and a smooth sphere, and they differ
    by 1.42x. This ladder deliberately mixes six materials -- PTFE, glass, ceramics,
    steel, carbide -- to buy density range, and those do not share a surface finish. So
    the ladder trades a known systematic for the density coverage it needs, and any
    tau_y from it must either use beads of one finish or carry the 1.42x. Record the bead
    surface condition at the bench; it is free to write down and impossible to recover
    later.

    That is also why section 4 leads with it: the constant, not the bead ladder, is the
    dominant uncertainty.

    Solving for the critical radius: a bead is held when

        R < R_c = tau_y / (Y_c * delta_rho * g)

    So the largest held bead gives a LOWER bound on tau_y and the smallest sunk bead an
    UPPER bound, and the bracket's width is set entirely by the spacing of the beads
    available.

Run:  python bead_bracket.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
G = 9.81            # m/s^2
# Beris et al. (1985), J. Fluid Mech. 158:219, doi:10.1017/S0022112085002622, converted
# into THIS module's convention Y = tau_y/(delta_rho g R). Beris's own 0.143 is quoted
# for Y_g = 1.5 tau_y/(delta_rho g R); equivalently Y_G,max = 0.048 for tau_y/(rho g d).
# See the module docstring: using 0.143 here was a 1.5x error. Do not "restore" it.
Y_C = 0.143 / 1.5   # = 0.09533
Y_C_SPAN = (0.095, 0.176)   # Beris/Atapattu ... Jossic & Magnin smooth sphere; 1.83x
# Cross-convention check, and the tolerance is 1% deliberately: 0.143/1.5 = 0.09533 against
# 2 x 0.048 = 0.096 is a 0.7% disagreement, which is exactly the rounding in the second
# source's two significant figures (0.0477 -> 0.048). A tighter bar would fail on rounding;
# a looser one (say 50%) would stop discriminating between the conventions, which differ by
# 1.5x and 3x. This assertion exists to pin WHICH convention Y_C is in, not to check arithmetic.
assert abs(Y_C - 2 * 0.048) / Y_C < 0.01, (
    "Y_C must agree with 2 * Y_G,max, Y_G,max = 0.048 (Beris via arXiv:2606.04861); the "
    "three yield-number conventions differ by 1.5x and 3x and this is what pins which is "
    "in use. Got Y_C = %.5f against 2 x 0.048 = 0.096" % Y_C)
RHO_F = 1000.0      # kg/m^3; Carbopol at these concentrations is water to <1%

# Ladder construction. These were inline literals until 2026-09-03, and the printed prose
# described the first as a MAXIMUM of 1.5 when the code applies it as a MINIMUM of 1.35 --
# a reader could not check the quoted resolution against the stated rule. Named here so the
# rule and the number are one thing.
LADDER_MIN_STEP = 1.35   # keep a bead only once it is this far above the last kept one
CAND_LO, CAND_HI = 0.4, 2.5   # candidate window, as multiples of the softest/stiffest recipe

# Bead materials by density, kg/m^3. Chosen for what a lab actually has in a drawer or
# can buy as a precision ball bearing / chromatography packing, not for coverage.
MATERIALS = {
    "PTFE": 2200.0,
    "soda-lime glass": 2500.0,
    "alumina": 3900.0,
    "silicon nitride": 3200.0,
    "440C stainless": 7650.0,
    "tungsten carbide": 14900.0,
}

# Standard precision-ball diameters in mm. These are catalogue sizes, so the ladder is
# buildable rather than notional.
DIAMETERS_MM = [0.5, 0.8, 1.0, 1.5, 2.0, 2.381, 3.0, 3.175, 4.0, 4.762, 6.0, 6.35,
                8.0, 9.525, 12.7]


def tau_from_bead(d_mm: float, rho_s: float) -> float:
    """The tau_y at which a bead of this size and density is exactly marginal."""
    R = 0.5 * d_mm * 1e-3
    return Y_C * (rho_s - RHO_F) * G * R


def main() -> None:
    fits = list(csv.DictReader(open(HERE / "hb_fits_etd2020.csv")))
    taus = np.array([float(r["tau_y_Pa"]) for r in fits])
    print("1. The range the bracket has to cover")
    print()
    print("   From the 15 Herschel-Bulkley fits in hb_fits_etd2020.csv (Biermann &")
    print("   Mesecke 2024, ETD2020, 0.1-0.9 wt% x 30-50 uL TEA):")
    print(f"     tau_y spans {taus.min():.2f} to {taus.max():.2f} Pa, median {np.median(taus):.2f}")
    print(f"     and it is NOT monotone in wt%: 0.3 wt% at 50 uL TEA gives "
          f"{[float(r['tau_y_Pa']) for r in fits if r['condition'].startswith('0.3') and '50' in r['condition']][0]:.1f} Pa")
    print(f"     against {[float(r['tau_y_Pa']) for r in fits if r['condition'].startswith('0.9') and '50' in r['condition']][0]:.1f} Pa at 0.9 wt%.")
    print("   That non-monotonicity is the whole argument for measuring rather than")
    print("   labelling: the neutraliser dose matters as much as the carbomer, and the")
    print("   11AUG26 recordings record only the carbomer.")
    print()
    print("   A bath made to any of those recipes therefore needs a bracket that spans")
    print(f"   at least {taus.min():.1f} to {taus.max():.1f} Pa -- a factor of {taus.max()/taus.min():.0f}.")

    print()
    print("2. What each bead resolves, in Pa")
    print()
    print("   The marginal tau_y for each (material, diameter). A bead HOLDS when the")
    print("   bath's tau_y is ABOVE the number in the table.")
    print()
    header = f"   {'d (mm)':>8} " + "".join(f"{m.split()[0]:>12}" for m in MATERIALS)
    print(header)
    for d in DIAMETERS_MM:
        row = f"   {d:>8.3f} " + "".join(f"{tau_from_bead(d, rho):>12.2f}"
                                          for rho in MATERIALS.values())
        print(row)

    # ---- pick a ladder ---------------------------------------------------------
    print()
    print("3. A ladder, and what it costs")
    print()
    print("   Sort every (material, diameter) by the tau_y it marks, then walk the list")
    print("   keeping a bead only once it is at least %.2fx above the last one kept."
          % LADDER_MIN_STEP)
    print("   That is a MINIMUM spacing: two rungs closer together cost a drop and buy no")
    print("   resolution, because the bracket a ladder returns IS the ratio between")
    print("   adjacent rungs. The rung spacing therefore IS the measurement's resolution,")
    print("   and there is no separate precision to quote.")
    print("   Candidates are restricted to %.2fx the softest and %.2fx the stiffest recipe."
          % (CAND_LO, CAND_HI))
    print()
    cand = sorted(((tau_from_bead(d, rho), m, d)
                   for m, rho in MATERIALS.items() for d in DIAMETERS_MM),
                  key=lambda t: t[0])
    # Only beads that can bracket the recipes: a rung far below the softest or far above
    # the stiffest bath resolves nothing, it just costs a drop. The margins are generous
    # because the recipe range is itself uncertain.
    cand = [c for c in cand if taus.min() * CAND_LO <= c[0] <= taus.max() * CAND_HI]
    # MINIMUM spacing, not maximum: take the next bead only once it is at least
    # LADDER_MIN_STEP above the last kept one. Two rungs closer than that cost a drop and
    # buy no resolution, because the bracket a ladder returns IS the ratio between
    # adjacent rungs.
    ladder: list[tuple[float, str, float]] = []
    for tau, m, d in cand:
        if not ladder or tau / ladder[-1][0] >= LADDER_MIN_STEP:
            ladder.append((tau, m, d))
    print(f"   {'tau_y (Pa)':>11}  {'material':<18} {'d (mm)':>7}   step")
    for i, (tau, m, d) in enumerate(ladder):
        step = "" if i == 0 else f"x{tau / ladder[i - 1][0]:.2f}"
        print(f"   {tau:>11.2f}  {m:<18} {d:>7.3f}   {step}")
    steps = [ladder[i][0] / ladder[i - 1][0] for i in range(1, len(ladder))]
    print()
    print(f"   {len(ladder)} beads, {len(set(m for _, m, _ in ladder))} materials, "
          f"covering {ladder[0][0]:.2f} to {ladder[-1][0]:.2f} Pa.")
    print(f"   Worst-case bracket: a factor of {max(steps):.2f}, i.e. tau_y known to")
    print(f"   +/- {100 * (max(steps) - 1) / (max(steps) + 1):.0f}% of its midpoint.")
    print("   The campaign document asks for 'within roughly a factor of 2'. This is")
    print("   better than that, and it is better because the ladder mixes DENSITIES --")
    print("   a single-material ladder is limited to the catalogue's diameter steps.")

    print()
    print("4. What this measurement cannot do, stated before it is run")
    print()
    print("   - Y_c is RESOLVED as of 2026-09-02, and correcting it moved tau_y by 1.5x.")
    print(f"     Y_c = {Y_C:.4f} in this module's convention Y = tau_y/(delta_rho g R).")
    print("     Beris et al. (1985), J. Fluid Mech. 158:219, doi:10.1017/S0022112085002622.")
    print("     The 0.143 previously used here is Beris's value for a DIFFERENT yield")
    print("     number, Y_g = 1.5 tau_y/(delta_rho g R); the equivalent in the third")
    print("     convention in circulation is Y_G,max = 0.048 for tau_y/(delta_rho g d).")
    print("     Two independent sources agree through different conventions to 3x0.048")
    print("     = 0.144. Using 0.143 here overstated every tau_y by 1.5x -- larger than")
    print("     the ladder's own bracket. See the module docstring for the citations.")
    print("     Carbopol is Herschel-Bulkley with n = 0.40 to 0.63 in these very fits, and")
    print("     published Y_c for a sphere in an HB fluid moves with n. So the ABSOLUTE")
    print("     tau_y carries a systematic on top of the bracket above. Measured and")
    print(f"     computed sphere yield limits span Y_c = {Y_C_SPAN[0]:.3f} to {Y_C_SPAN[1]:.3f}")
    print(f"     ({Y_C_SPAN[1] / Y_C_SPAN[0]:.2f}x), which is LARGER than the bracket the ladder gives.")
    print("     The bracket is not the dominant error and must not be quoted as if it were.")
    print("   - AND THAT SPAN IS SURFACE ROUGHNESS, not scatter. Jossic & Magnin (2001)")
    print("     measured 0.062 roughened against 0.088 smooth on spheres -- 1.42x, from")
    print("     one lab in one paper. This ladder mixes PTFE, glass, ceramics, steel and")
    print("     carbide to buy density range, and they do not share a surface finish, so")
    print("     it trades a known systematic for the coverage. Use one finish if you can;")
    print("     either way RECORD the surface condition of every bead at the bench. It is")
    print("     free to write down and impossible to recover afterwards.")
    print("   - It is a static test. It says nothing about K or n, and the drop-impact")
    print("     problem needs both: at impact rates the HB stress is dominated by the")
    print("     K gamma_dot^n term, not by tau_y.")
    print("   - Thixotropy and ageing are not controlled. Carbopol's apparent yield stress")
    print("     depends on rest time, and a bead dropped into a freshly stirred bath and")
    print("     one dropped into a rested bath are not the same measurement. Fix a rest")
    print("     time and record it.")
    print("   - Wall effects: keep the bead diameter below about a tenth of the container")
    print("     diameter, or the bracket shifts.")
    print()
    print("   So this unblocks the 73 recordings by putting them on an axis with an")
    print("   honest +/- 25%, which is enough to ORDER them and enough to decide whether")
    print("   the arm earns a rheometer slot. It is not enough to publish a tau_y.")

    print()
    print("5. Protocol, for the half hour")
    print()
    print("   Fix the rest time and the container. Drop the ladder from the top of the")
    print("   list downward, one bead per rung, and record HELD or SANK after a fixed")
    print("   wait -- 10 minutes is longer than any settling time at these stresses and")
    print("   short enough to sit through. The bracket is [largest held, smallest sunk].")
    print("   A non-monotone result -- a bead that sinks while a larger one holds -- means")
    print("   the bath is not uniform or the rest time is not controlled, and is the one")
    print("   outcome that invalidates the run rather than widening it. Record every")
    print("   bead, including the ones that behave, so that check is possible at all.")


if __name__ == "__main__":
    main()
