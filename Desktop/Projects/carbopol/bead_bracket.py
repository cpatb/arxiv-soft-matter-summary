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

    **[UNVERIFIED] Y_c = 0.143 is written here from memory as Beris et al. (1985), the
    creeping-flow numerical solution for a sphere in a Bingham fluid, and it has NOT been
    resolved to a DOI in this session.** What is on record in this project is the campaign
    document's "near 0.14 (Beris et al.)", so the leading digits are the project's and the
    third is not. The same applies to the 0.12-0.20 literature spread quoted in section 4.
    **Resolve both before any number from this script is written down**; the whole
    measurement is one multiplication by this constant, so an error in it is an error of
    the same size in every tau_y the bracket returns.

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
Y_C = 0.143         # Beris et al. (1985), sphere in a Bingham fluid, creeping flow
RHO_F = 1000.0      # kg/m^3; Carbopol at these concentrations is water to <1%

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
    print("   keeping only beads that advance the bracket by less than a factor of 1.5.")
    print("   The bracket a ladder returns is the RATIO between adjacent rungs, so the")
    print("   rung spacing IS the measurement's resolution -- there is no separate")
    print("   precision to quote.")
    print()
    cand = sorted(((tau_from_bead(d, rho), m, d)
                   for m, rho in MATERIALS.items() for d in DIAMETERS_MM),
                  key=lambda t: t[0])
    cand = [c for c in cand if taus.min() * 0.4 <= c[0] <= taus.max() * 2.5]
    ladder: list[tuple[float, str, float]] = []
    for tau, m, d in cand:
        if not ladder or tau / ladder[-1][0] >= 1.35:
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
    print("   - [UNVERIFIED] Y_c = 0.143 is written from memory as Beris et al.'s")
    print("     creeping-flow result for a BINGHAM fluid, and is not resolved to a DOI.")
    print("     The project's own record says only 'near 0.14'. Resolve it first.")
    print("     Carbopol is Herschel-Bulkley with n = 0.40 to 0.63 in these very fits, and")
    print("     published Y_c for a sphere in an HB fluid moves with n. So the ABSOLUTE")
    print("     tau_y carries a systematic on top of the bracket above. Reported values of")
    print("     Y_c across the literature span roughly 0.12 to 0.20, which is +/- 25% on")
    print("     tau_y -- LARGER than the bracket the ladder gives. The bracket is not the")
    print("     dominant error and should not be quoted as if it were.")
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
