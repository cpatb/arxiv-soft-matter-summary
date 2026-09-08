"""shear_vs_gust_moment.py -- the shear-vs-gust null, made moment-consistent.

WHY THIS EXISTS
---------------
FRISBEE_GENERALIZATION.md section 7 reported

    "The gust-induced rolling moment is 37 to 122 times the shear-induced one"

from the ratio `sigma_u / DeltaU`, with `DeltaU = (dU/dz)*d`. A physics referee called
this FATAL and it is (2026-09-04):

    A ROLLING MOMENT is produced by a velocity DIFFERENCE across the span. `DeltaU` is
    such a difference. `sigma_u` is NOT -- it is the RMS of the velocity at a point. A
    spatially uniform gust changes the disc's lift and drag; it exerts no roll couple at
    all. So the table divided a couple by a point velocity and labelled the quotient a
    ratio of moments. The two are not the same kind of object and the 37-122x is not a
    moment ratio.

THE MOMENT-CONSISTENT COMPARISON
--------------------------------
Put both channels on the same footing: the velocity difference each imposes ACROSS THE
DISC DIAMETER `d`.

    shear channel:  DeltaU_shear = |dU/dz| * d          (a MEAN gradient: smooth by
                                                         construction, so linear in d)
    gust  channel:  DeltaU_gust  = f(d)                  (a TURBULENT increment: the
                                                         scaling depends on the regime)

WHICH GUST SCALING -- AND THE FIRST VERSION OF THIS SCRIPT GOT IT WRONG
-----------------------------------------------------------------------
The turbulent velocity increment across a separation `d` has two limits:

    smooth / viscous  (d << eta):   DeltaU ~ sigma_u * (d / L)
    inertial subrange (eta << d << L): DeltaU ~ sqrt(C2) * sigma_u * (d / L)^(1/3)

The first version used the SMOOTH form and defended it in this docstring as "the standard
eddy-gradient estimate", claiming the only unmeasured thing was an order-one prefactor.
The unmeasured thing was the EXPONENT, and it is worth a factor of (L/d)^(2/3) ~ 8.

The regime is decided by a computed number, not by assertion. With eps ~ sigma_u^3 / L
the Kolmogorov scale is eta = (nu^3/eps)^(1/4) = 0.58-0.85 mm against d = 269 mm, so
d/eta = 315-468. The disc is FIRMLY in the inertial subrange and the (d/L)^(1/3) form is
the applicable one. This script prints d/eta beside the answer so the choice can be
checked rather than trusted.

The consequence is that the ORIGINAL 37-122x was nearer the truth than the smooth-limit
"correction" of 1.7-4.9x, even though the original derivation (dividing a couple by a
point RMS) was wrong. Both are reported below so the size of the correction is visible.

`C2` is the Kolmogorov constant for the second-order longitudinal structure function,
`D_LL(r) = C2 (eps r)^(2/3)`. It is ~2.0 (Pope, *Turbulent Flows*, ch. 6; Sreenivasan
1995 reports 2.0 +/- 0.1 across flows). It is a real constant rather than a fitted knob,
but the estimate remains a SCALING -- so the verdict is quoted as an order of magnitude.

Every input is a disc-band quantity read from the generated macro file, so the ratio
stays like-for-like and regenerates with the CFD.

    python shear_vs_gust_moment.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MACROS = HERE / "SeanEarl_DiscFlight" / "P5_Synthesis" / "vertical_generated.tex"

# Potts & Crowther Table 3-1 disc diameter, as used by the flight model.
# P1_AnalyticalModel/coefficients.py:516
D_DISC_M = 0.269

# Kinematic viscosity of air at 15 C, 1 atm (Cengel & Ghajar, Table A-15). Only enters
# through eta ~ nu^(3/4), so a 20 % error in nu moves d/eta by 15 % and cannot change
# the regime verdict at d/eta ~ 300-500.
NU_AIR = 1.5e-5          # m^2/s

# Kolmogorov constant for the second-order longitudinal structure function,
# D_LL(r) = C2 (eps r)^(2/3). Sreenivasan (1995) reports 2.0 +/- 0.1 across flows;
# Pope, Turbulent Flows, ch. 6. NOT a fitted knob.
C2_KOLMOGOROV = 2.0

# Above this d/eta the inertial-subrange form is used. The classical requirement is
# eta << d << L; at the loose end of the literature the inertial range is taken to open
# around d/eta ~ 30-50. Set at 50 (the conservative end) and reported beside the answer,
# because the whole point of the 2026-09-04 retraction is that this choice was made by
# assertion the first time.
D_OVER_ETA_INERTIAL = 50.0

HEADINGS = [("90", "HNineZero"), ("210", "HTwoOneZero"), ("270", "HTwoSevenZero")]


def load_macros() -> dict:
    txt = MACROS.read_text(encoding="utf-8")
    out = {}
    for name, val in re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{([^}]*)\}", txt):
        try:
            out[name] = float(val)
        except ValueError:
            pass          # ranges like "5.7--6.7" are not scalars; skip
    return out


def main() -> int:
    if not MACROS.is_file():
        print("macro file not found:", MACROS)
        return 1
    m = load_macros()

    print("=" * 78)
    print("Shear vs gust as ROLLING-MOMENT drivers: both as a velocity difference")
    print("across the disc diameter d = %.3f m" % D_DISC_M)
    print("=" * 78)
    print("Retired (FRISBEE_GENERALIZATION.md s7, 2026-09-03): sigma_u / DeltaU_shear,")
    print("which divides a point RMS by a span difference. Not a moment ratio.")
    print("Retired (2026-09-04): the SMOOTH-limit gust sigma_u*(d/L), which assumes the")
    print("disc is smaller than the Kolmogorov scale. The d/eta column below is the")
    print("number that decides which form applies, and it is printed so it can be")
    print("checked rather than trusted.")
    print()
    print("%6s %9s %8s %7s %9s %8s %12s %12s"
          % ("theta", "dU/dz", "sigma_u", "L_int", "eta", "d/eta", "dU_shear", "dU_gust"))
    print("%6s %9s %8s %7s %9s %8s %12s %12s"
          % ("(deg)", "(1/s)", "(m/s)", "(m)", "(mm)", "-", "(m/s)", "(m/s)"))

    rows = []
    for label, key in HEADINGS:
        shear = m["VertShearBand" + key]
        sigma = m["VertSigmaBand" + key]
        lint = m["VertLintBand" + key]

        # Regime test. eps ~ sigma_u^3 / L is the standard inertial-range estimate of the
        # dissipation rate; eta = (nu^3/eps)^(1/4) is the Kolmogorov scale.
        eps = sigma ** 3 / lint
        eta = (NU_AIR ** 3 / eps) ** 0.25
        d_over_eta = D_DISC_M / eta

        du_shear = abs(shear) * D_DISC_M
        if du_shear <= 0.0:
            print("%6s  shear is exactly zero; ratio undefined" % label)
            return 1

        du_gust_smooth = sigma * (D_DISC_M / lint)
        du_gust_inertial = (C2_KOLMOGOROV ** 0.5) * sigma * (D_DISC_M / lint) ** (1.0 / 3.0)
        inertial = d_over_eta >= D_OVER_ETA_INERTIAL
        du_gust = du_gust_inertial if inertial else du_gust_smooth

        rows.append({
            "label": label, "shear": shear, "sigma": sigma, "lint": lint,
            "eta": eta, "d_over_eta": d_over_eta, "inertial": inertial,
            "du_shear": du_shear, "du_gust": du_gust,
            "ratio": du_gust / du_shear,
            "ratio_smooth": du_gust_smooth / du_shear,
            "ratio_pointrms": sigma / du_shear,
        })
        print("%6s %9.4f %8.3f %7.2f %9.3f %8.0f %12.6f %12.6f"
              % (label, shear, sigma, lint, eta * 1e3, d_over_eta, du_shear, du_gust))

    print()
    print("%6s %14s %14s %14s" % ("theta", "ratio (used)", "smooth limit", "retired RMS"))
    for r in rows:
        print("%6s %14.1f %14.1f %14.0f"
              % (r["label"], r["ratio"], r["ratio_smooth"], r["ratio_pointrms"]))

    n_inertial = sum(r["inertial"] for r in rows)
    lo = min(r["ratio"] for r in rows)
    hi = max(r["ratio"] for r in rows)
    slo = min(r["ratio_smooth"] for r in rows)
    shi = max(r["ratio_smooth"] for r in rows)
    plo = min(r["ratio_pointrms"] for r in rows)
    phi = max(r["ratio_pointrms"] for r in rows)

    print()
    print("REGIME")
    print("------")
    print("d/eta = %.0f to %.0f, against the inertial-subrange threshold of %.0f."
          % (min(r["d_over_eta"] for r in rows),
             max(r["d_over_eta"] for r in rows), D_OVER_ETA_INERTIAL))
    if n_inertial == len(rows):
        print("All %d headings are FIRMLY inertial-subrange: the (d/L)^(1/3) second-order"
              % len(rows))
        print("structure function applies and the smooth-limit (d/L) form does not.")
    elif n_inertial == 0:
        print("No heading reaches the inertial subrange; the smooth limit is used.")
        print("This has never happened for this body and would be worth checking.")
    else:
        print("MIXED regime on %d of %d headings -- the verdict below mixes two scalings"
              % (n_inertial, len(rows)))
        print("and should not be quoted as a single range.")

    print()
    print("VERDICT")
    print("-------")
    print("Moment-consistent ratio: %.0f to %.0f x." % (lo, hi))
    print("  smooth-limit form (RETIRED, too small by (L/d)^(2/3)): %.1f to %.1f x"
          % (slo, shi))
    print("  point-RMS form    (RETIRED, not a moment ratio):       %.0f to %.0f x"
          % (plo, phi))
    print("The smooth limit understates the inertial one by %.1f to %.1f x."
          % (min(r["ratio"] / r["ratio_smooth"] for r in rows),
             max(r["ratio"] / r["ratio_smooth"] for r in rows)))
    print()
    # The verdict text is SELECTED BY the number above, never written beside it.
    if lo >= 10.0:
        print("=> The null SURVIVES with room to spare: the gust channel exceeds the")
        print("   shear channel by more than 10x at every heading, so the order-one")
        print("   prefactor in the structure-function estimate cannot reverse it.")
        print("   Note the direction of the two retractions: the ORIGINAL 37-122x was")
        print("   nearer the truth than the smooth-limit 'correction' of %.1f-%.1f x,"
              % (slo, shi))
        print("   even though the original derivation was wrong. Being closer to the")
        print("   answer is not the same as being right, and neither is a defence.")
    elif lo >= 3.0:
        print("=> The null survives but is no longer overwhelming.")
    elif lo > 1.0:
        print("=> The null survives ONLY IN ORDER OF MAGNITUDE, and at the weakest")
        print("   heading (%.1fx) it is NOT resolved against the order-one prefactor."
              % lo)
    else:
        print("=> The null FAILS: shear equals or exceeds gust at some heading.")

    print()
    print("STILL TRUE, AND UNAFFECTED BY EITHER RETRACTION")
    print("-----------------------------------------------")
    print("Caveat 1: this is the SHELTERED case. In open fetch dU/dz at 2 m is an order")
    print("of magnitude larger, which at %.0f-%.0f x here leaves gust ahead by only"
          % (lo, hi))
    print("%.1f-%.1f x there -- the null holds in open fetch too, but by a single-figure"
          % (lo / 10.0, hi / 10.0))
    print("margin that the order-one prefactor CAN reach at the weakest heading.")
    print("Caveat 2: the trajectory's own metre-scale height excursion samples a much")
    print("larger mean-wind difference than the disc's %.3f m span, and that term is"
          % D_DISC_M)
    print("still not in this table. It is deterministic, so it is a learnable bias.")
    print("Caveat 3: sqrt(C2) is a real constant (C2 = %.1f), but D_LL = C2 (eps r)^(2/3)"
          % C2_KOLMOGOROV)
    print("is an isotropic-turbulence result and eps ~ sigma_u^3/L is itself a scaling,")
    print("so this remains an ORDER-OF-MAGNITUDE statement.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
