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

    shear channel:  DeltaU_shear = |dU/dz| * d
    gust  channel:  DeltaU_gust  ~ sigma_u * (d / L_int)

The gust expression is the standard eddy-gradient estimate: an eddy of integral scale
`L_int` carries a velocity contrast of order `sigma_u` across its own width, so across a
body of size `d << L_int` the contrast it presents scales as `sigma_u * d / L_int`. This
is a SCALING, not an identity -- the prefactor is order one and unmeasured here, which
is why the verdict below is quoted as an order of magnitude and not to two figures.

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
    print()
    print("%6s %10s %10s %8s %14s %14s %10s %10s"
          % ("theta", "dU/dz", "sigma_u", "L_int", "dU_shear", "dU_gust",
             "ratio", "OLD"))
    print("%6s %10s %10s %8s %14s %14s %10s %10s"
          % ("(deg)", "(1/s)", "(m/s)", "(m)", "(m/s)", "(m/s)", "new", "sig/dU"))

    ratios, old_ratios = [], []
    for label, key in HEADINGS:
        shear = m["VertShearBand" + key]
        sigma = m["VertSigmaBand" + key]
        lint = m["VertLintBand" + key]

        du_shear = abs(shear) * D_DISC_M
        du_gust = sigma * D_DISC_M / lint
        if du_shear <= 0.0:
            print("%6s  shear is exactly zero; ratio undefined" % label)
            return 1
        ratio = du_gust / du_shear
        old = sigma / du_shear
        ratios.append(ratio)
        old_ratios.append(old)
        print("%6s %10.4f %10.3f %8.2f %14.6f %14.6f %10.2f %10.0f"
              % (label, shear, sigma, lint, du_shear, du_gust, ratio, old))

    lo, hi = min(ratios), max(ratios)
    olo, ohi = min(old_ratios), max(old_ratios)
    print()
    print("VERDICT")
    print("-------")
    print("Moment-consistent ratio: %.1f to %.1f x  (was reported as %.0f to %.0f x)."
          % (lo, hi, olo, ohi))
    print("The overstatement factor is d/L_int = %.4f to %.4f, i.e. %.0f-%.0f x too big."
          % (D_DISC_M / max(m["VertLintBand" + k] for _, k in HEADINGS),
             D_DISC_M / min(m["VertLintBand" + k] for _, k in HEADINGS),
             1.0 / (D_DISC_M / min(m["VertLintBand" + k] for _, k in HEADINGS)),
             1.0 / (D_DISC_M / max(m["VertLintBand" + k] for _, k in HEADINGS))))
    print()
    if lo >= 10.0:
        print("=> The null SURVIVES with room to spare: gust dominates shear by >10x")
        print("   at every heading.")
    elif lo >= 3.0:
        print("=> The null survives but is no longer overwhelming.")
    elif lo > 1.0:
        print("=> The null survives ONLY IN ORDER OF MAGNITUDE. Gust still exceeds shear")
        print("   at every heading, but by a single-figure factor -- and the gust")
        print("   estimate carries an unmeasured order-one prefactor, so the margin at")
        print("   the weakest heading (%.1fx) is NOT resolved against that prefactor."
              % lo)
        print("   'Shear is not the hazard' should be stated as 'shear is the smaller of")
        print("   the two and not the place to look first', not as a 100-fold null.")
        print()
        print("   This also SHARPENS the existing caveat 1: the document already says")
        print("   that in open fetch dU/dz is an order of magnitude larger and the two")
        print("   'become comparable'. At %.1fx sheltered, one order of magnitude in" % lo)
        print("   shear puts shear AHEAD of gust in open fetch, not merely level.")
    else:
        print("=> The null FAILS: shear equals or exceeds gust at some heading.")
    print()
    print("Caveat 2 of the original section is UNAFFECTED and becomes relatively more")
    print("important: the trajectory's own metre-scale height excursion samples a much")
    print("larger mean-wind difference than the disc's 0.269 m span, and that term was")
    print("never in this table.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
