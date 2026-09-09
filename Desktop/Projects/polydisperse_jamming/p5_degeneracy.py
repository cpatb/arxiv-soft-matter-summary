"""p5_degeneracy.py — does the P5 shape effect separate the mechanism from its rival?

THE QUESTION (NEXT_STEPS_2026-09-09.md section 2, item 2)
--------------------------------------------------------
P5 measured that the CV effect on `phi_c` is a distribution SHAPE effect, not a moment
effect: the zero-skew UNIFORM family gives a much smaller `dphi_c` than the LOGNORMAL
one at matched `CV` and matched `d32`. The handoff asks whether that breaks the standing
degeneracy with the simpler rival:

    RIVAL (plain suspension rheology): widening the size distribution raises the maximum
    flowable packing fraction `phi_m`, so `eta ~ (phi_m - phi)^-2` falls, the drop is
    less viscous and it spreads. No dilatancy, no pore pressure, no critical state.

    MECHANISM: widening raises the CRITICAL-STATE fraction `phi_c`, so `phi - phi_c`
    changes sign, the pore pressure changes sign, and the drop liquefies.

The handoff's hope was: "if `phi_m` and `phi_c` respond DIFFERENTLY to skewness, P5
becomes a simulation-only discriminator". This script checks that, and it needs no new
simulation, because the rival's parameter has already been measured here.

WHY IT NEEDS NO NEW RUN
-----------------------
`phi_m` in the rival is the divergence point of the viscosity, i.e. a jamming fraction of
the same packings. The project has measured one such fraction for BOTH families already:
`phi_J`, from `phase2_phiJ.py` on an IMPOSED-PRESSURE ladder, in
`out/phase2_summary.csv`. So a candidate prediction for the P5 contrast is
`dphi_J(uniform)` against `dphi_J(lognormal)`, which is on disk. If that ratio is close
to the `phi_c` ratio P5 measured, the rival predicts P5 and P5 discriminates nothing.

⚠ **THE IDENTIFICATION `phi_m = phi_J` IS AN ASSUMPTION, AND IT IS NOT VERIFIED HERE.**
Raised by a referee 2026-09-09. In Wyart-Cates-type treatments the viscosity-divergence
fraction is STRESS-DEPENDENT, sliding between a frictionless and a frictional limit, and
is closer in kind to a shear-protocol quantity than to an isotropic jamming point.
[UNVERIFIED -- not resolved to a source in this session.] So this script cannot show
that the rival predicts P5; what it shows is weaker and still decisive for the
programme:

    the shape effect is present in a packing fraction measured WITHOUT any shear,
    on a different protocol, by a different code path,

so a discriminator built on P5 would have to argue that the rival's `phi_m` responds to
skewness DIFFERENTLY from every packing-fraction measure this project has (isotropic
`phi_J`, the dense-end E1, the pressure root E2, the isostatic crossing E4). Nothing
supports that, and nobody has shown it. **P5 therefore does not ESTABLISH a
discrimination**; the burden is on any claim that it does.

⚠ Read this beside PREREG_ROOT_PINNING_2026-09-09.md. If the `phi_c` estimator turns out
to be pinned to a ladder centred on `phi_J`, then the two ratios agreeing is guaranteed
by construction and this script's conclusion follows a fortiori.
"""

from __future__ import annotations

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(HERE, "out")

MU_P = 0.3          # the production friction
N_REQ = 1500        # the production system size
CV_LO, CV_HI = 0.0, 0.5
FAMILIES = ("lognormal", "uniform")

# The phi_c contrasts P5 measured, read from the committed verdict files rather than
# typed: see _read_dphi_c below. The pre-registered separation bar was 1.5x.
P5_BAR = 1.5


def _phi_J(family: str, cv: float) -> tuple[float, float]:
    path = os.path.join(OUT_DIR, "phase2_summary.csv")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}")
    with open(path, encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            if (row["family"] == family
                    and abs(float(row["cv_nom"]) - cv) < 1e-9
                    and abs(float(row["mu_p"]) - MU_P) < 1e-9
                    and int(row["n_requested"]) == N_REQ):
                return float(row["phi_J"]), float(row["phi_J_sd"])
    raise SystemExit(f"no phase2 row for {family} cv={cv} mu={MU_P} N={N_REQ}")


def _read_dphi_c(verdict: str) -> float:
    """The ungated E2 `delta(cv 0.00->0.50)` printed in a phase-8 verdict file."""
    path = os.path.join(OUT_DIR, verdict)
    if not os.path.exists(path):
        raise SystemExit(f"missing {path} -- run analyse_phase8.py first")
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    # both the old and the new (gate-set) layouts print a line containing
    # "d(cv 0.00->0.50) = " or "delta(cv 0.00->0.50) = " under the E2 heading
    seen_e2 = False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("E2:") or s == "E2:":
            seen_e2 = True
        elif seen_e2 and ("E3" in s[:4] or "dense" in s[:6]):
            break
        if seen_e2 and "(cv 0.00->0.50)" in s:
            if s.lstrip().startswith(("quasistatic", "steady", "both")):
                continue          # the UNGATED number is the like-for-like one
            head = s.split("=")[1].strip().split()[0]
            return float(head)
    raise SystemExit(f"could not find an E2 cv 0->0.5 contrast in {verdict}")


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    p("=== a SHEAR-FREE packing fraction: phi_J, imposed-pressure ladder, phase 2 ===")
    d_phi_J = {}
    for fam in FAMILIES:
        lo, lo_sd = _phi_J(fam, CV_LO)
        hi, hi_sd = _phi_J(fam, CV_HI)
        d_phi_J[fam] = hi - lo
        p(f"   {fam:10s} cv{CV_LO:.2f}={lo:.6f}+-{lo_sd:.6f}  "
          f"cv{CV_HI:.2f}={hi:.6f}+-{hi_sd:.6f}   dphi_J = {hi - lo:+.6f}")
    ratio_J = d_phi_J["lognormal"] / d_phi_J["uniform"]
    p(f"   its family separation  dphi_J(lognormal)/dphi_J(uniform) = "
      f"{ratio_J:.3f}x")

    p("\n=== the mechanism's parameter: phi_c, constant-volume shear, phase 8 ===")
    d_phi_c = {"lognormal": _read_dphi_c("phase8_verdict.txt"),
               "uniform": _read_dphi_c("p5_verdict.txt")}
    for fam in FAMILIES:
        p(f"   {fam:10s} dphi_c = {d_phi_c[fam]:+.4f}")
    ratio_c = d_phi_c["lognormal"] / d_phi_c["uniform"]
    p(f"   MECHANISM's family separation dphi_c(lognormal)/dphi_c(uniform) = "
      f"{ratio_c:.3f}x   (pre-registered bar {P5_BAR}x)")

    p("\n=== do the two respond DIFFERENTLY to skewness? ===")
    p(f"   phi_J, no shear   {ratio_J:.3f}x")
    p(f"   phi_c, shear      {ratio_c:.3f}x")
    rel = abs(ratio_c - ratio_J) / ratio_J
    p(f"   they differ by {rel * 100:.1f} % of the rival's value")
    # the increments themselves, not just their ratio -- a ratio can agree while the
    # magnitudes do not, and the rival is consumed by the MAGNITUDE too
    for fam in FAMILIES:
        rel_f = abs(d_phi_c[fam] - d_phi_J[fam]) / abs(d_phi_J[fam])
        p(f"   {fam:10s}: dphi_c = {d_phi_c[fam]:+.4f} against dphi_J = "
          f"{d_phi_J[fam]:+.4f}  -- {rel_f * 100:.1f} % apart")

    p("\n=== VERDICT ===")
    if rel < 0.25:
        p("   🔴 P5 DOES NOT ESTABLISH A DISCRIMINATION. A packing fraction measured")
        p("   with NO SHEAR AT ALL -- the isotropic phi_J, a different protocol and a")
        p("   different code path -- shows the same shape effect to within the figure")
        p("   above. So the shape effect is not specific to the critical state; it is a")
        p("   property of the packing geometry. To use P5 against the rheological rival")
        p("   one would have to show that the rival's phi_m responds to skewness")
        p("   DIFFERENTLY from every packing-fraction measure this project has, which")
        p("   nobody has shown. (And note phi_m = phi_J is itself an unverified")
        p("   identification -- see the module docstring.)")
        p("   NEXT_STEPS_2026-09-09 section 2 item 2 is answered NO, and section 1's")
        p("   'no simulation separates them' is if anything strengthened: the shape")
        p("   effect lives in the packing geometry both accounts read.")
        p("   The separators remain the three experimental ones in SYNTHESIS section")
        p("   3(i) -- pore-pressure SIGN, strain HISTORY, and a boundary in (U0, D0).")
    else:
        p("   the two parameters respond differently by more than 25 %; P5 may")
        p("   discriminate. Check the phi_c estimator is not pinned to a ladder")
        p("   centred on phi_J before believing it (PREREG_ROOT_PINNING_2026-09-09.md).")

    path = os.path.join(OUT_DIR, "p5_degeneracy.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
