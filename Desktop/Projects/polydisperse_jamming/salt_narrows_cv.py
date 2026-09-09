"""salt_narrows_cv.py — a consequence nobody asked for: the double layer NARROWS the PSD.

THE OBSERVATION
---------------
`phi_axis_offset.py` estimates that a charge-stabilised pack jams when the EFFECTIVE
spheres -- true diameter plus twice the balance separation `h = lambda_D ln(sigma*/sigma)`
-- reach the hard-sphere critical state. That estimate was made to ask whether the
0.117 phi gap is closable.

But adding a CONSTANT shell to every particle does something else, and it is a separate,
testable prediction: it changes the SHAPE of the distribution the mechanics sees.

    d_eff = d + 2h    for every particle, h the same for all of them
    =>  the standard deviation is UNCHANGED   (adding a constant does not move sigma)
        the mean RISES by 2h
    =>  CV_eff = sigma / (d10 + 2h) = CV * d10 / (d10 + 2h)   <  CV

**The double layer makes the pack look NARROWER than it is**, and it does so most at low
salt, where the layer is thickest. Since the whole mechanism runs on `dphi_c/dCV`, this
predicts:

    the CV effect should be WEAKER at low salt and STRONGER at high salt,
    by the factor  d10 / (d10 + 2 lambda_D ln(sigma*/sigma)).

⚠ WHAT THIS IS AND IS NOT. It is a consequence of the effective-diameter ESTIMATE, which
`phi_axis_offset.py` labels a reachability estimate and not a value of `phi_c`. So this is
a prediction of that model, not of the measurement. It is worth writing down because it is
CHEAP TO TEST -- the lab's salt ladder already exists -- and because it points the
OPPOSITE way to the other salt effect in the same model (screening also moves `phi_c^eff`
itself, 0.23 -> 0.56 across 0-100 mM), so the two must be separated rather than summed by
eye.

    python salt_narrows_cv.py
"""

from __future__ import annotations

import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
D10 = 800e-9              # lab silica number mean


def _dphi_c_dcv() -> tuple[float, str]:
    """The measured slope of phi_c against CV, READ from out/phase8_verdict.txt.

    🔴 Hardcoded as 0.1247 until 2026-09-09, which is how a number goes stale: the
    verdict file now prints the slope under FOUR gate sets and they span 0.088-0.125.
    The UNGATED one is returned, because that is the one the rest of this module's
    prose is written against, and the gate is named in the output so the reader can
    see which of the four it is. See analyse_phase8's docstring for why no single
    gate set is the right answer -- and PREREG_ROOT_PINNING_2026-09-09.md for the
    larger question of whether the pressure-root slope measures the material at all.
    """
    path = os.path.join(HERE, "out", "phase8_verdict.txt")
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}; run analyse_phase8.py first")
    in_e2 = False
    for line in open(path, encoding="utf-8"):
        t = line.strip()
        if t.startswith("E2:"):
            in_e2 = True
            continue
        if in_e2 and (t.startswith("E3:") or t.startswith("E4:")):
            break
        if in_e2 and t.startswith("ungated") and "slope" in t:
            return float(t.split("slope")[1].split("/CV")[0]), "ungated E2"
    raise SystemExit(f"could not read an ungated E2 slope from {path}")


DPHI_C_DCV, DPHI_C_SOURCE = _dphi_c_dcv()
SALT_LADDER = ((30.0, "0 mM"), (10.0, "1 mM"), (3.0, "10 mM"), (1.0, "100 mM"))
STRESS_RATIOS = (1e-2, 1e-1)     # sigma/sigma*


def main() -> int:
    lines: list[str] = []

    def p(m: str = "") -> None:
        print(m)
        lines.append(m)

    p("=== THE DOUBLE LAYER NARROWS THE EFFECTIVE SIZE DISTRIBUTION ===")
    p("  d_eff = d + 2h with the same h for every particle, so sigma is unchanged and")
    p("  the mean rises: CV_eff = CV * d10 / (d10 + 2h).")
    p(f"  d10 = {D10 * 1e9:.0f} nm, h = lambda_D ln(sigma*/sigma)")
    p()
    p(f"  {'salt':>8s} {'lam_D[nm]':>10s} {'sigma/sigma*':>13s} {'h[nm]':>8s} "
      f"{'CV_eff/CV':>10s} {'CV effect':>11s}")
    for lam_nm, salt in SALT_LADDER:
        for ratio in STRESS_RATIOS:
            h = lam_nm * 1e-9 * math.log(1.0 / ratio)
            fac = D10 / (D10 + 2 * h)
            p(f"  {salt:>8s} {lam_nm:10.0f} {ratio:13.0e} {h * 1e9:8.2f} "
              f"{fac:10.4f} {(1 - fac) * 100:9.1f} % weaker")
    p()
    p("  🔑 THE PREDICTION: the SAME powder, at the SAME true CV, should show a")
    p("  measurably SMALLER polydispersity effect at low salt than at high salt.")
    p()
    lo = min(D10 / (D10 + 2 * lam * 1e-9 * math.log(1.0 / r))
             for lam, _ in SALT_LADDER for r in STRESS_RATIOS)
    hi = max(D10 / (D10 + 2 * lam * 1e-9 * math.log(1.0 / r))
             for lam, _ in SALT_LADDER for r in STRESS_RATIOS)
    p(f"  Across the lab's own 0-100 mM ladder the factor runs {lo:.3f} to {hi:.3f},")
    p(f"  i.e. the CV effect is up to {(1 - lo) * 100:.0f} % weaker at the low-salt end.")
    p()
    p("  In the quantity the mechanism actually uses, with the MEASURED slope")
    p(f"  dphi_c/dCV = {DPHI_C_DCV:.4f} ({DPHI_C_SOURCE}, read from "
      f"out/phase8_verdict.txt, not hardcoded):")
    p(f"  {'true CV':>9s} {'dphi_c at 100 mM':>18s} {'dphi_c at 0 mM':>16s} "
      f"{'difference':>12s}")
    for cv in (0.2, 0.35, 0.5):
        d_hi = DPHI_C_DCV * cv * hi
        d_lo = DPHI_C_DCV * cv * lo
        p(f"  {cv:9.2f} {d_hi:18.4f} {d_lo:16.4f} {d_hi - d_lo:12.4f}")
    p()
    p("  ⚠ IS IT RESOLVABLE? The seed scatter on dphi_c is 0.003 (out/phase8_verdict.txt),")
    p("  so a difference of 0.008-0.02 is 3-7 seed scatters IN THE SIMULATION. Whether")
    p("  the BENCH can resolve it is a different question and is not answered here.")
    p()
    p("=== ⚠ AND IT POINTS THE OPPOSITE WAY TO THE OTHER SALT EFFECT ===")
    p("  The same effective-diameter model has screening move phi_c^eff ITSELF from")
    p("  0.23 (0 mM) to 0.56 (100 mM) -- a LARGE effect that pushes the material toward")
    p("  the contractive side as salt is added. The narrowing effect above is a SMALL")
    p("  effect that makes the CV lever weaker as salt is REMOVED. On the salt axis the")
    p("  two are confounded, and separating them needs the CV ladder run AT TWO SALTS")
    p("  rather than a single salt series.")
    p()
    p("  🔑 That is a concrete experimental design: a 2 x 2 of (CV low/high) x (salt")
    p("  low/high). The interaction term is the narrowing prediction; the salt main")
    p("  effect is the phi_c^eff shift. Neither is identifiable from a one-dimensional")
    p("  ladder, which is what the lab currently has.")
    p()
    p("  ⚠⚠ PROVENANCE: all of this is a consequence of the effective-diameter ESTIMATE,")
    p("  which phi_axis_offset.py labels a reachability estimate and NOT a value of")
    p("  phi_c. It is a prediction OF THAT MODEL. The one measured quantity it uses,")
    p(f"  dphi_c/dCV = {DPHI_C_DCV:.4f} ({DPHI_C_SOURCE}), is real -- but see")
    p("  PREREG_ROOT_PINNING_2026-09-09.md: the pressure-root estimator that produced")
    p("  it follows the bottom of its own fitting window, so the slope may be a")
    p("  restatement of the isotropic dphi_J/dCV rather than an independent number.")
    p("  Nothing on this page depends on WHICH of the two it is -- both are the same")
    p("  value -- but a reader should not cite it as two agreeing measurements.")

    out = os.path.join(HERE, "out", "salt_narrows_cv.txt")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwritten -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
