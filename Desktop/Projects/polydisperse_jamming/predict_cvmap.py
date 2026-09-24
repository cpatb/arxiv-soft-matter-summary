"""predict_cvmap.py -- inputs and blind predictions for PREREG_CVMAP_2026-09-24.md.

The project's FOUNDING observation (author, bench, 2026-08-31): at similar particle size,
raising the size-distribution width CV moves an impacting drop from purely dilatant to
spreading and inertial-dominated. The DEM + closure say CV acts through x = phi/phi_0 (it
raises phi_0 ~9.9 % while leaving chi = phi_m/phi_0 nearly flat), so at a FIXED WEIGHED phi
raising CV walks the suspension DOWN and OUT of the shear-jamming window. This has only ever
been COMPOSED (jam_map.py section 4); it has never been run through the drop solver.

This script reads every input from the graded artifacts (no typed chi/phi values):
  window endpoints  -- three_knobs.load() at the quoted rung, via jam_map.phi_endpoints
  delta, sigma* (CV 0.2)  -- out/actshape_cv02.txt, mu 0.3, p 3e-4, derjaguin
  delta, sigma* (CV 0.5)  -- out/grade_actshape_fs.txt F2/F3 (N = 4000, N-independent)
  sigma* (CV 0.5, d32)    -- out/grade_actshape.txt H2 (the size-convention bracket)
and writes out/predict_cvmap.txt plus out/cvmap_cells.txt (the launcher's cell list).
Run BEFORE any cell; must not be regenerated after the campaign lands.
"""
from __future__ import annotations

import math
import os
import re

import jam_map as jm
import three_knobs as tk

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
PHI_SCALE = 0.0021          # protocol run-to-run scatter (mixture_model.py), the grading scale
DNS_PHI0 = 0.62             # the DNS's own phi_0 (kmermap convention); x and chi transfer, phi does not
DELTA_REF, SSTAR_REF = 1.3169, 300.0   # cs_fit singlet reference (kmermap convention)
MP_EXPO = 2
X0S = (0.92, 0.932, 0.95)   # x at CV 0; the weighed phi is held fixed across CV
U0S = (2.0, 4.0)
CVS = (0.0, 0.2, 0.5)


def _grab(path, pattern):
    txt = open(os.path.join(OUT, path), encoding="utf-8").read()
    m = re.search(pattern, txt, re.M)
    if not m:
        raise SystemExit(f"pattern not found in {path}: {pattern}")
    return float(m.group(1))


def activation():
    """delta and sigma* RATIOS (CV / CV 0), mu_p 0.3, derjaguin, from the graded text."""
    d02 = _grab("actshape_cv02.txt", r"^\s*0\.3\s+3e-04\s+derjaguin\s+d1090(?:\s+\S+){5}\s+(\S+)")
    s02 = _grab("actshape_cv02.txt", r"^\s*0\.3\s+3e-04\s+derjaguin\s+sstar(?:\s+\S+){5}\s+(\S+)")
    d05 = _grab("grade_actshape_fs.txt", r"F2 CV0 \(N1000\) -> CV0\.5 \(N4000\): ratio (\S+),")
    s05 = _grab("grade_actshape_fs.txt", r"F3 sigma\*/kappa ratio CV0\.5\(N4000\)/CV0\(N1000\) at d10 = (\S+)\s")
    s05_32 = _grab("grade_actshape.txt", r"derjaguin mu 0\.3: .* at d32 (\S+)\s")
    return {0.0: (1.0, 1.0), 0.2: (d02, s02), 0.5: (d05, s05)}, s05_32


def theta(x, chi, n=MP_EXPO):
    """Thickening ratio [(1-x)/(1-x/chi)]^n, finite only for x < chi (outside the window)."""
    return ((1 - x) / (1 - x / chi)) ** n if x < chi else math.inf


def main() -> int:
    lines = []

    def say(s=""):
        print(s)
        lines.append(s)

    best = tk.load()
    ends = {}
    for cv in CVS:
        p0, s0, n0 = jm.phi_endpoints(best, cv, 0.0)
        pm, sm, nm = jm.phi_endpoints(best, cv, 0.3)
        ends[cv] = dict(phi0=p0, s0=s0, n0=n0, phim=pm, sm=sm, nm=nm, chi=pm / p0)
    act, s05_32 = activation()

    say("predict_cvmap -- inputs and blind predictions for PREREG_CVMAP_2026-09-24.md")
    say(f"window endpoints at p_hat = {jm.RUNG:g}, converged records, mu_m = 0.3 (the clean branch)")
    say("   CV   phi_0 (sd, n)            phi_m(mu .3) (sd, n)     chi      delta ratio  sigma* ratio (d10)")
    for cv in CVS:
        e = ends[cv]
        say(f"  {cv:.1f}  {e['phi0']:.5f} ({e['s0']:.5f}, {e['n0']})   "
            f"{e['phim']:.5f} ({e['sm']:.5f}, {e['nm']})   {e['chi']:.5f}   "
            f"{act[cv][0]:.3f}        {act[cv][1]:.3f}")
    say(f"   sigma* ratio CV 0.5 at fixed d32 instead: {s05_32:.3f}  (the convention bracket; H2)")
    say("   ⚠ phi_0(CV 0.5) has n = 1 converged seed (the frictionless arm caps); a capped cell reads")
    say("     LOW, so the CV shift below is if anything UNDER-stated.")

    # ---- 1. The founding observation as a closed-form inequality -------------------------
    r05 = ends[0.0]["phi0"] / ends[0.5]["phi0"]
    chi05 = ends[0.5]["chi"]
    say("")
    say("=== 1. CLOSED FORM: can a CV 0.5 suspension jam at the weighed phi of a CV 0 one? ===")
    say("   At fixed weighed phi, x(CV) = x0 * phi_0(0)/phi_0(CV). A CV-0 suspension that flows at rest")
    say("   has x0 < 1, so the CV-0.5 one can jam somewhere only if phi_0(0)/phi_0(0.5) > chi(0.5).")
    say(f"   phi_0(0)/phi_0(0.5) = {r05:.5f}   chi(0.5) = {chi05:.5f}   margin chi - ratio = {chi05 - r05:+.5f}")
    # uncertainty of the margin: PHI_SCALE on each of the three phi values it is built from
    dm = math.sqrt((PHI_SCALE / ends[0.5]["phi0"] * r05) ** 2          # phi_0(0.5) in the ratio
                   + (PHI_SCALE / ends[0.0]["phi0"] * r05) ** 2        # phi_0(0) in the ratio
                   + (PHI_SCALE / ends[0.5]["phi0"]) ** 2              # phi_m(0.5) in chi
                   + (PHI_SCALE / ends[0.5]["phi0"] * chi05) ** 2)     # phi_0(0.5) in chi
    say(f"   margin / (PHI_SCALE propagated, {dm:.5f}) = {(chi05 - r05) / dm:.2f}   "
        "(phi_0(0.5) enters both terms with the same sign, so this sum-in-quadrature is conservative)")
    if chi05 > r05:
        say("   => at mu_m 0.3 NO weighed phi exists at which a CV-0-fluid suspension, re-made at CV 0.5, can")
        say("      shear-jam: the whole CV-0 window maps below chi(0.5). The founding observation, in closed form.")
    else:
        x_band = chi05 / r05
        say(f"   => a band x0 in [{x_band:.4f}, 1) remains in which CV 0.5 can still jam.")
    for mu_lab, mu in (("mu_m 3", 3.0),):
        p0 = ends[0.5]["phi0"]
        pm3 = jm.phi_endpoints(best, 0.5, mu)[0]
        c3 = pm3 / p0
        say(f"   contrast at {mu_lab}: chi(0.5) = {c3:.5f} < ratio {r05:.5f}, so a band x0 in "
            f"[{c3 / r05:.4f}, 1) CAN still jam -- the claim is friction-scoped.")

    # ---- 2. The grid --------------------------------------------------------------------
    say("")
    say("=== 2. THE GRID (fixed weighed phi per x0; DNS phi_0 = 0.62, so PHI = x*0.62, PHIM = chi*0.62) ===")
    say("   x0     CV   x        chi      in window  f_c      Theta(n=2)  delta     sigma*")
    cells = []
    for x0 in X0S:
        for cv in CVS:
            x = x0 * ends[0.0]["phi0"] / ends[cv]["phi0"]
            chi = ends[cv]["chi"]
            inside = chi <= x < 1
            fc = (1 - x) / (1 - chi) if inside else 1.0
            d = DELTA_REF * act[cv][0]
            ss = SSTAR_REF * act[cv][1]
            say(f"  {x0:.3f}  {cv:.1f}  {x:.5f}  {chi:.5f}  {'yes' if inside else 'no ':3s}        "
                f"{fc:.4f}   {theta(x, chi):9.2f}   {d:.4f}   {ss:.1f}")
            for u0 in U0S:
                cells.append((f"c{cv:.1f}_x{x0:.3f}_U{u0:.1f}", x, d, ss, fc, chi, u0))
    # decomposition arms: at x0 0.932 and 0.95, U0 4, CV 0.5
    say("")
    say("=== 3. CHANNEL SPLIT (U0 4, CV 0 -> 0.5, every x0) ===")
    say("   W (window only): x, chi from CV 0.5; delta, sigma* from CV 0.")
    say("   A (activation only): x, chi from CV 0; delta, sigma* (d10) from CV 0.5.")
    say("   D (full, d32 convention): as the full CV 0.5 cell but sigma* at fixed d32.")
    for x0 in X0S:
        x5 = x0 * ends[0.0]["phi0"] / ends[0.5]["phi0"]
        chi0, chi5 = ends[0.0]["chi"], ends[0.5]["chi"]
        fc0 = (1 - x0) / (1 - chi0) if chi0 <= x0 < 1 else 1.0
        cells.append((f"W_x{x0:.3f}_U4.0", x5, DELTA_REF, SSTAR_REF, 1.0, chi5, 4.0))
        cells.append((f"A_x{x0:.3f}_U4.0", x0, DELTA_REF * act[0.5][0], SSTAR_REF * act[0.5][1],
                      fc0, chi0, 4.0))
        cells.append((f"D_x{x0:.3f}_U4.0", x5, DELTA_REF * act[0.5][0], SSTAR_REF * s05_32, 1.0, chi5, 4.0))
    with open(os.path.join(OUT, "cvmap_cells.txt"), "w", encoding="utf-8", newline="\n") as f:
        for lab, x, d, ss, fc, chi, u0 in cells:
            f.write(f"{lab} {x * DNS_PHI0:.5f} {d:.4f} {ss:.1f} {fc:.4f} {chi * DNS_PHI0:.5f} {u0:.1f}\n")
    say(f"\n{len(cells)} cells written to out/cvmap_cells.txt (label PHI DELTA SIGMA_STAR JAM_FHAT PHIM U0)")
    with open(os.path.join(OUT, "predict_cvmap.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
