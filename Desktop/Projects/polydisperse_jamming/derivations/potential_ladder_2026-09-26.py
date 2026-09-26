"""potential_ladder_2026-09-26.py -- does SURFACE POTENTIAL (salt) move the lab silica between ball and hybrid?

A scaling estimate on the project's poroelastic model (impact_stress/pore_relaxation via
derivations/ball_hybrid_scaling.groups); no new physics. Writes derivations/potential_ladder_2026-09-26.txt.

Salt enters twice:
  1. the double-layer network modulus K_p = P_conf a / (6 lambda_D) (pore_relaxation.network_modulus_edl; its
     silica_poro bracket IS the lab salt axis: lambda_D 30 nm at nominal 0 mM -> 1 nm at 100 mM);
  2. the onset stress sigma* (P_conf is set to sigma* by silica_poro's own rule). Two conventions, both printed:
     [fixed] sigma* independent of salt; [bench] sigma* scaled by the lab salt series' onset ratio, x0.31 at 100 mM
     relative to 0 mM (a reading of recorded bench numbers, derivations/multiplet_knobs_2026-09-22.md §9),
     interpolated log-linearly in lambda_D between those two ends (an ASSUMPTION, stated).
Ball corner: De >> 1 and S << 1 (silica 7e5, 0.004); hybrid: De ~ 1-30, S ~ 1-20 (cornstarch 7, 4).
"""
from __future__ import annotations

import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ball_hybrid_scaling as bhs  # noqa: E402

pr = bhs.pr
L = []


def say(s=""):
    print(s)
    L.append(s)


base = pr.silica_poro(phi=0.44)
a = base.sys.a
ss0 = bhs.SIGMA_STAR_SILICA[1]          # the value silica_poro uses for P_conf (retracted fit; kept for comparability)
LAM0, LAM100, R100 = 30e-9, 1e-9, 0.31
say("potential_ladder -- silica 800 nm, phi 0.44, lab impact point (D0 2.5 mm, U0 2.3 m/s)")
say(f"   sigma*(0 mM) = {ss0:g} Pa (silica_poro's P_conf); bench onset ratio at 100 mM = {R100}")
say("   salt/mM  lambda_D/nm | [fixed sigma*] K_p/Pa   De        S       | [bench sigma*] sigma*/Pa  K_p/Pa   De        S")
for label, lam in (("~0 (resid.)", 30e-9), ("1", pr.debye_length(1.0)), ("10", pr.debye_length(10.0)),
                   ("100", pr.debye_length(100.0)), ("(1 nm end)", 1e-9)):
    row = []
    for conv in ("fixed", "bench"):
        if conv == "fixed":
            ss = ss0
        else:
            f = min(max(math.log(LAM0 / lam) / math.log(LAM0 / LAM100), 0.0), 1.0)
            ss = ss0 * R100 ** f
        K = pr.network_modulus_edl(ss, a, lam)
        p = pr.PoroSystem(sys=base.sys, K_p=K, K_p_low=K, K_p_high=K, K_p_origin="salt ladder")
        g = bhs.groups(p, **bhs.LAB["silica"])
        row.append((ss, K, g["De"], g["S"]))
    say(f"   {label:11s}  {lam*1e9:7.2f}    | {row[0][1]:9.3g}  {row[0][2]:8.3g}  {row[0][3]:7.3g} | "
        f"{row[1][0]:8.3g}   {row[1][1]:9.3g}  {row[1][2]:8.3g}  {row[1][3]:7.3g}")
say("\nreading: across the whole lab salt range S stays << 1 and De >> 1 under both conventions -- salt moves the")
say("ONSET (sigma*, ~x0.31) and stiffens the network ~x10-30, but does not reach the hybrid corner. That is consistent")
say("with the bench reading that the ARRESTED state barely moves with salt (x1.28 against a x6.35 total range).")
with open(os.path.join(HERE, "potential_ladder_2026-09-26.txt"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
