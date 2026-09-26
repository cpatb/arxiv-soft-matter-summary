"""cornstarch_fraction_2026-09-25.py -- can FRACTIONATED cornstarch reach the silica 'ball' corner of (De, S)?

A scaling estimate, reusing the project's poroelastic machinery (impact_stress/pore_relaxation, via
derivations/ball_hybrid_scaling.groups) -- no new physics. Writes derivations/cornstarch_fraction_2026-09-25.txt.

Ladders (each printed, none chosen):
  1. size fraction at FIXED cornstarch chemistry: a (diameter) 5-20 um, Hertz K_p recomputed at each a with
     cornstarch_poro's own P_conf and E_p bracket (1e7-1e9 Pa);
  2. solvent viscosity at a = 10 um (water -> glycerol-water, mu_w 1e-3 .. 1 Pa s), K_p unchanged;
  3. the silica reference row from ball_hybrid_scaling (same operating point convention).
The cornstarch operating point is [M26]'s (D0 4 mm, U0 5.1 m/s); phi 0.42 as in cornstarch_poro.
Caveat carried from ball_hybrid_scaling: De uses the contact-patch drainage length 0.49 D0.
"""
from __future__ import annotations

import dataclasses
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ball_hybrid_scaling as bhs  # noqa: E402

pr, idm = bhs.pr, bhs.idm
L = []


def say(s=""):
    print(s)
    L.append(s)


def poro_at(a: float, mu_w: float | None = None, E: float | None = None) -> "pr.PoroSystem":
    base = pr.cornstarch_poro(phi=0.42)
    s = dataclasses.replace(base.sys, a=a) if mu_w is None else dataclasses.replace(base.sys, a=a, mu_w=mu_w)

    def K(Ep):
        return pr.network_modulus_affine(
            pr.contact_stiffness_hertz(Ep, 0.3, s.a, s.phi, pr.P_CONF_CORNSTARCH), s.a, s.phi)
    k = K(pr.E_STARCH_PA if E is None else E)
    return pr.PoroSystem(sys=s, K_p=k, K_p_low=K(1.0e7), K_p_high=K(1.0e9), K_p_origin="fraction ladder")


op = bhs.LAB["cornstarch"]
say("cornstarch_fraction -- (De, S) for fractionated cornstarch; silica 'ball' corner for reference")
sil = bhs.groups(pr.silica_poro(phi=0.44), **bhs.LAB["silica"])
say(f"  silica 800 nm phi 0.44 (lab point): De {sil['De']:.3g}, S {sil['S']:.3g}   <- the ball corner")
base = bhs.groups(pr.cornstarch_poro(phi=0.42), **op)
say(f"  cornstarch a 10 um phi 0.42 ([M26] point): De {base['De']:.3g}, S {base['S']:.3g}")

say("\n1. size fraction at fixed chemistry (Hertz K_p recomputed; E_p bracket 1e7 / default / 1e9 Pa)")
say("   a/um   kappa      K_p(default)   De [lo, def, hi E_p]            S [lo, def, hi]")
for a in (5e-6, 7e-6, 10e-6, 14e-6, 20e-6):
    p = poro_at(a)
    g = [bhs.groups(p, K_p=k, **op) for k in (p.K_p_low, p.K_p, p.K_p_high)]
    say(f"   {a*1e6:5.1f}  {p.sys.kappa:.2e}  {p.K_p:.3g}        "
        f"{g[0]['De']:.3g} / {g[1]['De']:.3g} / {g[2]['De']:.3g}        "
        f"{g[0]['S']:.3g} / {g[1]['S']:.3g} / {g[2]['S']:.3g}")

say("\n2. solvent viscosity at a = 10 um (K_p, S unchanged by mu_w in this model)")
say("   NOTE: bhs.groups() calls pore_relaxation_time WITHOUT mu_f (water default), so it is inert to mu_w;")
say("   tau_p is recomputed here with mu_f passed explicitly (checked: mu_w 1e-3 reproduces groups()).")
for mu in (1e-3, 1e-2, 1e-1, 1.0):
    p = poro_at(10e-6, mu_w=mu)
    s = p.sys
    Ld = pr.drainage_lengths(op["D0"], s.a, s.phi)["contact"]
    De = pr.drainage_number(pr.pore_relaxation_time(Ld, s.kappa, p.K_p, s.phi, mu_f=mu), op["D0"], op["U0"])
    g = bhs.groups(p, **op)
    if mu == 1e-3:
        assert abs(De / g["De"] - 1) < 1e-12, (De, g["De"])
    say(f"   mu_w {mu:7.3g} Pa s: De {De:.3g}, S {g['S']:.3g}")

say("\nreading: the ball corner is De >> 1 AND S << 1. Fractionation moves De as a^-2 at most and leaves S")
say("at O(1) (Hertz K_p is nearly size-free at fixed P_conf); a viscous solvent moves De but not S. Neither")
say("reaches S << 1, which in this model needs a SOFT contact network (EDL-like K_p), i.e. surface chemistry.")
with open(os.path.join(HERE, "cornstarch_fraction_2026-09-25.txt"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
