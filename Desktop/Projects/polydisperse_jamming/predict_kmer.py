"""predict_kmer.py -- blind predictions for PREREG_KMER_2026-09-22.md, from MEASURED endpoints.

Writes out/predict_kmer.txt. Written and run BEFORE any mixture cell; must not be
regenerated after the campaign lands (the file is committed with the prereg).

Endpoints are the shapechi records at CV = 0 (N = 1700 grains): alpha = 1 is the singlet
suspension and alpha = 2 (a TANGENT dumbbell) is the pure-doublet suspension of the same
primary -- packing fraction is scale-free, so the volume-preserving construction of the
shape axis and the same-primary construction of this axis agree on phi.

H_mix (ideal specific-volume mixing, no excess): 1/phi_J(w) = (1-w)/phi_J(0) + w/phi_J(1),
applied separately to phi_0 (mu_p = 0) and phi_m (mu_p = 3.0; and mu_p = 0.3 as the realistic
frictional branch). w = doublet SOLID-VOLUME fraction. Justification for "no excess": every
grain is built of the SAME primary, so there is no size-ratio filling (the body-size ratio is
2^(1/3) = 1.26, far below the ~2.4 at which filling turns on in this project's binaries).
Alternative registered beside it: a positive excess (mixing shapes packs denser than either),
the known behaviour of mildly aspherical grains near alpha ~ 1.4.

Isostatic counting (H_Z). A rigid body in 3D has 6 degrees of freedom; a frictionless body
feels only normal forces, and a tangent dumbbell's rotation about its own axis is unresisted,
so it has 5 effective dofs (a triangle has 6). Frictionless isostaticity gives
Z = 2 * <dof>; perfectly rough (mu -> inf) contacts carry 3 constraints each, so Z = 2*6/3 = 4
for EVERY shape. Measured endpoints: frictionless singlets 5.97 (iso 6), doublets 7.92 (iso 10,
HYPOSTATIC by 21 %, as aspherical grains are known to be); mu_p = 3: 4.04 and 3.88.
"""
import json
import glob
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")
RUNGS = (3e-5, 1e-3)
PHI_SCALE = 0.0021


def endpoints():
    e = {}
    for f in glob.glob(os.path.join(OUT, "shape", "phase2_shapechi*.jsonl")):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if r.get("failed") or r["cv"] != 0.0:
                continue
            e.setdefault((r["alpha"], r["mu_p"], r["p_hat_target"]), []).append(r)
    return {k: (float(np.mean([x["phi"] for x in v])), float(np.mean([x["z_nr"] for x in v])),
                len(v)) for k, v in e.items()}


def hmix(p1, p2, w):
    return 1.0 / ((1 - w) / p1 + w / p2)


def main() -> int:
    L = []

    def say(s=""):
        print(s)
        L.append(s)

    E = endpoints()
    say("predict_kmer -- blind predictions from the shapechi singlet (alpha 1) and doublet")
    say("(alpha 2, tangent) endpoints. w = doublet solid-volume fraction.")
    for rung in RUNGS:
        say("")
        say(f"=== rung p_hat = {rung:g} ===")
        say("  endpoint  phi(mu0)  phi(mu0.3)  phi(mu3)   Z(mu0)  Z(mu0.3)  Z(mu3)")
        for a, lab in ((1.0, "singlet"), (2.0, "doublet")):
            say(f"  {lab:8s}  " + "  ".join(f"{E[(a, m, rung)][0]:.4f}" for m in (0.0, 0.3, 3.0))
                + "   " + "  ".join(f"{E[(a, m, rung)][1]:.3f}" for m in (0.0, 0.3, 3.0)))
        say("  H_mix     w    phi0     phim(3)  chi(3)  window(3)  phim(.3)  chi(.3)  window(.3)")
        for w in (0.0, 0.25, 0.5, 0.75, 1.0):
            p0 = hmix(E[(1.0, 0.0, rung)][0], E[(2.0, 0.0, rung)][0], w)
            p3 = hmix(E[(1.0, 3.0, rung)][0], E[(2.0, 3.0, rung)][0], w)
            p03 = hmix(E[(1.0, 0.3, rung)][0], E[(2.0, 0.3, rung)][0], w)
            say(f"          {w:.2f}  {p0:.4f}   {p3:.4f}   {p3 / p0:.4f}  {1 - p3 / p0:.4f}"
                f"     {p03:.4f}    {p03 / p0:.4f}  {1 - p03 / p0:.4f}")
    # the critical doublet fraction that OPENS the window at a given x, at the quoted rung
    say("")
    say("=== w*(x): doublet fraction at which chi(w) = x (window opens), rung 3e-5, H_mix ===")
    r = 3e-5
    for mm, lab in ((3.0, "mu_m = 3"), (0.3, "mu_m = 0.3")):
        ws = np.linspace(0, 1, 100001)
        chi = hmix(E[(1.0, mm, r)][0], E[(2.0, mm, r)][0], ws) / \
            hmix(E[(1.0, 0.0, r)][0], E[(2.0, 0.0, r)][0], ws)
        # x is defined against the SINGLET phi_0 (what a bench weighs against)
        say(f"  {lab}: chi(0) = {chi[0]:.4f}, chi(1) = {chi[-1]:.4f}")
        for x in (0.80, 0.83, 0.85, 0.87, 0.89):
            ok = np.flatnonzero(chi <= x)
            say(f"     x = {x:.2f}: w* = {ws[ok[0]]:.3f}" if ok.size else
                f"     x = {x:.2f}: never (chi(1) > x)")
    say("")
    say("=== H_Z (isostatic counting), registered bars ===")
    say("  mu_p = 3: Z_nr = 4 +- 0.25 at EVERY composition (Z = 2*6/3 for any rigid body).")
    s0, d0 = E[(1.0, 0.0, r)][1], E[(2.0, 0.0, r)][1]
    say("  mu_p = 0: hypostatic; Z_nr interpolates LINEARLY IN NUMBER FRACTION n2 between")
    say(f"           {s0:.3f} (n2 = 0) and {d0:.3f} (n2 = 1), bar +- 0.30. Isostatic would be 6 + 4 n2.")
    say("")
    say("=== Onset (Love-Weber), mu_p = 0.3, p_hat 3e-5 ===")
    zs, zd = E[(1.0, 0.3, r)][1], E[(2.0, 0.3, r)][1]
    say(f"  external contacts PER PRIMARY: singlet {zs:.3f}, doublet {zd / 2:.3f} "
        f"(ratio {zd / 2 / zs:.3f}).")
    say("  At fixed p and phi the mean contact force scales as 1/(contacts per volume x mean")
    say("  lever arm). Fewer contacts per primary raises <f> by 1/ratio; the lever arm from a")
    say("  dumbbell's centroid to a lobe contact is ~sqrt(2) x a sphere's (a/2 lobe offset in a")
    say("  random direction), which lowers it. Net: <f> ratio in [1/(ratio*sqrt(2)), 1/ratio] ->")
    lo, hi = 1 / (zd / 2 / zs * math.sqrt(2)), 1 / (zd / 2 / zs)
    say(f"  <f>_doublet/<f>_singlet in [{lo:.3f}, {hi:.3f}], so sigma*(w=1)/sigma*(0) in "
        f"[{1 / hi:.3f}, {1 / lo:.3f}] (every contact is primary-primary, so F_rep is the same).")
    say("  REGISTERED: sigma*(w=1)/sigma*(0) < 1, and inside [0.55, 0.95].")
    with open(os.path.join(OUT, "predict_kmer.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
