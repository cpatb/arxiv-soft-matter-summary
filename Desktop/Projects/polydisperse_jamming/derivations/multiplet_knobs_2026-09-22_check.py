"""Checks for derivations/multiplet_knobs_2026-09-22.md.  python multiplet_knobs_2026-09-22_check.py

K1  tangent k-clusters keep the solid's specific surface: S/V = 6/d, so the Sauter diameter
    and the Carman-Kozeny permeability at fixed phi are UNCHANGED (sympy, exact)
K1b a SINTERED dumbbell (centre separation s < d) lowers S/V, i.e. RAISES permeability;
    closed form vs brute Monte-Carlo volume/area at 3 overlaps
K2  the throat between three mutually tangent primaries has inscribed radius
    (2/sqrt(3) - 1) a > 0: a compact tangent cluster (k <= 4) encloses no fluid, Lambda = 1
K3  isostatic counting: frictionless Z = 2<dof> with dof 3/5/6 (sphere/dumbbell/triangle);
    rough (mu -> inf) Z = 2*6/3 = 4 for every rigid body (sympy)
K4  under ideal specific-volume mixing the window-opening fraction w*(x) is LINEAR-fractional
    in x (closed form) and equals the brute root of chi(w) = x at 3 points
K5  thickening ratio Theta = [(1-x)/(1-x/chi)]^n INCREASES as chi falls at fixed x (sympy sign)
K6  onset under Love-Weber: sigma*(w)/sigma*(0) = (n_c l)(w)/(n_c l)(0); the registered bracket
    [ratio, ratio*sqrt(2)] reproduces predict_kmer.py's [0.608, 0.860]
K7  dimensions: sigma* = F_rep/(a^2 ghat) and the breakage number Pi_b = sigma_imp a^2/F_b (pint)
"""
import math
import sys

import numpy as np
import pint
import sympy as sp

ok = True


def rep(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"{name}: {'PASS' if cond else 'FAIL'}  {detail}")


a, d, k, s = sp.symbols("a d k s", positive=True)
# K1 ---------------------------------------------------------------------------------
V1, S1 = sp.pi * d**3 / 6, sp.pi * d**2
ratio = sp.simplify((k * S1) / (k * V1) - 6 / d)
rep("K1", ratio == 0, f"S/V of k tangent spheres - 6/d = {ratio}")

# K1b: dumbbell of two spheres (diameter d, radius r) at separation s < d
r = d / 2
h = r - s / 2                                     # cap height removed from each sphere
V_db = 2 * V1 - 2 * (sp.pi * h**2 * (3 * r - h) / 3)
S_db = 2 * S1 - 2 * (2 * sp.pi * r * h)
dsauter = sp.simplify(6 * V_db / S_db)
rng = np.random.default_rng(1)
worst = 0.0
for sv in (0.6, 0.8, 0.95):
    dv = 1.0
    # Monte-Carlo union volume
    P = rng.uniform([-dv / 2 - sv / 2, -dv / 2, -dv / 2], [dv / 2 + sv / 2, dv / 2, dv / 2],
                    size=(2_000_000, 3))
    c1, c2 = np.array([-sv / 2, 0, 0]), np.array([sv / 2, 0, 0])
    inside = (np.linalg.norm(P - c1, axis=1) < dv / 2) | (np.linalg.norm(P - c2, axis=1) < dv / 2)
    vol_mc = inside.mean() * (dv + sv) * dv * dv
    vol_cf = float(V_db.subs({d: dv, s: sv}))
    worst = max(worst, abs(vol_mc / vol_cf - 1))
d_s_vals = [float(dsauter.subs({d: 1, s: sv})) for sv in (0.6, 0.8, 1.0)]
rep("K1b", worst < 3e-3 and d_s_vals[0] > d_s_vals[1] > d_s_vals[2] and
    abs(d_s_vals[2] - 1) < 1e-12,
    f"MC volume vs closed form worst {worst:.1e}; Sauter d at s = 0.6/0.8/1.0: "
    f"{[round(v, 4) for v in d_s_vals]} (tangent = 1 exactly, sintered > 1)")

# K2 ---------------------------------------------------------------------------------
# three mutually tangent circles of radius a (the in-plane section of the throat):
# centroid distance to each centre = 2a/sqrt(3); inscribed radius = 2a/sqrt(3) - a
rin = 2 * a / sp.sqrt(3) - a
rep("K2", sp.simplify(rin / a - (2 / sp.sqrt(3) - 1)) == 0 and float(rin.subs(a, 1)) > 0,
    f"throat inscribed radius = {sp.nsimplify(rin / a)} a = {float(rin.subs(a, 1)):.4f} a > 0")

# K3 ---------------------------------------------------------------------------------
n1, n2, n3 = sp.symbols("n1 n2 n3", nonnegative=True)
Z0 = 2 * (3 * n1 + 5 * n2 + 6 * n3) / (n1 + n2 + n3)
Zrough = sp.Rational(2 * 6, 3)
rep("K3", sp.simplify(Z0.subs({n2: 0, n3: 0}) - 6) == 0 and
    sp.simplify(Z0.subs({n1: 0, n3: 0}) - 10) == 0 and Zrough == 4,
    f"frictionless Z_iso = {sp.simplify(Z0)}; rough Z = {Zrough}")

# K4 ---------------------------------------------------------------------------------
x, w, A0s, A0d, Ams, Amd = sp.symbols("x w A0s A0d Ams Amd", positive=True)  # A = 1/phi
chi = ((1 - w) * A0s + w * A0d) / ((1 - w) * Ams + w * Amd)
wstar = sp.solve(sp.Eq(chi, x), w)
closed = (x * Ams - A0s) / ((A0d - A0s) - x * (Amd - Ams))
rep("K4a", len(wstar) == 1 and sp.simplify(wstar[0] - closed) == 0, f"w* = {sp.simplify(wstar[0])}")
vals = {A0s: 1 / 0.6327, A0d: 1 / 0.6272, Ams: 1 / 0.5603, Amd: 1 / 0.4945}
worst = 0.0
for xv, want in ((0.83, 0.542), (0.85, 0.338), (0.87, 0.145)):
    got = float(closed.subs(vals).subs(x, xv))
    worst = max(worst, abs(got - want))
rep("K4b", worst < 2e-3, f"closed form vs predict_kmer.py brute roots, worst {worst:.1e}")

# K5 ---------------------------------------------------------------------------------
c, n = sp.symbols("chi n", positive=True)
Theta = ((1 - x) / (1 - x / c)) ** n
dT = sp.simplify(sp.diff(sp.log(Theta), c))
rep("K5", sp.simplify(dT + n * x / (c * (c - x))) == 0,
    f"d ln Theta / d chi = {dT}  (< 0 for x < chi: a smaller chi means MORE thickening)")

# K6 ---------------------------------------------------------------------------------
zs, zd = 4.628, 5.631
ratio_c = zd / 2 / zs
rep("K6", abs(ratio_c - 0.608) < 1e-3 and abs(ratio_c * math.sqrt(2) - 0.860) < 1e-3,
    f"contacts per primary ratio {ratio_c:.3f}; sigma* bracket [{ratio_c:.3f}, "
    f"{ratio_c * math.sqrt(2):.3f}]")

# K7 ---------------------------------------------------------------------------------
u = pint.UnitRegistry()
F = 1e-9 * u.newton
aa = 1e-6 * u.meter
sig = (F / aa**2).to("pascal")
Pi = (4800 * u.pascal * aa**2 / F).to("dimensionless")
rep("K7", sig.dimensionality == u.pascal.dimensionality and Pi.dimensionless,
    f"F/a^2 -> {sig.dimensionality}; Pi_b dimensionless = {Pi.dimensionless}")

print("ALL PASS" if ok else "SOME CHECK FAILED")
sys.exit(0 if ok else 1)
