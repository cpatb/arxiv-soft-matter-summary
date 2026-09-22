"""Checks for derivations/jam_stress_2026-09-22.md. Run: python jam_stress_2026-09-22_check.py

C1  WC window: phi >= phiJ(f)  <=>  f >= f_c = (1-x)/(1-chi)                 (SymPy)
C2  Weibull activation: f_ss(sigma) = f_c  <=>  sigma = sigma* L^(-1/delta)  (SymPy)
C3  d ln sigma_jam / d delta = ln(L)/delta^2, so the sign flips at f_c = 1/e (SymPy)
C4  limits: delta -> inf, delta -> 0, f_c -> 1, f_c -> 0                    (SymPy)
C5  numeric: brute-force threshold on a sigma grid vs the closed form, 4 points
C6  matched-kappa form sigma_jam/kappa = 1/Q_{1-f_c}(g) vs brute force on a sample
C7  dimensions with pint, and the size exponent sigma* ~ d^(m-2)             (SymPy)
C8  composition: sigma_jam at the edge of the window diverges as (1-f_c)^(-1/delta)
"""

import math
import sys

import numpy as np
import pint
import sympy as sp

sys.path.insert(0, "..")
import actshape as A  # noqa: E402

ok = True


def report(name, cond, detail=""):
    global ok
    ok &= bool(cond)
    print(f"{name}: {'PASS' if cond else 'FAIL'}  {detail}")


x, chi, f, phi0 = sp.symbols("x chi f phi0", positive=True)
sig, sstar, delta, fc = sp.symbols("sigma sigma_star delta f_c", positive=True)

# C1 -------------------------------------------------------------------------------
phiJ = phi0 - (phi0 - chi * phi0) * f            # phiM = chi*phi0
fsol = sp.solve(sp.Eq(x * phi0, phiJ), f)
report("C1", len(fsol) == 1 and sp.simplify(fsol[0] - (1 - x) / (1 - chi)) == 0,
       f"f at phi = phiJ(f): {sp.simplify(fsol[0])}")

# C2 -------------------------------------------------------------------------------
L = sp.log(1 / fc)
sjam = sstar * L ** (-1 / delta)
# on the domain 0 < f_c < 1 write f_c = exp(-Lp) with Lp > 0 (SymPy cannot infer the
# sign of ln(1/f_c) from f_c > 0 alone -- a representation fix, not a physics change)
Lp = sp.symbols("L", positive=True)
resid = sp.simplify(sp.exp(-(sstar / sjam.subs(fc, sp.exp(-Lp))) ** delta) - sp.exp(-Lp))
report("C2", resid == 0, f"exp(-(s*/s_jam)^d) - f_c with f_c = e^-L = {resid}")
# and f_ss is increasing in sigma, so f_ss >= f_c exactly for sigma >= sigma_jam
dfs = sp.diff(sp.exp(-(sstar / sig) ** delta), sig)
report("C2b", sp.simplify(dfs * sig / (delta * (sstar / sig) ** delta
                                       * sp.exp(-(sstar / sig) ** delta))) == 1,
       "d f_ss/d sigma = (delta/sigma)(s*/sigma)^delta f_ss > 0")

# C3 -------------------------------------------------------------------------------
dl = sp.simplify(sp.diff(sp.log(sjam / sstar), delta))
report("C3", sp.simplify(dl - sp.log(L) / delta ** 2) == 0, f"d ln(s_jam/s*)/d delta = {dl}")
signs = [(v, float(dl.subs({fc: v, delta: 1.3}))) for v in (0.1, 0.25, math.exp(-1), 0.75)]
report("C3b", signs[0][1] > 0 and signs[1][1] > 0 and abs(signs[2][1]) < 1e-15
       and signs[3][1] < 0, f"sign at f_c = 0.1, 0.25, 1/e, 0.75: {signs}")

# C4 -------------------------------------------------------------------------------
r = sjam / sstar
lim_inf = [sp.limit(r.subs(fc, v), delta, sp.oo) for v in (sp.Rational(1, 4), sp.Rational(3, 4))]
lim0_deep = sp.limit(r.subs(fc, sp.Rational(1, 4)), delta, 0, "+")
lim0_edge = sp.limit(r.subs(fc, sp.Rational(3, 4)), delta, 0, "+")
epsq = sp.symbols("epsq", positive=True)
lim_f1 = sp.limit(r.subs(delta, sp.Rational(13, 10)).subs(fc, 1 - epsq), epsq, 0, "+")
lim_f0 = sp.limit(r.subs(delta, sp.Rational(13, 10)), fc, 0, "+")
report("C4", lim_inf == [1, 1] and lim0_deep == 0 and lim0_edge == sp.oo
       and lim_f1 == sp.oo and lim_f0 == 0,
       f"delta->inf {lim_inf}; delta->0 deep {lim0_deep}, edge {lim0_edge}; "
       f"f_c->1 {lim_f1}; f_c->0 {lim_f0}")

# C5 -------------------------------------------------------------------------------
grid = np.logspace(-3, 3, 2_000_001)
worst = 0.0
for fcv, dv in ((0.1, 0.7), (0.25, 1.3169), (0.5, 2.6), (0.9, 1.0)):
    fss = np.exp(-(1.0 / grid) ** dv)                     # sigma* = 1
    brute = grid[np.argmax(fss >= fcv)]
    closed = A.weibull_sigma_jam_ratio(fcv, dv)
    worst = max(worst, abs(brute / closed - 1))
report("C5", worst < 1e-5, f"worst |brute/closed - 1| over 4 points = {worst:.2e}")

# C6 -------------------------------------------------------------------------------
rng = np.random.default_rng(22)
g = rng.lognormal(0.0, 0.8, 200_000)                      # deliberately NOT Weibull
worst = 0.0
for fcv in (0.1, 0.37, 0.75):
    s = np.logspace(-2, 2, 400_001)
    frac = A.activation_curve(g, s)
    brute = s[np.argmax(frac >= fcv)]
    worst = max(worst, abs(brute / A.sigma_jam_over_kappa(g, fcv) - 1))
report("C6", worst < 1e-3, f"lognormal g, worst |brute/closed - 1| = {worst:.2e}")

# C7 -------------------------------------------------------------------------------
u = pint.UnitRegistry()
kappa = 1e-3 * u.newton / u.meter           # F_rep = kappa R*
fN = 1e-6 * u.newton
Rs = 1e-6 * u.meter
sigma = 10 * u.pascal
gq = (fN / (Rs * sigma)).to_base_units()
sj = (kappa / gq).to("pascal")
report("C7a", gq.dimensionality == u.meter.dimensionality
       and sj.dimensionality == u.pascal.dimensionality,
       f"g = f/(R* sigma) has {gq.dimensionality}; kappa/Q(g) has {sj.dimensionality}")
d, m, km, gh, rh = sp.symbols("d m kappa_m ghat rhat", positive=True)
# contact force f = sigma d^2 ghat; reduced radius R* = d rhat; threshold kappa_m R*^m
sig_thr = sp.solve(sp.Eq(sig * d**2 * gh, km * (d * rh) ** m), sig)[0]
expo = sp.simplify(sp.diff(sp.log(sig_thr), d) * d)
report("C7b", sp.simplify(expo - (m - 2)) == 0, f"d ln sigma_thr / d ln d = {expo}")

# C8 -------------------------------------------------------------------------------
eps = sp.symbols("epsilon", positive=True)
lead = sp.limit((r.subs(fc, 1 - eps)) * eps ** (1 / delta), eps, 0, "+")
report("C8", sp.simplify(lead - 1) == 0, f"lim (1-f_c)^(1/d) s_jam/s* = {lead}")

print("ALL PASS" if ok else "SOME CHECK FAILED")
sys.exit(0 if ok else 1)
