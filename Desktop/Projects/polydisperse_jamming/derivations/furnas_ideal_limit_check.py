"""Closed forms for the binary linear packing model in the ideal (r -> infinity) limit,
and the finite-ratio corrections -- SymPy, limiting cases, spot checks against the
shipped `mixture_model` functions.

Run:  python derivations/furnas_ideal_limit_check.py > derivations/furnas_ideal_limit_numbers.txt

The model (mixture_model.packing_density) for two classes, coarse (index 1) and fine
(index 2), fines volume fraction x, u = d_fine/d_coarse, monodisperse density beta:

    gamma_1 = beta / (1 - (1 - a) x)                          coarse-dominant branch
    gamma_2 = beta / (1 - [1 - beta + b beta (1 - 1/beta)] (1 - x))
            = beta / (1 - (1 - beta)(1 - b) (1 - x))          fine-dominant branch
    phi     = min(gamma_1, gamma_2)

with a = a(u), b = b(u) in [0, 1], a(0) = b(0) = 0, a(1) = b(1) = 1.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import mixture_model as mm  # noqa: E402  the SHIPPED functions, not a private copy

x, beta, a, b = sp.symbols("x beta a b", positive=True)

g1 = beta / (1 - (1 - a) * x)
g2 = beta / (1 - (1 - beta) * (1 - b) * (1 - x))

# --- 1. algebraic identity: the fine-branch bracket ----------------------------------
# The model's coefficient is 1 - beta + b*beta*(1 - 1/beta). Show it equals (1-beta)(1-b).
bracket_model = 1 - beta + b * beta * (1 - 1 / beta)
assert sp.simplify(bracket_model - (1 - beta) * (1 - b)) == 0
print("1. fine-branch coefficient  1 - beta + b beta (1 - 1/beta)  ==  (1 - beta)(1 - b)   [sympy: 0]")

# --- 2. branch crossing x* and phi_max, general a, b ---------------------------------
# gamma_1 increases with x, gamma_2 decreases with x (both denominators shrink/grow
# monotonically), so min(gamma_1, gamma_2) is maximised where they cross.
xs = sp.solve(sp.Eq(g1, g2), x)
assert len(xs) == 1
xstar = sp.simplify(xs[0])
phimax = sp.simplify(g1.subs(x, xstar))
print("2. crossing        x* =", xstar)
print("   maximum   phi_max =", phimax)

# --- 3. ideal limit a = b = 0 ----------------------------------------------------------
xstar0 = sp.simplify(xstar.subs({a: 0, b: 0}))
phimax0 = sp.simplify(phimax.subs({a: 0, b: 0}))
print("3. ideal limit (a = b = 0):  x* =", xstar0, "   phi_max =", sp.factor(phimax0))
assert sp.simplify(xstar0 - (1 - beta) / (2 - beta)) == 0
assert sp.simplify(phimax0 - beta * (2 - beta)) == 0
print("   checks: x* == (1-beta)/(2-beta), phi_max == beta(2-beta)   [sympy: 0, 0]")

# --- 4. limiting cases ----------------------------------------------------------------
# (i) beta -> 1: nothing to fill. x* -> 0, phi_max -> 1.
print("4. limits: beta->1: x* =", sp.limit(xstar0, beta, 1), " phi_max =", sp.limit(phimax0, beta, 1))
# (ii) beta -> 0: a vanishingly dense pack; every void is fillable, x* -> 1/2.
print("          beta->0: x* =", sp.limit(xstar0, beta, 0), " phi_max/beta =", sp.limit(phimax0 / beta, beta, 0))
# (iii) a = b = 1 (u = 1, same size): both branches equal beta at every x.
assert sp.simplify(g1.subs(a, 1) - beta) == 0 and sp.simplify(g2.subs(b, 1) - beta) == 0
print("          a = b = 1: gamma_1 == gamma_2 == beta for every x   [sympy]")
# (iv) phi_max as a function of a alone at b = 0, and of b alone at a = 0
print("          phi_max(a, b=0) =", sp.simplify(phimax.subs(b, 0)))
print("          phi_max(a=0, b) =", sp.simplify(phimax.subs(a, 0)))
# (v) phi_max is decreasing in both a and b on (0,1): the sign of each partial.
dpa = sp.simplify(sp.diff(phimax, a))
dpb = sp.simplify(sp.diff(phimax, b))
print("          d phi_max/da =", dpa)
print("          d phi_max/db =", dpb)
# signs on the open unit square, numerically (the closed forms are not sign-definite by
# inspection): phi_max falls in both a and b; x* rises in a and falls in b.
fpa = sp.lambdify((a, b, beta), dpa)
fpb = sp.lambdify((a, b, beta), dpb)
fxa = sp.lambdify((a, b, beta), sp.diff(xstar, a))
fxb = sp.lambdify((a, b, beta), sp.diff(xstar, b))
A, B = np.meshgrid(np.linspace(0.001, 0.999, 60), np.linspace(0.001, 0.999, 60))
for bval in (0.5947, 0.6413, 0.55):
    print(f"          beta {bval}: max dphi_max/da {fpa(A, B, bval).max():+.2e}, "
          f"max dphi_max/db {fpb(A, B, bval).max():+.2e}, "
          f"min dx*/da {fxa(A, B, bval).min():+.2e}, max dx*/db {fxb(A, B, bval).max():+.2e}")
    assert fpa(A, B, bval).max() < 0 and fpb(A, B, bval).max() < 0
    assert fxa(A, B, bval).min() > 0 and fxb(A, B, bval).max() < 0

# --- 5. dimensional check ---------------------------------------------------------------
# Every symbol is a pure number (a volume fraction, a ratio of diameters, an exponent):
# no dimension can enter. Recorded, not tested.
print("5. dimensions: x, beta, a, b dimensionless -> x*, phi_max dimensionless. Nothing to check.")

# --- 6. numeric spot checks against the SHIPPED code -----------------------------------
print("6. spot checks against mixture_model.furnas_ideal / binary_density")
worst = 0.0
for bval in (0.5947, 0.6413, 0.55):
    xst = float(xstar0.subs(beta, bval))
    pm = float(phimax0.subs(beta, bval))
    ship = mm.furnas_ideal(xst, bval)
    # and the maximum of the shipped ideal curve on a fine grid must not exceed pm
    grid = np.linspace(0, 1, 100001)
    vals = np.array([mm.furnas_ideal(v, bval) for v in grid])
    igrid = int(np.argmax(vals))
    err = abs(ship - pm)
    worst = max(worst, err, abs(vals[igrid] - pm), abs(grid[igrid] - xst))
    print(f"   beta {bval}: x* {xst:.6f} phi_max {pm:.6f} | shipped at x*: {ship:.6f} "
          f"(|diff| {err:.1e}); shipped grid max {vals[igrid]:.6f} at x {grid[igrid]:.5f}")
print(f"   worst discrepancy {worst:.2e}")
assert worst < 2e-5

# finite a, b: the general crossing against a root-find on the shipped binary_density
print("   finite-ratio crossing against the shipped binary_density (PA, PB = 0.961, 2.214):")
pa_, pb_ = 0.9613854761022634, 2.214433963895545
bval = 0.5946865273380492
for r in (2.0, 3.0, 4.0, 5.0):
    u = 1.0 / r
    av = float(np.sqrt(1 - (1 - u) ** pa_))      # a(u) as in mixture_model.a_loosening
    bv = float(1 - (1 - u) ** pb_)
    xst = float(xstar.subs({a: av, b: bv, beta: bval}))
    pm = float(phimax.subs({a: av, b: bv, beta: bval}))
    # shipped: the crossing is where the two branch values are equal; bracket by the
    # argmax of the shipped curve
    grid = np.linspace(1e-4, 1 - 1e-4, 20001)
    vals = np.array([mm.binary_density(r, v, bval, pa_, pb_) for v in grid])
    i = int(np.argmax(vals))
    print(f"   r {r}: a {av:.4f} b {bv:.4f} -> x* {xst:.4f} phi_max {pm:.4f} | shipped argmax "
          f"x {grid[i]:.4f} phi {vals[i]:.4f}  (|dphi| {abs(vals[i]-pm):.1e}, |dx| {abs(grid[i]-xst):.1e})")
    assert abs(vals[i] - pm) < 1e-5 and abs(grid[i] - xst) < 1e-3
print("ALL CHECKS PASSED")
