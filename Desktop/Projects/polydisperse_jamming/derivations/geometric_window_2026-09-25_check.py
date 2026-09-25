"""geometric_window_2026-09-25_check.py -- checks for derivations/geometric_window_2026-09-25.md.

Writes derivations/geometric_window_2026-09-25_check.txt. Sections:
  C1  isostatic counting Z_iso = 2 n_f / n_c for every contact class (SymPy)
  C2  the window from the Voronoi map phi = Z/(Z + 2 sqrt3 w): closed form, limits, dchi/dw sign (SymPy)
  C3  clump lever arm <e>/r_eq = (pi/4)(alpha-1)(d_s/d_eq): SymPy average + Monte Carlo of the geometry
  C4  polyhedral contact eccentricity (faceted-grain proxy): Monte Carlo over random orientations
  C5  drainage: Kozeny-Carman with sphericity, size vs shape leverage (numbers)
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from shape import sub_sphere_diameter  # noqa: E402

L = []


def say(s=""):
    print(s)
    L.append(s)


rng = np.random.default_rng(20260925)

# ---------------------------------------------------------------- C1
say("C1  isostatic counting, Z_iso = 2 n_f / n_c")
nf, nc = sp.symbols("n_f n_c", positive=True)
Ziso = 2 * nf / nc
cases = [("frictionless sphere (rotations are null modes)", 3, 1),
         ("frictionless general body", 6, 1),
         ("frictional sphere, sliding blocked (mu -> inf)", 6, 3),
         ("+ rolling blocked", 6, 5),
         ("+ rolling and twisting blocked", 6, 6)]
for name, f, c in cases:
    say(f"   {name:48s} n_f {f} n_c {c}: Z_iso = {Ziso.subs({nf: f, nc: c})}")
xs, xr, xt = sp.symbols("x_s x_r x_t", nonnegative=True)
nc_mob = 1 + 2 * (1 - xs) + 2 * (1 - xr) + (1 - xt)
say(f"   partial mobilisation: n_c = {sp.expand(nc_mob)} (x_* = fraction of contacts at that plastic limit)")
assert sp.simplify(nc_mob.subs({xs: 0, xr: 1, xt: 1}) - 3) == 0
assert sp.simplify(nc_mob.subs({xs: 0, xr: 0, xt: 0}) - 6) == 0
say("   check: x_r = x_t = 1 (free rolling/twisting), x_s = 0 -> n_c = 3 (Z 4); all blocked -> n_c = 6 (Z 2): OK")

# ---------------------------------------------------------------- C2
say("\nC2  window from the Voronoi map phi = Z/(Z + 2 sqrt3 w)  [SWM 2008 doi:10.1038/nature06981 at w = 1, Z in [4,6]]")
Z, Zm, Z0, w = sp.symbols("Z Z_m Z_0 w", positive=True)
a = 2 * sp.sqrt(3)
phi = Z / (Z + a * w)
chi = phi.subs(Z, Zm) / phi.subs(Z, Z0)
chi_s = sp.simplify(chi)
say(f"   chi = {chi_s}")
target = Zm * (Z0 + a * w) / (Z0 * (Zm + a * w))
assert sp.simplify(chi_s - target) == 0
dchi_dw = sp.simplify(sp.diff(chi_s, w))
say(f"   dchi/dw = {dchi_dw}")
num = sp.factor(sp.numer(sp.together(dchi_dw)))
say(f"   numerator factored: {num}  -> negative when Z_m < Z_0: better filling (smaller w) WIDENS? sign check below")
for zm, z0, ww in ((4, 6, 1.0), (4, 6, 0.9), (2.4, 6, 1.0)):
    v = float(dchi_dw.subs({Zm: zm, Z0: z0, w: ww}))
    say(f"      dchi/dw at Z_m {zm}, Z_0 {z0}, w {ww}: {v:+.4f}")
say(f"   limits: w -> 0: chi -> {sp.limit(chi_s, w, 0)}; w -> oo: chi -> {sp.limit(chi_s, w, sp.oo)}; Z_m -> Z_0: chi -> {sp.simplify(chi_s.subs(Zm, Z0))}")
for zm in (4.0, 3.18, 2.5, 2.4, 2.0):
    say(f"   w = 1: Z_m {zm}: phi_m {float(phi.subs({Z: zm, w: 1})):.4f}, chi vs Z_0 6: {float(chi_s.subs({Zm: zm, Z0: 6, w: 1})):.4f}")
say("   ⚠ SWM state the map for Z in [4, 6]; Z_m < 4 (rolling blocked) is an EXTRAPOLATION of it.")

# ---------------------------------------------------------------- C3
say("\nC3  clump lever arm")
psi = sp.symbols("psi", positive=True)
mean_sin = sp.integrate(sp.sin(psi) * sp.sin(psi), (psi, 0, sp.pi)) / sp.integrate(sp.sin(psi), (psi, 0, sp.pi))
say(f"   <|sin psi|> over isotropic normals = {sp.simplify(mean_sin)} = {float(mean_sin):.6f}")
assert sp.simplify(mean_sin - sp.pi / 4) == 0


def lever_mc(alpha, n=400_000):
    """Distance from the centroid of a two-sphere clump (members of diameter ds at +-c on
    the z axis) to the normal line through a member centre, over isotropic normals, in r_eq."""
    ds = sub_sphere_diameter(1.0, alpha)
    c = (alpha - 1.0) * ds / 2.0
    nvec = rng.normal(size=(n, 3))
    nvec /= np.linalg.norm(nvec, axis=1)[:, None]
    sgn = np.where(rng.random(n) < 0.5, 1.0, -1.0)
    centre = np.zeros((n, 3))
    centre[:, 2] = sgn * c
    # only normals on the member's OUTER hemisphere are exposed; exposure does not change |sin psi|
    # symmetry, so weight all (the formula's assumption, printed)
    e = np.linalg.norm(np.cross(centre, nvec), axis=1)
    return float(e.mean() / 0.5), ds


say("   alpha   d_s/d_eq   formula (pi/4)(a-1)(ds/deq)   Monte Carlo   rel err")
for alpha in (1.1, 1.2, 1.3, 1.5, 1.75, 2.0):
    mc, ds = lever_mc(alpha)
    f = (math.pi / 4) * (alpha - 1) * ds
    say(f"   {alpha:<6} {ds:.4f}     {f:.4f}                       {mc:.4f}        {abs(mc - f) / f:.2e}")
    assert abs(mc - f) / f < 5e-3
say("   ⚠ assumes contact normals isotropic over each member, including the occluded neck; a heuristic, not a derivation of eta")

# ---------------------------------------------------------------- C4
say("\nC4  polyhedral contact eccentricity: vertex-on-face contacts, <|v_support(n) x n|>/r_eq, isotropic n")
phi_g = (1 + 5 ** 0.5) / 2


def solid(name):
    if name == "cube":
        V = np.array([[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)], float)
    elif name == "octahedron":
        V = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float)
    elif name == "dodecahedron":
        V = [[x, y, z] for x in (-1, 1) for y in (-1, 1) for z in (-1, 1)]
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                V += [[0, s1 / phi_g, s2 * phi_g], [s1 / phi_g, s2 * phi_g, 0], [s1 * phi_g, 0, s2 / phi_g]]
        V = np.array(V, float)
    elif name == "icosahedron":
        V = []
        for s1 in (-1, 1):
            for s2 in (-1, 1):
                V += [[0, s1, s2 * phi_g], [s1, s2 * phi_g, 0], [s1 * phi_g, 0, s2]]
        V = np.array(V, float)
    return V


def hull_volume(V):
    from scipy.spatial import ConvexHull
    return ConvexHull(V).volume


for name in ("cube", "octahedron", "dodecahedron", "icosahedron"):
    V = solid(name)
    V = V - V.mean(axis=0)
    r_eq = (3 * hull_volume(V) / (4 * math.pi)) ** (1 / 3)
    n = rng.normal(size=(200_000, 3))
    n /= np.linalg.norm(n, axis=1)[:, None]
    idx = np.argmax(n @ V.T, axis=1)            # the support vertex in direction n
    e = np.linalg.norm(np.cross(V[idx], n), axis=1)
    circ = np.linalg.norm(V, axis=1).max()
    say(f"   {name:13s}: <e>/r_eq = {e.mean() / r_eq:.3f} (circumradius/r_eq {circ / r_eq:.3f})")
# rounded polyhedra (Minkowski sum with a sphere of radius rho): lever arm of the core, radius grows
V = solid("icosahedron")
V = V - V.mean(axis=0)
r_core = (3 * hull_volume(V) / (4 * math.pi)) ** (1 / 3)
n = rng.normal(size=(200_000, 3))
n /= np.linalg.norm(n, axis=1)[:, None]
e_core = np.linalg.norm(np.cross(V[np.argmax(n @ V.T, axis=1)], n), axis=1).mean()
for rho_frac in (0.0, 0.25, 0.5, 0.75, 0.9):
    rho = rho_frac / (1 - rho_frac) * r_core if rho_frac < 1 else 0
    say(f"   rounded icosahedron, rounding radius = {rho_frac:.2f} of total: <e>/r ~ {e_core / (r_core + rho):.3f}")
say("   (a Minkowski-rounded body's normal passes through the core support point; lever arm unchanged, size grows)")

# ---------------------------------------------------------------- C5
say("\nC5  drainage leverage: Kozeny-Carman kappa = psi^2 d_v^2 (1-phi)^3 / (180 phi^2)  [form UNVERIFIED here; ratios only]")
for dv_cs, psi_cs in ((10e-6, 1.0), (14e-6, 0.85), (10e-6, 0.8)):
    ratio = (psi_cs ** 2 * dv_cs ** 2) / (1.0 * (0.8e-6) ** 2)
    say(f"   cornstarch d_v {dv_cs*1e6:.0f} um psi {psi_cs}: kappa ratio to 0.8 um silica spheres (same phi) = {ratio:.0f}; shape part psi^2 = {psi_cs**2:.2f}")
say("   size carries ~150-300x, sphericity <= 1.56x (psi 0.8): shape is a minor drainage lever")

with open(os.path.join(HERE, "geometric_window_2026-09-25_check.txt"), "w", encoding="utf-8", newline="\n") as fh:
    fh.write("\n".join(L) + "\n")
