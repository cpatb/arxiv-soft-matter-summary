#!/bin/bash
# probe_clump_s0.sh -- CAPABILITY PROBE ONLY. Runs no production cell, writes no shard.
#
# scripts/probe_clump.sh asked four questions and found two defects (utils.porosity()
# and clump mass both double-count the member overlap). It asked whether clump members
# contact EACH OTHER. It did not ask the mirror question:
#
#   does an EXTERNAL neighbour contact BOTH members of the same clump, and therefore
#   appear TWICE in the contact graph?
#
# That question decides three things the shape campaign rests on:
#
#   Q1  whether alpha = 1 built as a clump of two COINCIDENT spheres is the sphere the
#       project has already run, or a body with every contact doubled -- i.e. whether
#       gate S0 (bit identity at alpha = 1) is achievable AS WRITTEN or must be run with
#       the degenerate single-sphere body;
#   Q2  whether z_all and n_contacts along the alpha ladder mix shape with a
#       contact-counting artefact, and how big that is per rung;
#   Q3  whether the timestep moves, because PWaveTimeStep reads the MEMBER radius and
#       the MEMBER density, and the density is exactly what the mass correction changes.
#
# Q4 is separate and cheap: can phi computed analytically reproduce 1 - utils.porosity()
# BITWISE for plain spheres? If it cannot, S0's bit-identity bar is unachievable for a
# reason that has nothing to do with clumps.
#
# Sourced BEFORE any `set -e` (lessons.md: sourcing under set -e dies silently).
source "$HOME/pj/scripts/env.sh"
set -u

SRC=/mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming
mkdir -p ~/pj/out ~/pj/scripts
for f in shape.py psd.py measure.py packing.py; do
  cp "$SRC/$f" ~/pj/"$f" || exit 1
  tr -d '\r' < ~/pj/"$f" > ~/pj/"$f".tmp && mv ~/pj/"$f".tmp ~/pj/"$f"
done
cd ~/pj || exit 1
echo "--- md5 in WSL (after CR strip) ---"
for f in shape.py psd.py measure.py packing.py; do echo "  $(md5sum "$f")"; done

python3 - <<'PYEOF'
import math
import sys

import numpy as np
import yade
from yade import O, utils

sys.path.insert(0, ".")
from shape import ALPHA_LADDER, alpha_to_u, clump_centres, span_factor, sub_sphere_diameter

MAT = dict(young=1.0, poisson=0.3, frictionAngle=0.3)


def fresh(L=8.0, density=1.0):
    O.reset()
    O.periodic = True
    O.cell.hSize = yade.Matrix3(L, 0, 0, 0, L, 0, 0, 0, L)
    return O.materials.append(utils.FrictMat(density=density, **MAT))


def engines():
    O.engines = [
        utils.ForceResetter(),
        utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb()],
                                    allowBiggerThanPeriod=False),
        utils.InteractionLoop(
            [utils.Ig2_Sphere_Sphere_ScGeom()],
            [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
            [utils.Law2_ScGeom_FrictPhys_CundallStrack()]),
        utils.NewtonIntegrator(damping=0.4),
    ]


def add_clump(centre, d_eq, alpha, axis, mat):
    """Two-sphere clump. At alpha = 1 the centres COINCIDE -- deliberately built the
    naive way here, because whether that is the sphere is exactly the question."""
    ds = sub_sphere_diameter(d_eq, alpha)
    off = clump_centres(d_eq, alpha, axis)
    ids = [O.bodies.append(utils.sphere(tuple(np.asarray(centre) + off[j]),
                                        ds / 2.0, material=mat)) for j in range(2)]
    cid = O.bodies.clump(ids)
    return cid, ids


def count_real():
    n = 0
    pairs = {}
    for i in O.interactions:
        if not i.isReal:
            continue
        n += 1
        pairs[(i.id1, i.id2)] = pairs.get((i.id1, i.id2), 0) + 1
    return n, pairs


print("=" * 78)
print("Q1 -- does an EXTERNAL neighbour contact BOTH members of one clump?")
print("=" * 78)
# A clump at the origin and ONE plain sphere overlapping it head-on. A sphere pair in
# contact gives exactly one interaction; anything more is the double count.
for alpha in ALPHA_LADDER:
    mat = fresh()
    d_eq = 1.0
    # axis along z; the neighbour is placed along +x (the WAIST direction), which is
    # where both members are equidistant and the double count is most likely.
    cid, ids = add_clump((4.0, 4.0, 4.0), d_eq, alpha, (0, 0, 1), mat)
    ds = sub_sphere_diameter(d_eq, alpha)
    # neighbour of the same sub-sphere size, just touching the waist
    gap = 0.98 * (ds / 2.0 + ds / 2.0)
    nb = O.bodies.append(utils.sphere((4.0 + gap, 4.0, 4.0), ds / 2.0, material=mat))
    engines()
    O.dt = 1e-6
    O.run(2, True)
    n, pairs = count_real()
    with_nb = {k: v for k, v in pairs.items() if nb in k}
    print(f"  alpha={alpha:.2f}  ds={ds:.6f}  real interactions={n}  "
          f"distinct member-neighbour pairs={len(with_nb)}  -> "
          + ("DOUBLE COUNT" if len(with_nb) > 1 else "single contact"))

print()
print("=" * 78)
print("Q2 -- in a REAL packing, how many external contacts are duplicates?")
print("=" * 78)
N = 200
for alpha in ALPHA_LADDER:
    mat = fresh(L=9.0)
    rng = np.random.default_rng(7)
    axes = rng.normal(size=(N, 3))
    axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    clump_of = {}
    for i in range(N):
        c = rng.uniform(0.5, 8.5, size=3)
        cid, ids = add_clump(c, 1.0, alpha, axes[i], mat)
        for b in ids:
            clump_of[b] = cid
    engines()
    O.dt = 0.3 * utils.PWaveTimeStep()
    O.run(200, True)
    n, pairs = count_real()
    # map each member-member interaction to the CLUMP pair behind it
    cpairs = {}
    internal = 0
    for (a, b), c in pairs.items():
        ca, cb = clump_of.get(a, -1), clump_of.get(b, -2)
        if ca == cb:
            internal += c
            continue
        key = (min(ca, cb), max(ca, cb))
        cpairs[key] = cpairs.get(key, 0) + c
    uniq = len(cpairs)
    dup = n - internal - uniq
    frac = dup / (n - internal) if (n - internal) else float("nan")
    print(f"  alpha={alpha:.2f}  real={n}  internal={internal}  unique clump pairs={uniq}"
          f"  duplicates={dup}  dup_frac={frac:.4f}  z_raw={2*(n-internal)/N:.3f}"
          f"  z_dedup={2*uniq/N:.3f}")

print()
print("=" * 78)
print("Q3 -- does PWaveTimeStep move with the member density?")
print("=" * 78)
# The mass correction scales the MATERIAL density by 1/(2*(ds/d_eq)^3) so the clump's
# mass is the union mass. PWaveTimeStep reads the member radius and member density.
mat = fresh(density=1.0)
for i in range(50):
    O.bodies.append(utils.sphere((0.1 + 0.15 * i, 4.0, 4.0), 0.5, material=mat))
engines()
dt_sphere = utils.PWaveTimeStep()
print(f"  plain spheres, r=0.5, density=1.0        : PWaveTimeStep = {dt_sphere:.12e}")
for alpha in ALPHA_LADDER:
    ds = sub_sphere_diameter(1.0, alpha)
    rho = 1.0 / (2.0 * (ds / 1.0) ** 3)
    mat = fresh(density=rho)
    rng = np.random.default_rng(3)
    axes = rng.normal(size=(20, 3))
    axes /= np.linalg.norm(axes, axis=1, keepdims=True)
    for i in range(20):
        add_clump(rng.uniform(1, 7, size=3), 1.0, alpha, axes[i], mat)
    engines()
    dt = utils.PWaveTimeStep()
    print(f"  alpha={alpha:.2f} ds={ds:.6f} density={rho:.6f} : "
          f"PWaveTimeStep = {dt:.12e}   ratio to a d=1 sphere at density 1 = "
          f"{dt / (0.5 * math.sqrt(1.0 / 1.0)):.6f}")

print()
print("=" * 78)
print("Q4 -- can an ANALYTIC phi reproduce 1 - utils.porosity() BITWISE for spheres?")
print("=" * 78)
rng = np.random.default_rng(11)
radii = rng.uniform(0.2, 0.6, size=400)
L = 9.0
mat = fresh(L=L)
for r in radii:
    O.bodies.append(utils.sphere(tuple(rng.uniform(0, L, 3)), float(r), material=mat))
phi_yade = 1.0 - utils.porosity()
v_np = float(np.sum(4.0 / 3.0 * math.pi * radii ** 3))
phi_np = v_np / L ** 3
v_loop = 0.0
for r in radii:
    v_loop += 4.0 / 3.0 * math.pi * float(r) ** 3
phi_loop = v_loop / L ** 3
print(f"  1 - utils.porosity()      = {phi_yade!r}")
print(f"  numpy sum / L^3           = {phi_np!r}   bitwise equal: {phi_np == phi_yade}")
print(f"  python loop sum / L^3     = {phi_loop!r}   bitwise equal: {phi_loop == phi_yade}")
print(f"  |numpy - yade| = {abs(phi_np - phi_yade):.3e}, "
      f"|loop - yade| = {abs(phi_loop - phi_yade):.3e}")
print(f"  utils.porosity() cell volume check: O.cell.volume = {O.cell.volume!r}, "
      f"L^3 = {L**3!r}")

print()
print("=" * 78)
print("Q5 -- span and the collider clearance per alpha (S2)")
print("=" * 78)
for alpha in ALPHA_LADDER:
    print(f"  alpha={alpha:.2f}  span/d_eq={span_factor(alpha):.6f}  "
          f"u={alpha_to_u(alpha):.3f}  ds/d_eq={sub_sphere_diameter(1.0, alpha):.6f}")
print("=== probe done ===")
PYEOF
echo "PROBE_CLUMP_S0 DONE rc=$?"
