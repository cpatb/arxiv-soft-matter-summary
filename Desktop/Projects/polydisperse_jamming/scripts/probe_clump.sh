#!/bin/bash
# probe_clump.sh -- CAPABILITY PROBE ONLY. Runs no production cell, writes no shard.
#
# scripts/probe_shape.sh established that a clump can be MADE. That is not the same as
# "the clump route works in THIS protocol", and four specific things can break between
# the two. Each is asked here, once, on a tiny scene, BEFORE a campaign is designed:
#
#   1. does a clump survive a PERIODIC cell with PeriTriaxController (the protocol),
#      and does the collider accept it with allowBiggerThanPeriod = False?
#   2. does `utils.porosity()` DOUBLE-COUNT the overlap between clump members? If it
#      sums member sphere volumes, then at alpha < 2 the members overlap and every phi
#      the project's central observable reports would be too LARGE, silently. This is
#      the one that would have poisoned the whole axis.
#   3. do clump members interact with THEMSELVES (an internal contact would add fake
#      contacts to n_contacts, z and the fabric)?
#   4. what do mass and inertia come out as -- the same double-count question, on the
#      dynamics rather than on the measurement.
#
# Sourced BEFORE any `set -e` (lessons.md: sourcing under set -e dies silently).
source "$HOME/pj/scripts/env.sh"
set -u

cp /mnt/c/Users/compl/Desktop/Projects/polydisperse_jamming/shape.py ~/pj/shape.py
tr -d '\r' < ~/pj/shape.py > ~/pj/shape.tmp && mv ~/pj/shape.tmp ~/pj/shape.py
cd ~/pj || exit 1

python3 - <<'PYEOF'
import math
import sys

import numpy as np
import yade
from yade import O, utils

sys.path.insert(0, ".")
from shape import clump_centres, corrugation, span_factor, sub_sphere_diameter

ALPHA = 1.5
D_EQ = 1.0
N = 60

print(f"=== probe: alpha = {ALPHA}, d_eq = {D_EQ}, N = {N} clumps ===")
ds = sub_sphere_diameter(D_EQ, ALPHA)
print(f"sub-sphere diameter {ds:.6f}, span/d_eq {span_factor(ALPHA):.6f}, "
      f"corrugation {corrugation(D_EQ, ALPHA):.6f}")

O.reset()
O.periodic = True
L = 8.0
O.cell.hSize = yade.Matrix3(L, 0, 0, 0, L, 0, 0, 0, L)
mat = O.materials.append(utils.FrictMat(young=1.0, poisson=0.3, frictionAngle=0.3,
                                        density=1.0))

rng = np.random.default_rng(0)
axes = rng.normal(size=(N, 3))
axes /= np.linalg.norm(axes, axis=1, keepdims=True)
made = 0
for i in range(N):
    c = rng.uniform(0.1 * L, 0.9 * L, size=3)
    off = clump_centres(D_EQ, ALPHA, axes[i])
    ids = [O.bodies.append(utils.sphere(tuple(c + off[j]), ds / 2.0, material=mat))
           for j in range(2)]
    try:
        O.bodies.clump(ids)
        made += 1
    except Exception as e:
        print("CLUMP FAILED:", type(e).__name__, e)
        break
print(f"1. clumps made in a PERIODIC cell: {made}/{N}")

nclump = sum(1 for b in O.bodies if b is not None and b.isClump)
nmemb = sum(1 for b in O.bodies if b is not None and b.isClumpMember)
nstand = sum(1 for b in O.bodies if b is not None
             and not b.isClump and not b.isClumpMember)
print(f"   bodies: clumps {nclump}, members {nmemb}, standalone {nstand}")

# --- 2. THE POROSITY QUESTION -------------------------------------------------------
v_true = made * (math.pi / 6.0) * D_EQ ** 3
v_members = made * 2 * (math.pi / 6.0) * ds ** 3
phi_true = v_true / L ** 3
phi_members = v_members / L ** 3
phi_yade = 1.0 - utils.porosity()
print("2. POROSITY")
print(f"   phi from the analytic UNION volume : {phi_true:.8f}")
print(f"   phi if member volumes are SUMMED   : {phi_members:.8f}"
      f"   (ratio {phi_members/phi_true:.6f})")
print(f"   phi from utils.porosity()          : {phi_yade:.8f}")
if abs(phi_yade - phi_true) < 1e-9:
    print("   -> utils.porosity() gives the UNION volume. Safe.")
elif abs(phi_yade - phi_members) < 1e-9:
    print("   -> utils.porosity() DOUBLE-COUNTS the member overlap. `packing.measure`")
    print("      must compute phi from the analytic union volume for clump runs.")
else:
    print("   -> utils.porosity() matches NEITHER. It must not be used here at all.")

# --- 3. internal contacts ------------------------------------------------------------
O.engines = [
    utils.ForceResetter(),
    utils.InsertionSortCollider([utils.Bo1_Sphere_Aabb()],
                                allowBiggerThanPeriod=False),
    utils.InteractionLoop(
        [utils.Ig2_Sphere_Sphere_ScGeom()],
        [utils.Ip2_FrictMat_FrictMat_FrictPhys()],
        [utils.Law2_ScGeom_FrictPhys_CundallStrack()]),
    utils.NewtonIntegrator(damping=0.3),
]
O.dt = 1e-3
try:
    O.run(50, True)
    print(f"3. 50 steps ran. O.iter = {O.iter}")
except Exception as e:
    print("3. RUN FAILED:", type(e).__name__, e)

clump_of = {}
for b in O.bodies:
    if b is not None and b.isClumpMember:
        clump_of[b.id] = b.clumpId
internal = 0
real = 0
for i in O.interactions:
    if not i.isReal:
        continue
    real += 1
    if clump_of.get(i.id1, -1) == clump_of.get(i.id2, -2):
        internal += 1
print(f"   real interactions {real}, of which INTERNAL to a clump: {internal}")
print("   -> " + ("internal contacts are present and must be filtered out of "
                  "n_contacts, z and the fabric" if internal else
                  "no internal contacts; the contact graph is clean as written"))

# --- 4. mass and inertia --------------------------------------------------------------
cl = [b for b in O.bodies if b is not None and b.isClump]
if cl:
    m = cl[0].state.mass
    m_union = (math.pi / 6.0) * D_EQ ** 3        # density 1
    m_members = 2 * (math.pi / 6.0) * ds ** 3
    print(f"4. clump mass {m:.8f};  union {m_union:.8f};  summed members {m_members:.8f}"
          f"   (ratio to union {m/m_union:.6f})")
print("=== probe done ===")
PYEOF
