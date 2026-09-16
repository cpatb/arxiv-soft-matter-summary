#!/bin/bash
# probe_shape.sh -- CAPABILITY PROBE ONLY. Runs no cell and writes no shard.
#
# Answers one question before any shape-axis design is written: what does this Yade
# install actually support for NON-SPHERICAL particles? The project is spheres-only
# today, and "add shape" means one of three very different amounts of work depending
# on the answer:
#
#   clumps present      -> rigid sphere-aggregates, aspect ratio by construction.
#                          Cheapest real shape axis; works with the existing
#                          Ig2_Sphere_Sphere_ScGeom contact law untouched.
#   PotentialParticles  -> true convex shapes, different contact law, big change.
#   neither             -> shape is out of reach on this install.
#
# Sourced BEFORE any `set -e`, because lessons.md records that sourcing a bashrc under
# `set -e` kills the script with no output and exit 1.

source "$HOME/pj/scripts/env.sh"
set -u

python3 - <<'PYEOF'
import sys
try:
    import yade
    from yade import utils, wrapper
except Exception as e:
    print("IMPORT FAILED:", type(e).__name__, e)
    sys.exit(1)

print("yade version:", getattr(yade, "version", "?"))

# 1. clumps -- the cheap route
have = {}
for name in ("clump", "clumpCloud"):
    have[name] = hasattr(utils, name)
o = yade.O
for name in ("addClump", "clumps"):
    have["O." + name] = hasattr(o, name) or hasattr(o.bodies, name)
print("clump API:", have)

# 2. non-spherical shape CLASSES actually compiled in
shapes = []
for name in ("Sphere", "Facet", "Box", "Clump", "PotentialParticle", "Polyhedra",
             "Ellipsoid", "Level Set", "LevelSet", "Tetra", "GridNode"):
    if hasattr(wrapper, name.replace(" ", "")):
        shapes.append(name)
print("shape classes present:", shapes)

# 3. the contact laws that would be needed for each
laws = [n for n in dir(wrapper)
        if n.startswith(("Ig2_", "Law2_")) and
        any(k in n for k in ("Polyhedra", "Potential", "Ellipsoid", "LevelSet", "Tetra"))]
print("non-sphere contact laws:", laws[:12] if laws else "NONE")

# 4. can a clump actually be MADE? the only question that decides the cheap route.
try:
    yade.O.reset()
    from yade import wrapper as w
    ids = [yade.O.bodies.append(utils.sphere((0, 0, 0), 0.5)),
           yade.O.bodies.append(utils.sphere((0.6, 0, 0), 0.5))]
    cid = yade.O.bodies.clump(ids)
    b = yade.O.bodies[cid]
    print("CLUMP MADE: id", cid, "isClump", b.isClump,
          "members", len(b.shape.members) if hasattr(b.shape, "members") else "?")
    print("-> the cheap shape route is AVAILABLE")
except Exception as e:
    print("CLUMP CREATION FAILED:", type(e).__name__, e)
    print("-> the cheap shape route is NOT available as written")
PYEOF
