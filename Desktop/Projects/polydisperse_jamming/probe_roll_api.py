"""probe_roll_api.py — what does THIS Yade build actually expose for rolling resistance?

Written rather than recalled: "the API cannot do X" usually means "the API I tried cannot
do X", and the mirror is just as bad -- writing an incantation from memory that this build
does not have. Run in WSL under scripts/env.sh; prints and exits.
"""
from yade import utils


def show(name, obj, want):
    have = [a for a in dir(obj) if not a.startswith("_")]
    print(f"\n--- {name} ---")
    for w in want:
        print(f"   {w:28s} {'PRESENT' if w in have else 'ABSENT'}"
              + (f"   default = {getattr(obj, w)!r}" if w in have else ""))


m = utils.CohFrictMat()
show("CohFrictMat", m, ["etaRoll", "etaTwist", "momentRotationLaw", "normalCohesion",
                        "shearCohesion", "frictionAngle", "young", "poisson",
                        "alphaKr", "alphaKtw", "isCohesive"])

ip = utils.Ip2_CohFrictMat_CohFrictMat_CohFrictPhys()
show("Ip2_CohFrictMat_CohFrictMat_CohFrictPhys", ip,
     ["setCohesionNow", "setCohesionOnNewContacts", "etaRoll", "etaTwist"])

law = utils.Law2_ScGeom6D_CohFrictPhys_CohesionMoment()
show("Law2_ScGeom6D_CohFrictPhys_CohesionMoment", law,
     ["always_use_moment_law", "useIncrementalForm", "creep_viscosity"])

fm = utils.FrictMat()
show("FrictMat (the production material, for comparison)", fm,
     ["frictionAngle", "young", "poisson", "density"])

print("\n--- what the rolling-resistance path would need ---")
print("   CohFrictMat(normalCohesion=0, shearCohesion=0, momentRotationLaw=True,")
print("               etaRoll=<eta>, alphaKr=<kr>)")
print("   + Ig2_Sphere_Sphere_ScGeom6D, Ip2_CohFrictMat_..., Law2_ScGeom6D_...")
print("   TRAP: with zero cohesion nothing is BONDED, so if the moment law only acts on")
print("   bonded contacts the knob is INERT. `always_use_moment_law` is the flag that")
print("   decides that, and its presence/default is printed above.")
