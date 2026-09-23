"""Direct geometric check of the inferred deep overlaps (run in WSL from ~/pjdt).
Builds the bad cell (mu 3, w2 = 1, p 1e-3, seed 1), then reads i.geom.penetrationDepth and
i.phys.kn per inter-grain contact, and compares fn/kn with the geometric depth."""
import numpy as np
from yade import O
import packing

st = packing.build_and_compress(cv=0.0, mu_p=3.0, p_hat=1e-3, n=1700, seed=1, verbose=False,
                                alpha=1.0, doublet_w=1.0)
print("phi", st.phi, "conv", st.converged)
gid = {}
for b in O.bodies:
    if b is None or b.isClump:
        continue
    gid[b.id] = b.clumpId if b.isClumpMember else b.id
rows = []
for i in O.interactions:
    if not i.isReal or gid.get(i.id1) == gid.get(i.id2):
        continue
    rows.append((i.phys.normalForce.norm(), i.phys.kn, i.geom.penetrationDepth,
                 i.phys.shearForce.norm()))
a = np.array(rows)
fn, kn, pen, fs = a.T
print("n contacts", len(a), "kn unique", np.unique(np.round(kn, 6))[:5])
rel = pen / 0.5
print("geometric u/r: median %.2e max %.3f frac>0.05 %.4f" % (np.median(rel), rel.max(), (rel > 0.05).mean()))
print("fn/kn vs pen: max |fn/kn - pen| = %.3e" % np.abs(fn / kn - pen).max())
top = np.argsort(fn)[::-1][:5]
for t in top:
    print("  fn %.4f kn %.3f pen %.4f (u/r %.3f) fs/fn %.3f" % (fn[t], kn[t], pen[t], pen[t] / 0.5,
                                                               fs[t] / max(fn[t], 1e-30)))
