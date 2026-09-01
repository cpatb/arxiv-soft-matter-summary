"""diag_rattler.py -- why is the rattler fraction coming back 1.000?"""
import numpy as np
from measure import iterated_rattlers as ir_meas
from packing import build_and_compress, contact_graph

st = build_and_compress(cv=0.0, mu_p=0.3, p_hat=1e-3, n=500, seed=1, verbose=False)
from yade import O  # noqa: E402  (Yade must be imported after the scene is built)
counts, neigh, forces = contact_graph(O)
z = np.array(list(counts.values()))
print(f"bodies={len(counts)}  real contacts={sum(counts.values())//2}")
print(f"Z: mean={z.mean():.3f} min={z.min()} max={z.max()} "
      f"frac(Z<4)={np.mean(z<4):.3f} frac(Z==0)={np.mean(z==0):.3f}")
print(f"phi={st.phi:.4f} p_hat={st.p_hat:.3e} n_contacts={st.n_contacts}")
r_meas = ir_meas(counts, neigh)
print(f"iterated_rattlers  measure.py: {len(r_meas)}/{len(counts)}")
one_pass = {b for b,c in counts.items() if c < 4}
print(f"one-pass          : {len(one_pass)}/{len(counts)}")
