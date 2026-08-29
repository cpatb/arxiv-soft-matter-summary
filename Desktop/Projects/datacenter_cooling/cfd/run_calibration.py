"""Sweep the CFD recirculation fraction phi over containment eta (and two heat loads),
saving each converged point to cfd_recirc.json as it is computed (timeout-safe)."""
import json, os, numpy as np
from aisle import AisleCell

OUT = "cfd_recirc.json"
results = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(r["eta"], r["Ra"]) for r in results}

ETAS = [0.0, 0.2, 0.4, 0.6, 0.8, 0.95]
RAS = [4e3, 8e3]          # two rack heat loads (Ra ~ dissipated power)
N = 40

for Ra in RAS:
    for eta in ETAS:
        if (eta, Ra) in done:
            continue
        a = AisleCell(N=N, H_room=1.0, eta=eta, Ra=Ra)
        h = min(a.dx, a.dy)
        dt = 0.35 * min(h*h/(4*max(a.Pr, 1.0)), h/(np.sqrt(a.Ra*a.Pr)+1.0))
        prev, phi = 0.0, 0.0
        for n in range(1, 40001):
            a.step(dt)
            if n % 1000 == 0:
                phi = a.recirc_fraction()
                if n >= 6000 and abs(phi - prev) < 8e-4:
                    break
                prev = phi
        rec = dict(eta=eta, Ra=Ra, phi=round(phi, 4),
                   phi_tracer=round(a.recirc_fraction_tracer(), 4), steps=n)
        results.append(rec)
        json.dump(results, open(OUT, "w"), indent=2)
        print(f"Ra={Ra:.0e} eta={eta:.2f}: phi={phi:.3f}  ({n} steps)")
print("done ->", OUT)
