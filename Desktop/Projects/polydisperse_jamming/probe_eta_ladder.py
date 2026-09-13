"""probe_eta_ladder.py — how far does rolling resistance move phi, and where does it
saturate?

Scoping only, before PREREG_ROLLING_2026-09-13.md fixes a grid. Friction saturates
(phi_J moves 0.016 between mu_p 1 and 3 and then stops), so the first thing to know about
the new knob is whether it does too, and over what range of `etaRoll` the lever is real.

phi at a single p_hat = 3e-4, not an extrapolated phi_J: this is a SCOPING run and the
number it reports is a proxy. One seed per rung, N = 600, a few seconds each.
"""
from __future__ import annotations

import json
import os

from packing import build_and_compress

ETAS = [None, 0.0, 0.05, 0.1, 0.3, 1.0, 3.0, 10.0]
CVS = [0.0, 0.5]
SEED = 1
N = 600
P_HAT = 3e-4
MU_P = 0.3

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")


def main() -> int:
    rows = []
    for cv in CVS:
        for eta in ETAS:
            st = build_and_compress(cv=cv, mu_p=MU_P, p_hat=P_HAT, n=N, seed=SEED,
                                    eta_roll=eta, verbose=False)
            d = json.loads(st.to_json())
            rows.append({"cv": cv, "eta": eta, "phi": d["phi"], "z_nr": d["z_nr"],
                         "rat": d["rattler_frac"], "conv": d["converged"],
                         "wall": d["wall_s"]})
            print(f"cv {cv}  eta {str(eta):>5}  phi {d['phi']:.4f}  "
                  f"z_nr {d['z_nr']:.2f}  rat {d['rattler_frac']:.3f}  "
                  f"conv {d['converged']}  {d['wall_s']:.1f}s", flush=True)

    print("\n=== the lever ===")
    for cv in CVS:
        rr = [r for r in rows if r["cv"] == cv and r["conv"]]
        base = [r for r in rr if r["eta"] == 0.0]
        if not base:
            continue
        b = base[0]["phi"]
        print(f"cv {cv}: phi(etaRoll=0) = {b:.4f}")
        for r in rr:
            if r["eta"] in (None, 0.0):
                continue
            print(f"        eta {r['eta']:>5}: phi {r['phi']:.4f}  "
                  f"({r['phi'] - b:+.4f})")
    print("\n=== the CV contrast at each eta (what the campaign measures) ===")
    for eta in ETAS:
        a = [r for r in rows if r["cv"] == 0.0 and r["eta"] == eta and r["conv"]]
        c = [r for r in rows if r["cv"] == 0.5 and r["eta"] == eta and r["conv"]]
        if a and c:
            print(f"   eta {str(eta):>5}: phi(CV0) {a[0]['phi']:.4f}  "
                  f"phi(CV0.5) {c[0]['phi']:.4f}  dphi {c[0]['phi'] - a[0]['phi']:+.4f}")
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "eta_ladder_probe.json"), "w", encoding="utf-8") as fh:
        json.dump(rows, fh, indent=1)
    print("\n⚠ ONE SEED, and phi at a single pressure rather than an extrapolated phi_J. "
          "Scoping only.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
