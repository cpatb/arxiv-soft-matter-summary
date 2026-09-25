"""zbysize_dump.py -- per-grain (diameter, contact count, rattler flag) for one isotropic cell.

Tests derivations/rattler_closure_2026-09-24.md eq. 3 DIRECTLY (z(a) against a(a m1 + m2)), which no
existing record allows: every record stores only pack-level aggregates. The packing is built by the
production `packing.build_and_compress` with the SAME arguments `phase2_phiJ.py` uses (lognormal, default
truncation, N, seed, p_hat, cap from PJ_MAX_ITER), so its phi must reproduce the existing phase-2 record
of the same key (printed as a wiring check, not assumed).

    yade-batch-free: python3 zbysize_dump.py --cv 0.5 --mu 0.3 --seed 1 --phat 3e-5 --n 1500
    -> out/zbysize/zb_cv{cv}_mu{mu}_s{seed}_p{phat}.npz  (d, k, rattler, phi, key)
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np

import packing
from measure import iterated_rattlers
from packing import build_and_compress, contact_graph


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cv", type=float, required=True)
    ap.add_argument("--mu", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--phat", type=float, required=True)
    ap.add_argument("--n", type=int, default=1500)
    a = ap.parse_args()
    st = build_and_compress(cv=a.cv, mu_p=a.mu, p_hat=a.phat, n=a.n, seed=a.seed, family="lognormal")
    from yade import O  # noqa: PLC0415 - WSL only
    counts, neigh, _ = contact_graph(O)
    ids = sorted(counts)
    d = np.array([2.0 * float(O.bodies[i].shape.radius) for i in ids])
    k = np.array([counts[i] for i in ids], dtype=int)
    rat = iterated_rattlers(counts, neigh)
    r = np.array([i in rat for i in ids], dtype=bool)
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out", "zbysize")
    os.makedirs(out, exist_ok=True)
    rec = json.loads(st.to_json())
    fn = os.path.join(out, f"zb_cv{a.cv}_mu{a.mu}_s{a.seed}_p{a.phat:g}.npz")
    pos = {i: n for n, i in enumerate(ids)}
    pairs = np.array(sorted({(min(pos[i], pos[j]), max(pos[i], pos[j])) for i in ids for j in neigh[i]}),
                     dtype=int).reshape(-1, 2)
    assert 2 * len(pairs) == int(k.sum()), (len(pairs), int(k.sum()))
    np.savez(fn, d=d, k=k, rattler=r, pairs=pairs, phi=rec["phi"], converged=rec["converged"],
             z_all=rec.get("z_all", np.nan), rattler_frac=rec.get("rattler_frac", np.nan),
             max_iter=packing.max_iter())
    print(json.dumps({"file": fn, "phi": rec["phi"], "converged": rec["converged"], "N": len(ids),
                      "z_all_rec": rec.get("z_all"), "z_all_here": float(k.mean()),
                      "fr_rec": rec.get("rattler_frac"), "fr_here": float(r.mean())}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
