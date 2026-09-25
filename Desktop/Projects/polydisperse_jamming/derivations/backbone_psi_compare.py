"""backbone_psi_compare.py -- rattler_closure_2026-09-24.md §4b under both lambda shapes (C and a^2).

    python derivations/backbone_psi_compare.py  ->  out/backbone_psi_compare.txt

phi_B = phi (1 - psi), psi = <a^3 P(k<2 | lambda(a))>/<a^3>, lambda normalised to each record's z_all,
isotropic phase2 records at the p_hat 3e-5 rung, seed means per (family, mu, CV).
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
from scipy.stats import poisson

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, ROOT)
from rattler_closure_check import radii_family  # noqa: E402

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def shape(a, which):
    m1, m2 = a.mean(), (a * a).mean()
    s = a * (a * m1 + m2) if which == "C" else a * a
    return s / s.mean()


def psi(a, Z, which):
    pr = poisson.cdf(1, Z * shape(a, which))
    return float(np.sum(a**3 * pr) / np.sum(a**3))


def main() -> int:
    recs = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase2*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if ("error" in r or not r.get("converged") or r.get("family") not in ("lognormal", "uniform", "bidisperse")
                    or r.get("alpha", 1.0) != 1.0 or r.get("doublet_w", 0) or r.get("triplet_w", 0)
                    or "z_all" not in r or abs(r.get("p_hat_target", 0) - 3e-5) > 1e-9):
                continue
            recs[r["key"]] = r
    groups = {}
    for r in recs.values():
        groups.setdefault((r["family"], r["mu_p"], round(r["cv"], 3)), []).append(r)
    say("backbone phi_B = phi (1 - psi_model), p_hat 3e-5, seed means; d = change vs CV 0 of the same (family, mu)")
    say("  family      mu   cv    n   phi     dphi    | C: psi   dphi_B  | a2: psi  dphi_B")
    base = {}
    for (fam, mu, cv), rs in sorted(groups.items()):
        a = radii_family(fam, cv)
        phi = np.mean([r["phi"] for r in rs])
        out = {}
        for w in ("C", "a2"):
            ps = np.mean([psi(a, r["z_all"], w) for r in rs])
            pb = np.mean([r["phi"] * (1 - psi(a, r["z_all"], w)) for r in rs])
            out[w] = (ps, pb)
        if cv == 0.0:
            base[(fam, mu)] = (phi, out["C"][1], out["a2"][1])
        b = base.get((fam, mu), (np.nan, np.nan, np.nan))
        say(f"  {fam:10s} {mu:.1f} {cv:.2f} {len(rs):4d}  {phi:.4f}  {phi - b[0]:+.4f}  | {out['C'][0]:.4f}  "
            f"{out['C'][1] - b[1]:+.4f} | {out['a2'][0]:.4f}  {out['a2'][1] - b[2]:+.4f}")
    with open(os.path.join(ROOT, "out", "backbone_psi_compare.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
