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
    store = {}
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
        store[(fam, mu, cv)] = (phi, out["a2"][0], out["a2"][1], out["C"][0])
    # d ln phi = d ln phi_B - d ln(1 - psi): backbone vs rattler-loading share of each rise (a^2 psi), and the
    # implied change of chi = phi(mu)/phi(0) between a frictional arm and the frictionless arm of one family.
    say("\nsplit of d ln phi (CV 0 -> CV) into backbone d ln phi_B and rattler loading -d ln(1-psi), a^2 psi")
    say("  family      mu   cv   dln phi   backbone  rattler   | d ln chi vs mu 0 (backbone part, rattler part)")
    for (fam, mu, cv), (phi, ps, pb, psC) in sorted(store.items()):
        if cv == 0.0 or (fam, mu, 0.0) not in store:
            continue
        phi0, ps0, pb0, psC0 = store[(fam, mu, 0.0)]
        dl, db, dr = np.log(phi / phi0), np.log(pb / pb0), -np.log((1 - ps) / (1 - ps0))
        extra = ""
        if mu > 0 and (fam, 0.0, cv) in store and (fam, 0.0, 0.0) in store:
            f1, s1, b1, c1 = store[(fam, 0.0, cv)]
            f0, s0, b0, c0 = store[(fam, 0.0, 0.0)]
            dl0, db0, dr0 = np.log(f1 / f0), np.log(b1 / b0), -np.log((1 - s1) / (1 - s0))
            drC = -np.log((1 - psC) / (1 - psC0)) + np.log((1 - c1) / (1 - c0))
            extra = f"| {dl - dl0:+.4f} ({db - db0:+.4f}, {dr - dr0:+.4f}) | C-closure rattler part {drC:+.4f}"
        say(f"  {fam:10s} {mu:.1f} {cv:.2f}  {dl:+.4f}   {db:+.4f}   {dr:+.4f}   {extra}")
    with open(os.path.join(ROOT, "out", "backbone_psi_compare.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
