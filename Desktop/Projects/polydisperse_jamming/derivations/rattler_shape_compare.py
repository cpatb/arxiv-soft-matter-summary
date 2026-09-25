"""rattler_shape_compare.py -- which size dependence of lambda(a) predicts rattler COUNTS best?

    python derivations/rattler_shape_compare.py  ->  out/rattler_shape_compare.txt

Same records, same Poisson rule, same z_all normalisation as rattler_closure_check.py closure (i); only the
shape of lambda(a) changes: C = a(a m1 + m2), a^1, a^2, a^3, size-blind. Added after the per-grain test
(out/zbysize_check.txt) preferred a^2 over C in one cell.
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
from rattler_closure_check import radii, radii_family  # noqa: E402

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def shape(a, which):
    m1, m2 = a.mean(), (a * a).mean()
    s = {"C": a * (a * m1 + m2), "a1": a, "a2": a * a, "a3": a ** 3, "blind": np.ones_like(a)}[which]
    return s / s.mean()


WH = ("C", "a1", "a2", "a3", "blind")


def fr(a, Z, which):
    return float(np.mean(poisson.cdf(1, Z * shape(a, which))))


def main() -> int:
    sheared = []
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase9_actshape_*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if "error" not in r and r["steady"] and r["servo_held"] and r["quasistatic"]:
                sheared.append((radii(r["cv"]), r["z_ss"], r["rattler_frac_ss"], r["cv"]))
    iso = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase2*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if ("error" in r or not r.get("converged") or r.get("family") not in ("lognormal", "uniform", "bidisperse")
                    or r.get("alpha", 1.0) != 1.0 or r.get("doublet_w", 0) or r.get("triplet_w", 0)
                    or "z_all" not in r or "rattler_frac" not in r):
                continue
            iso[r["key"]] = r
    isol = [(radii_family(r["family"], round(r["cv"], 3)), r["z_all"], r["rattler_frac"], r["cv"]) for r in iso.values()]
    for name, data in (("sheared actshape (52)", sheared), (f"isotropic phase2 ({len(isol)})", isol)):
        say(f"=== {name} ===")
        for w in WH:
            e = np.array([fr(a, Z, w) - f for a, Z, f, _ in data])
            e5 = np.array([fr(a, Z, w) - f for a, Z, f, cv in data if cv >= 0.4])
            say(f"  {w:5s}: rms {np.sqrt(np.mean(e**2)):.4f} mean {e.mean():+.4f} | CV>=0.4: rms "
                f"{np.sqrt(np.mean(e5**2)):.4f} mean {e5.mean():+.4f}")
    with open(os.path.join(ROOT, "out", "rattler_shape_compare.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
