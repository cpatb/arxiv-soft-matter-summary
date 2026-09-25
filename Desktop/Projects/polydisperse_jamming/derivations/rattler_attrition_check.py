"""rattler_attrition_check.py -- convergence attrition of the isotropic phase2 records by (family, mu, CV),
and whether the dropped (non-converged) records rattle more. Deduplicated by key (last record wins, as in
rattler_closure_check.isotropic_block). Same family/alpha/multiplet filters, WITHOUT the converged filter.

    python derivations/rattler_attrition_check.py  ->  out/rattler_attrition.txt
"""
from __future__ import annotations

import glob
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def main() -> int:
    recs = {}
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase2*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            if ("error" in r or r.get("family") not in ("lognormal", "uniform", "bidisperse")
                    or r.get("alpha", 1.0) != 1.0 or r.get("doublet_w", 0) or r.get("triplet_w", 0)
                    or "z_all" not in r or "rattler_frac" not in r):
                continue
            recs[r["key"]] = r
    rs = list(recs.values())
    say(f"records after dedup: {len(rs)}")
    say("by CV (all families, mu):")
    for cv in sorted({round(r["cv"], 3) for r in rs}):
        g = [r for r in rs if round(r["cv"], 3) == cv]
        c = sum(bool(r.get("converged")) for r in g)
        say(f"  CV {cv:.2f}: {c}/{len(g)} converged ({100 * c / len(g):.1f} %)")
    say("by (family, CV >= 0.4): converged / total, mean f_r converged vs not")
    for fam in ("lognormal", "uniform", "bidisperse"):
        for cv in (0.4, 0.5):
            g = [r for r in rs if r["family"] == fam and round(r["cv"], 3) == cv]
            if not g:
                continue
            cc = [r["rattler_frac"] for r in g if r.get("converged")]
            nc = [r["rattler_frac"] for r in g if not r.get("converged")]
            say(f"  {fam:10s} CV {cv}: {len(cc)}/{len(g)} converged; f_r {np.mean(cc) if cc else float('nan'):.4f} "
                f"(conv) vs {np.mean(nc) if nc else float('nan'):.4f} (not conv, n {len(nc)})")
    with open(os.path.join(ROOT, "out", "rattler_attrition.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
