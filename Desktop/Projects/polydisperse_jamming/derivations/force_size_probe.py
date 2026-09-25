"""force_size_probe.py -- EXPLORATORY: how the Derjaguin-reduced contact force g = f_n/(R* p) depends on the
pair size product a_i a_j in the sheared actshape dumps, and whether a single exponent beta in
<g | aa> ~ (aa)^beta reproduces the measured force-size correlation C = <g aa>/(<g><aa>).

    python derivations/force_size_probe.py  ->  out/force_size_probe.txt

beta per cell = OLS slope of ln g on ln(a_i a_j) over all contacts of the 10 snapshots (fn > 0). C_pred =
<(aa)^{1+beta}>/(<(aa)^beta><aa>) over the same contacts (so it tests only the power-law form, not the pair
statistics). CV 0 cells have aa constant and are skipped.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import actshape as A  # noqa: E402

lines: list[str] = []


def say(s: str = "") -> None:
    print(s)
    lines.append(s)


def main() -> int:
    say("mu   p       CV   seed  N     beta    C_meas   C_pred(beta)   C_pred(beta=0)")
    rows = []
    for f in sorted(glob.glob(os.path.join(ROOT, "out", "phase9_actshape_*.jsonl"))):
        for line in open(f, encoding="utf-8"):
            r = json.loads(line)
            if "error" in r or r["cv"] == 0.0 or not (r["steady"] and r["servo_held"] and r["quasistatic"]):
                continue
            path = os.path.join(ROOT, "out", "actshape_dump", r["key"] + ".npz")
            if not os.path.exists(path):
                continue
            d = np.load(path)
            E = float(d["young"])
            edges = np.concatenate([[0], np.cumsum(d["n"])])
            gs, aas = [], []
            for k in range(len(d["n"])):
                sl = slice(edges[k], edges[k + 1])
                fn, r1, r2 = d["fn"][sl], d["r1"][sl], d["r2"][sl]
                ok = fn > 0
                gs.append(A.g_values(fn, r1, r2, float(d["p_hat"][k]) * E, "derjaguin"))
                aas.append((r1 * r2)[ok])
            g, aa = np.concatenate(gs), np.concatenate(aas)
            beta = float(np.polyfit(np.log(aa), np.log(g), 1)[0])
            C = float(np.mean(g * aa) / (g.mean() * aa.mean()))
            Cp = float(np.mean(aa ** (1 + beta)) / (np.mean(aa ** beta) * aa.mean()))
            rows.append((r["mu_p"], r["p_hat_target"], r["cv"], beta, C, Cp))
            say(f"{r['mu_p']:.1f}  {r['p_hat_target']:.0e}  {r['cv']:.1f}  {r['seed']:3d}  {r['n_placed']:5d}  "
                f"{beta:+.4f}  {C:.4f}   {Cp:.4f}         1.0000")
    b = np.array([x[3] for x in rows])
    say(f"\nbeta over {len(rows)} cells: mean {b.mean():+.4f}, sd {b.std(ddof=1):.4f}, range [{b.min():+.4f}, {b.max():+.4f}]")
    for cv in (0.2, 0.5):
        for mu in (0.0, 0.3):
            s = [x for x in rows if x[2] == cv and x[0] == mu]
            if s:
                say(f"  CV {cv} mu {mu}: beta {np.mean([x[3] for x in s]):+.4f}  C_meas {np.mean([x[4] for x in s]):.4f}  "
                    f"C_pred {np.mean([x[5] for x in s]):.4f}")
    with open(os.path.join(ROOT, "out", "force_size_probe.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
