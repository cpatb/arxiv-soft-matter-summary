"""diag_kladder_descriptor.py -- what single quantity, if any, organises the K-ladder
residual?

**POST-HOC. NOT A PRE-REGISTERED GATE.** Written after `grade_kladder.py` ran, and it
says so wherever it prints. Nothing here may be quoted as a blind prediction, and the
registered verdicts in `out/grade_kladder.txt` do not depend on a line of it.

WHAT IT DOES. `grade_kladder.py` writes one row per admitted cell to
`out/kladder_cells.csv` (label, K, y_S, measured phi, the additive prediction, the
residual, and the finite-size diagnostics). This script ranks candidate descriptors by
how well each orders that residual, and it deliberately includes the two NULL candidates
that would mean the result is an artifact:

  * `n_min`, the smallest design class count -- a residual rank-correlated with a
    sample-size variable is a finite-size artifact until a second N says otherwise
    (lessons.md 2026-09-14). Both N are pooled here and the split is printed.
  * `d_max_over_L`, the periodic-cell finite-size diagnostic.

Any descriptor that beats those two is reported with the null's own number beside it, so
a reader can see the comparison rather than take the ranking on trust.

WHAT IT CANNOT DO. Ranking descriptors on the same 18 cells that motivated the ranking
is a search, not a test. A descriptor that wins here is a HYPOTHESIS for the next block,
and the file says so in its own output.
"""

from __future__ import annotations

import csv
import os

import numpy as np
from scipy import stats

import kladder as kl

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")
CSV = os.path.join(OUT, "kladder_cells.csv")


def descriptors(d_rel, y):
    """Computable summaries of a discrete PSD, each a candidate control parameter."""
    p = kl.psd_for(d_rel, y)
    d = np.asarray(p.diameters, float)
    yv = np.asarray(p.volume_fractions, float)
    k = d.size
    pairs = kl.pair_table(d)
    # volume-product weight per pair, and the same weighted by how strongly a pair
    # fills (1 - u is the size contrast; u -> 0 is a perfectly filling pair)
    wsum = 0.0
    fill = 0.0
    for di, dj, u in pairs:
        i = int(np.argmin(np.abs(d - di)))
        j = int(np.argmin(np.abs(d - dj)))
        w = float(yv[i] * yv[j])
        wsum += w
        fill += w * (1.0 - u)
    trip = 0.0
    for a in range(k):
        for b in range(a + 1, k):
            for c in range(b + 1, k):
                trip += float(yv[a] * yv[b] * yv[c])
    return {
        "K": float(k),
        "n_pairs": float(k * (k - 1) / 2),
        "n_triples": float(k * (k - 1) * (k - 2) / 6),
        "y_S": float(yv[0]),
        "CV": float(p.cv),
        "d43_over_d32": float(p.d43 / p.d32),
        "pair_weight": wsum,
        "fill_weight": fill,
        "triple_weight": trip,
    }


def main() -> int:
    lines: list[str] = []

    def say(s: str = "") -> None:
        lines.append(s)
        print(s)

    if not os.path.exists(CSV):
        print(f"{CSV} not found -- run grade_kladder.py first")
        return 2
    rows = list(csv.DictReader(open(CSV, encoding="utf-8")))
    say("=== diag_kladder_descriptor.py -- POST-HOC, NOT A PRE-REGISTERED GATE ===")
    say("   Ranking descriptors on the cells that motivated the ranking is a SEARCH.")
    say("   Anything that wins here is a hypothesis for the next block, not a result.")
    say()

    cells = {c["label"]: c for c in kl.build_cells(kl.N_BIG)[0]}
    for n in sorted({int(r["N"]) for r in rows}):
        sub = [r for r in rows if int(r["N"]) == n and r["arm"] in ("K", "SPLIT")]
        if len(sub) < 5:
            continue
        R = np.array([float(r["R"]) for r in sub])
        say(f"=== N = {n}: {len(sub)} cells (arms K and SPLIT) ===")
        say(f"   {'descriptor':>16}{'spearman':>11}{'p':>9}{'pearson':>10}   note")
        cand: dict[str, np.ndarray] = {}
        for r in sub:
            for kname, v in descriptors(cells[r["label"]]["d_rel"],
                                        cells[r["label"]]["y"]).items():
                cand.setdefault(kname, []).append(v)  # type: ignore[union-attr]
        cand = {a: np.asarray(b, float) for a, b in cand.items()}
        cand["n_min"] = np.array([float(r["n_min"]) for r in sub])
        cand["d_max_over_L"] = np.array([float(r["d_max_over_L"]) for r in sub])
        out = {}
        for nm, v in cand.items():
            rho = stats.spearmanr(v, R)
            pr = float(np.corrcoef(v, R)[0, 1])
            out[nm] = (float(rho.statistic), float(rho.pvalue), pr)  # type: ignore[attr-defined]
        for nm in sorted(out, key=lambda x: -abs(out[x][0])):
            rho, p, pr = out[nm]
            note = ("<-- NULL: a sample-size variable" if nm in
                    ("n_min", "d_max_over_L") else "")
            say(f"   {nm:>16}{rho:+11.4f}{p:9.4f}{pr:+10.4f}   {note}")
        say(f"   best NULL |spearman| = "
            f"{max(abs(out['n_min'][0]), abs(out['d_max_over_L'][0])):.4f}; any "
            f"descriptor below that has explained nothing the sample size does not.")
        say()

    # pooled across both N, which is the only way the sample-size nulls can be broken
    allsub = [r for r in rows if r["arm"] in ("K", "SPLIT")]
    ns = sorted({int(r["N"]) for r in allsub})
    if len(ns) > 1:
        say(f"=== pooled over both N ({len(allsub)} cells) -- the split that breaks the "
            "sample-size nulls ===")
        R = np.array([float(r["R"]) for r in allsub])
        for nm in ("n_min", "d_max_over_L"):
            v = np.array([float(r[nm]) for r in allsub])
            rho = stats.spearmanr(v, R)
            say(f"   spearman(R, {nm}) pooled = {float(rho.statistic):+.4f} "  # type: ignore[attr-defined]
                f"(p {float(rho.pvalue):.4f})")  # type: ignore[attr-defined]
        say("   Each cell appears at BOTH N with identical composition, so a descriptor")
        say("   that survives pooling is not reading the grain count.")
        say()
        say("   per-cell N sensitivity (R at N = 3000 minus R at N = 1500):")
        by = {}
        for r in allsub:
            by.setdefault(r["label"], {})[int(r["N"])] = float(r["R"])
        ds = []
        for lbl, v in sorted(by.items()):
            if len(v) == 2:
                d = v[ns[-1]] - v[ns[0]]
                ds.append(d)
                say(f"      {lbl:>14}  {v[ns[0]]:+7.2f} -> {v[ns[-1]]:+7.2f}   "
                    f"{d:+7.2f}")
        if ds:
            say(f"   mean shift {float(np.mean(ds)):+.2f}, rms "
                f"{float(np.sqrt(np.mean(np.square(ds)))):.2f} scale over "
                f"{len(ds)} cells")

    say()
    say("=== what this is worth ===")
    say("   A descriptor that orders the residual here has NOT been tested; it has been")
    say("   selected. The test is a block designed on it, with cells placed where the")
    say("   descriptor and its rivals DISAGREE -- which is the same move that turned the")
    say("   bidisperse 'no descriptor works' negative into the filling functional.")
    with open(os.path.join(OUT, "diag_kladder_descriptor.txt"), "w",
              encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
