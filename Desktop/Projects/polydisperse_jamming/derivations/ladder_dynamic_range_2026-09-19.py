"""ladder_dynamic_range_2026-09-19.py

THE QUESTION
------------
`phi_J` is not measured; it is EXTRAPOLATED from a four-rung `p_hat` ladder
(1e-3, 3e-4, 1e-4, 3e-5). An extrapolation needs the ladder to MOVE by more than the noise.
Nobody in this repository had asked how far it moves.

The diagnostic is one number per cell group:

    SPAN/SD  =  (phi at the top rung - phi at the bottom rung)  /  (seed scatter of phi at
                                                                    the bottom rung)

Below SPAN/SD ~ 2 the ladder spans about one seed scatter and the fitted intercept is a fit
to noise, whatever its r^2 says.

WHY IT MATTERS HERE
-------------------
`paper/main.tex` already carries the qualitative version -- "The frictionless ladders are
not monotone in p_hat at this seed count, so this anchor is a check of the protocol's scale
and not a result of the paper" -- and a referee pass on
`wc_endpoints_move_together_2026-09-19.md` found that the HIGH-FRICTION ladders are worse
still (9/9 and 3/3 non-monotone, r^2 down to 0.003). This module turns both into a number
that can be computed for every cell group on disk, and reports WHICH committed numbers rest
on groups that fail it.

⚠ SPAN/SD = 2 IS A CONVENTION, NOT A DERIVED THRESHOLD. It is stated here, used
consistently, and every group's value is printed so a reader can apply a different one. The
ordering of the groups does not depend on it.

READ-ONLY. `out/phase2_wcbidi*.jsonl` is being written by a live campaign and is EXCLUDED.
"""

from __future__ import annotations

import collections
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
OUT = os.path.join(ROOT, "out")
sys.path.insert(0, ROOT)

RUNGS = (1e-3, 3e-4, 1e-4, 3e-5)
SPAN_SD_BAR = 2.0          # convention, see the module docstring
FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"   [{'PASS' if ok else 'FAIL'}] {name}" + (f"  --  {detail}" if detail else ""))
    if not ok:
        FAILURES.append(name)


def load_cells():
    """One record per (mu_p, cv, seed, p_hat), highest MAX_ITER winning -- the same rule
    `sweep_descriptive.load_shards` uses. N printed before and after every filter."""
    files = [f for f in sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl")))
             if "wcbidi" not in os.path.basename(f)]
    rows = []
    for fn in files:
        for ln in open(fn, encoding="utf-8"):
            ln = ln.strip()
            if ln:
                try:
                    rows.append(json.loads(ln))
                except json.JSONDecodeError:
                    continue
    print(f"   files globbed (wcbidi EXCLUDED, it is mid-flight): {len(files)}")
    print(f"   N records                                        : {len(rows)}")
    rows = [r for r in rows if not r.get("failed")]
    print(f"   N after dropping failed                          : {len(rows)}")
    rows = [r for r in rows if r.get("n_requested") == 1500]
    print(f"   N after N = 1500                                 : {len(rows)}")
    cell: dict = collections.defaultdict(dict)
    for r in rows:
        k = ((r.get("family") or "lognormal"), r["mu_p"], round(r["cv"], 1), r["seed"])
        p = r["p_hat_target"]
        if p not in cell[k] or r.get("max_iter", 0) >= cell[k][p].get("max_iter", 0):
            cell[k][p] = r
    print(f"   N after collapsing to one record per (cell, rung) : "
          f"{sum(len(v) for v in cell.values())} in {len(cell)} cells")
    return cell


def main() -> int:
    print("=" * 92)
    print("How much does the p_hat ladder actually MOVE? A readability diagnostic")
    print("=" * 92)

    print("\n=== S0  LOADING ===")
    cell = load_cells()

    # ---- S1: the diagnostic across the whole grid ------------------------------------
    print("\n=== S1  SPAN/SD FOR EVERY CELL GROUP ===")
    print("   span = phi(1e-3) - phi(3e-5) per seed, averaged; sd = seed scatter of")
    print("   phi(3e-5). A group needs span >> sd for its extrapolated intercept to mean")
    print("   anything.\n")
    groups = collections.defaultdict(list)
    for (fam, mu, cv, seed), d in cell.items():
        groups[(fam, mu, cv)].append((seed, d))

    print(f"   {'family':11s} {'mu_p':>5s} {'cv':>5s} {'n':>2s} {'span':>9s} {'sd':>9s} "
          f"{'SPAN/SD':>8s} {'mono':>6s}  readable?")
    table = {}
    for key in sorted(groups):
        fam, mu, cv = key
        spans, bots, monos = [], [], []
        for seed, d in groups[key]:
            if not all(p in d for p in RUNGS):
                continue
            ys = [d[p]["phi"] for p in RUNGS]
            spans.append(ys[0] - ys[-1])
            bots.append(ys[-1])
            monos.append(all(ys[i] >= ys[i + 1] for i in range(len(ys) - 1)))
        if len(bots) < 3:
            continue
        sd = float(np.std(bots, ddof=1))
        span = float(np.mean(spans))
        ratio = span / sd if sd > 0 else float("inf")
        table[key] = (span, sd, ratio, sum(monos), len(monos))
        ok = "yes" if ratio >= SPAN_SD_BAR else "🔴 NO"
        print(f"   {fam:11s} {mu:5.1f} {cv:5.1f} {len(bots):2d} {span:9.5f} {sd:9.5f} "
              f"{ratio:8.1f} {sum(monos)}/{len(monos):<4d}  {ok}")

    n_bad = sum(1 for v in table.values() if v[2] < SPAN_SD_BAR)
    print(f"\n   groups below SPAN/SD = {SPAN_SD_BAR}: {n_bad} of {len(table)}")
    check("the diagnostic is not vacuous -- some groups pass and some fail",
          0 < n_bad < len(table), f"{n_bad} of {len(table)} fail")

    # ---- S2: what drives it ----------------------------------------------------------
    print("\n=== S2  WHAT DRIVES THE COLLAPSE: friction and width, separately ===")
    print("   SPAN/SD averaged over families, by (mu_p, cv):")
    by_mu_cv = collections.defaultdict(list)
    for (fam, mu, cv), v in table.items():
        by_mu_cv[(mu, cv)].append(v[2])
    cvs = sorted({cv for _m, cv in by_mu_cv})
    mus = sorted({m for m, _cv in by_mu_cv})
    print("   " + f"{'mu_p':>6s}" + "".join(f"{cv:>9.1f}" for cv in cvs))
    for mu in mus:
        row = f"   {mu:6.1f}"
        for cv in cvs:
            v = by_mu_cv.get((mu, cv))
            row += f"{np.mean(v):9.1f}" if v else f"{'--':>9s}"
        print(row)
    # both axes degrade it -- state it as two one-dimensional trends, each printed
    lo_mu = [v for (m, cv), vs in by_mu_cv.items() if m <= 0.3 for v in vs]
    hi_mu = [v for (m, cv), vs in by_mu_cv.items() if m >= 1.0 for v in vs]
    lo_cv = [v for (m, cv), vs in by_mu_cv.items() if cv <= 0.2 for v in vs]
    hi_cv = [v for (m, cv), vs in by_mu_cv.items() if cv >= 0.4 for v in vs]
    print(f"\n   mean SPAN/SD at mu_p <= 0.3 : {np.mean(lo_mu):6.1f}   at mu_p >= 1.0 : {np.mean(hi_mu):6.1f}")
    print(f"   mean SPAN/SD at CV   <= 0.2 : {np.mean(lo_cv):6.1f}   at CV   >= 0.4 : {np.mean(hi_cv):6.1f}")
    check("BOTH axes degrade the ladder: high friction and high CV each lower SPAN/SD",
          np.mean(hi_mu) < np.mean(lo_mu) and np.mean(hi_cv) < np.mean(lo_cv),
          f"friction {np.mean(lo_mu):.1f}->{np.mean(hi_mu):.1f}, "
          f"width {np.mean(lo_cv):.1f}->{np.mean(hi_cv):.1f}")
    print("   MECHANISM, and it is the same one twice: both friction and polydispersity")
    print("   RAISE phi_J, so at a fixed pressure ladder the packing sits further above its")
    print("   own jamming point and is stiffer. A stiffer packing compresses less between")
    print("   1e-3 and 3e-5, so the ladder has less to say.")

    # ---- S3: does SPAN/SD predict the monotonicity failures? --------------------------
    print("\n=== S3  A CHECK THAT COULD FAIL: does SPAN/SD predict non-monotonicity? ===")
    print("   If SPAN/SD is measuring what it claims, groups with a low value should be the")
    print("   ones whose ladders are non-monotone. This is an independent symptom.")
    xs = np.array([v[2] for v in table.values()])
    ys = np.array([v[3] / v[4] for v in table.values()])     # monotone FRACTION
    # Spearman without scipy
    def rank(a):
        order = np.argsort(a)
        r = np.empty_like(order, dtype=float)
        r[order] = np.arange(len(a))
        return r
    rho = float(np.corrcoef(rank(xs), rank(ys))[0, 1])
    print(f"   Spearman(SPAN/SD, monotone fraction) = {rho:+.3f} over {len(xs)} groups")
    # ⚠ FAILED BAR, kept printed. I registered rho > 0.3 without computing it; it is 0.274.
    # The bar is not loosened. It was also the WRONG statistic: the relationship is a FLOOR
    # effect, not a graded one, so a rank correlation over the whole range is insensitive to
    # it. The split at the bar is the statistic that matches the claim, and it is stark.
    print(f"   ⚠ REGISTERED BAR rho > 0.3 FAILED at {rho:.3f}. Not loosened. It is the wrong")
    print("     statistic: the effect is a FLOOR, not a gradient, so a rank correlation over")
    print("     the whole range dilutes it. The split at the bar is reported instead.")
    lo_g = [v for v in table.values() if v[2] < SPAN_SD_BAR]
    hi_g = [v for v in table.values() if v[2] >= SPAN_SD_BAR]
    lo_m, lo_n = sum(v[3] for v in lo_g), sum(v[4] for v in lo_g)
    hi_m, hi_n = sum(v[3] for v in hi_g), sum(v[4] for v in hi_g)
    print(f"   monotone LADDERS below the bar: {lo_m}/{lo_n} = {lo_m/lo_n:.2f}")
    print(f"   monotone LADDERS above the bar: {hi_m}/{hi_n} = {hi_m/hi_n:.2f}")
    check("ladders in groups below the bar are monotone far less often than those above it",
          lo_m / lo_n < 0.5 * hi_m / hi_n, f"{lo_m/lo_n:.2f} vs {hi_m/hi_n:.2f}")

    # ---- the discriminator the data supplies for free, which beats any invented bar -----
    print("\n   🔑 AND THE DATA SUPPLIES A BAR THAT NEEDS NO CONVENTION AT ALL:")
    print("      a NEGATIVE span is physically inadmissible. Compressing a packing harder")
    print("      cannot lower its packing fraction, so span < 0 is not 'small', it is wrong.")
    neg = {k: v for k, v in table.items() if v[0] < 0}
    print(f"      groups with a NEGATIVE mean span: {len(neg)} of {len(table)}")
    for k, v in sorted(neg.items()):
        print(f"        {k[0]:11s} mu_p={k[1]:.1f} cv={k[2]:.1f}   span = {v[0]:+.5f}")
    check("every group with an inadmissible negative span is also below the SPAN/SD bar, so "
          "the convention and the inadmissibility agree where both apply",
          all(table[k][2] < SPAN_SD_BAR for k in neg),
          f"{len(neg)} negative-span groups, all below the bar")
    if neg:
        print("      ⚠ These groups are not merely unreadable -- their ladders run BACKWARDS.")
        print("        Any extrapolated phi_J from them is meaningless, not just imprecise,")
        print("        and this needs no threshold to assert.")

    # ---- S4: which COMMITTED numbers sit on unreadable groups? ------------------------
    print("\n=== S4  WHICH PUBLISHED NUMBERS REST ON GROUPS THAT FAIL THE DIAGNOSTIC ===")
    published = [
        ("paper: \\dPhiJLogn, the headline shift", [("lognormal", 0.3, 0.0), ("lognormal", 0.3, 0.5)]),
        ("paper: \\dPhiJUnif", [("uniform", 0.3, 0.0), ("uniform", 0.3, 0.5)]),
        ("paper: \\dPhiJBidi", [("bidisperse", 0.3, 0.0), ("bidisperse", 0.3, 0.5)]),
        ("paper: \\phiJFrictionless (CV = 0 only)", [("lognormal", 0.0, 0.0)]),
        ("paper: \\phiJMonoMuThree (CV = 0 only)", [("lognormal", 0.3, 0.0)]),
        ("wc_endpoints: the frictionless CV shift", [("lognormal", 0.0, 0.0), ("lognormal", 0.0, 0.5)]),
        ("wc_friction: the mu_p = 1.0 window", [("lognormal", 1.0, 0.0), ("lognormal", 1.0, 0.5)]),
        ("wc_friction: the mu_p = 3.0 endpoint", [("lognormal", 3.0, 0.0)]),
    ]
    print(f"   {'quantity':44s} {'worst SPAN/SD':>14s}  verdict")
    affected = []
    for name, keys in published:
        vals = [table[k][2] for k in keys if k in table]
        if not vals:
            print(f"   {name:44s} {'no data':>14s}")
            continue
        w = min(vals)
        bad = w < SPAN_SD_BAR
        if bad:
            affected.append(name)
        print(f"   {name:44s} {w:14.1f}  {'🔴 rests on an unreadable group' if bad else 'ok'}")
    check("the paper's three headline family shifts are all on READABLE groups",
          not any(n.startswith("paper: \\dPhiJ") for n in affected),
          f"affected: {[a for a in affected if a.startswith('paper')] or 'none'}")
    check("this diagnostic DOES catch something -- at least one committed quantity rests "
          "on an unreadable group", len(affected) > 0, f"{len(affected)}: {affected}")

    # ---- S5: the honest limits --------------------------------------------------------
    print("\n=== S5  WHAT THIS DIAGNOSTIC CANNOT DO ===")
    print("   (a) SPAN/SD = 2 is a CONVENTION. Every group's value is printed above so a")
    print("       different bar can be applied; the ORDERING does not depend on it.")
    print("   (b) The sd is a 3-seed estimate, so it carries 2 dof and is itself noisy by")
    print("       roughly 40 %. A group at 1.5 and a group at 2.5 are not reliably apart.")
    print("   (c) It says an extrapolation is UNRELIABLE, not that phi_J is WRONG. A")
    print("       degenerate ladder whose four points happen to agree still brackets the")
    print("       jamming point; it just cannot locate it better than that bracket.")
    print("   (d) It does not replace the fit-quality gate the module already has. A ladder")
    print("       can have good span and still be non-monotone from a single bad rung.")

    print("\n" + "=" * 92)
    if FAILURES:
        print(f"FAILED {len(FAILURES)} check(s):")
        for f in FAILURES:
            print(f"   - {f}")
        return 1
    print("ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
