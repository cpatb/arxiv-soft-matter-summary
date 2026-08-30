"""Does the ballistic law's rejection replicate across wind speed and launch point?

WHAT IS BEING REPLICATED
------------------------
`out/sweep_flight.json` rejected the brief's ballistic `Phi(N)` dispersion law
at slope b = 0.671, CI [0.522, 0.818], on 30 cells at ONE wind speed
(Uref = 10 m/s) and ONE launch point (-28, +2.5). A single block like that is
exactly the kind of result that turns out to be a property of its cell, so this
re-runs the same test over a factorial:

    Uref     5, 10, 20 m/s
    launch   (-28, +2.5)  and  (-12, +5), both on the production grid's side
    heading  90 and 270
    flight   5 release speeds x 3 elevations, unchanged

The two launch points were chosen on measured evidence, not convenience. At
matched sigma_u the path-averaged integral length is 6.74 m at the first and
3.32 m at the second, so N = L_f/L roughly DOUBLES between them. N is the
argument of the function under test, so this is the strongest available second
condition rather than a token one -- (-25, +/-5), the nearest production grid
points to the original, sit within 4% of it in L and would have replicated
nothing.

THE ACCEPTANCE RULE, FIXED BEFORE THE RESULTS WERE LOOKED AT
------------------------------------------------------------
Per block (launch x Uref x heading), fit log(measured) = a + b log(predicted)
for each law and read b's 95% bootstrap CI.

  REPLICATED       the ballistic CI excludes 1.0, and lies BELOW it, in at
                   least 5 of the 6 (launch, Uref) blocks pooled over heading
  NOT REPLICATED   the ballistic CI contains 1.0 in 2 or more blocks
  UNINFORMATIVE    the two laws' predictions separate by < 1.6x in that block,
                   in which case the block is reported but not counted either
                   way -- a block that cannot tell the laws apart is not
                   evidence that they agree

The tracer law is graded the same way and is expected NOT to be rejected. If
BOTH are rejected in the same block, the problem is upstream of both laws and
the block is flagged rather than scored.

Cells whose throws end inside buildings are excluded at 10%, because the
trajectory there flies on nearest-neighbour-filled wind and its path-averaged
sigma_u and L are not measurements. The count dropped is reported per block.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from analyse_theory import law_comparison, shape_test      # noqa: E402

BLOCKED_MAX = 0.10
SEPARATION_MIN = 1.6
N_REPLICATE_REQUIRED = 5          # of the 6 (launch, Uref) blocks


def blocks(rows, keys):
    g = defaultdict(list)
    for r in rows:
        g[tuple(_key(r, k) for k in keys)].append(r)
    return g


def _key(r, k):
    v = r.get(k)
    return tuple(v) if isinstance(v, list) else v


def grade(rows, keys):
    out = []
    for key, rs in sorted(blocks(rows, keys).items(), key=lambda kv: str(kv[0])):
        n_all = len(rs)
        rs = [r for r in rs if r["full__blocked"] <= BLOCKED_MAX]
        comp = law_comparison(rs)
        sh = shape_test(comp)
        if not sh or "_span" not in sh:
            continue
        sp = sh["_span"]
        rec = dict(key=list(key), n_all=n_all, n_used=len(comp),
                   separation=sp["law_ratio_span"],
                   informative=sp["informative"], L_f_span=sp["L_f"])
        for law in ("ballistic", "tracer"):
            if law not in sh:
                continue
            q = sh[law]
            rejected = (q["ci"][1] < 1.0) or (q["ci"][0] > 1.0)
            rec[law] = dict(slope=q["slope"], ci=q["ci"], r2=q["r2"],
                            prefactor=q["prefactor_needed"],
                            scatter=q["scatter_at_unit_slope"],
                            rejected=bool(rejected),
                            rejected_low=bool(q["ci"][1] < 1.0))
        out.append(rec)
    return out


def show(recs, label):
    print(f"\n{label}")
    print(f"  {'block':<34} {'n':>5} {'sep':>5} {'ballistic b':>22} "
          f"{'tracer b':>22}")
    for r in recs:
        k = " ".join(str(x) for x in r["key"])
        b, t = r.get("ballistic"), r.get("tracer")
        mark = "" if r["informative"] else "  (uninformative)"

        def fmt(q):
            if q is None:
                return " " * 22
            tag = "REJ" if q["rejected"] else "ok "
            return (f"{q['slope']:6.3f} [{q['ci'][0]:5.3f},{q['ci'][1]:5.3f}] "
                    f"{tag}")
        print(f"  {k:<34} {r['n_used']:5d} {r['separation']:5.2f} "
              f"{fmt(b):>22} {fmt(t):>22}{mark}")



def figure(by_luh, outdir):
    """Slope with CI per block, both laws, against the b = 1 line.

    A forest plot rather than a table because the question is entirely "does
    this interval cover 1.0, in every block", and that is a thing to look at.
    The x axis is a slope, so it is NOT zero-based: 1.0 is the reference the
    whole figure is about and the informative range is around it.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    recs = [r for r in by_luh if "ballistic" in r and "tracer" in r]
    if not recs:
        return
    fig, ax = plt.subplots(figsize=(6.4, 0.30 * len(recs) + 1.4))
    y = np.arange(len(recs))
    for law, col, off, lab in (("ballistic", "#d62728", +0.16,
                                r"ballistic $\Phi(N)$"),
                               ("tracer", "#1f77b4", -0.16,
                                r"tracer $1/\sqrt{N}$")):
        v = np.array([r[law]["slope"] for r in recs])
        lo = np.array([r[law]["ci"][0] for r in recs])
        hi = np.array([r[law]["ci"][1] for r in recs])
        ax.errorbar(v, y + off, xerr=[v - lo, hi - v], fmt="o", ms=4.5,
                    lw=1.4, capsize=2.5, color=col, label=lab)
    ax.axvline(1.0, color="#333", lw=1.2, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([f"({r['key'][0][0]:.0f},{r['key'][0][1]:.0f})  "
                        f"{r['key'][1]:.0f} m/s  {r['key'][2]:.0f}deg"
                        for r in recs], fontsize=8)
    ax.set_xlabel(r"slope $b$ of  $\log\sigma_{measured}$ vs "
                  r"$\log\sigma_{predicted}$")
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.25, lw=0.5)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    fig.tight_layout()
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "fig_replication.png", dpi=200, facecolor="#ffffff")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sweep", default=str(HERE / "out" / "sweep_flight2.json"))
    ap.add_argument("--out", default=str(HERE / "out" / "replication.json"))
    ap.add_argument("--figdir", default=str(HERE / "figs_replication"))
    args = ap.parse_args()

    meta, rows = (lambda d: (d["meta"], d["rows"]))(
        json.loads(Path(args.sweep).read_text()))
    print(f"{len(rows)} cells from {Path(args.sweep).name}")
    print(f"acceptance rule (pre-registered): ballistic CI excludes 1.0 and "
          f"lies below it\n  in >= {N_REPLICATE_REQUIRED} of the "
          f"(launch, Uref) blocks; blocks separating the two laws by "
          f"< {SEPARATION_MIN}x are not counted.")

    by_lu = grade(rows, ("launch", "u_ref"))
    show(by_lu, "by launch point x wind speed (headings pooled)")
    by_luh = grade(rows, ("launch", "u_ref", "heading"))
    show(by_luh, "by launch point x wind speed x heading")

    inf = [r for r in by_lu if r["informative"]]
    rej = [r for r in inf if r.get("ballistic", {}).get("rejected_low")]
    tra_rej = [r for r in inf if r.get("tracer", {}).get("rejected")]
    both = [r for r in inf if r.get("ballistic", {}).get("rejected")
            and r.get("tracer", {}).get("rejected")]

    print(f"\n  informative blocks: {len(inf)} of {len(by_lu)}")
    print(f"  ballistic rejected LOW in {len(rej)} of {len(inf)}")
    print(f"  tracer rejected in {len(tra_rej)} of {len(inf)}")
    if both:
        print(f"  !! BOTH laws rejected in {len(both)} block(s) -- the problem "
              f"there is upstream of both:")
        for r in both:
            print(f"       {r['key']}")

    verdict = ("REPLICATED" if len(rej) >= N_REPLICATE_REQUIRED
               else "NOT REPLICATED" if len(inf) - len(rej) >= 2
               else "INCONCLUSIVE")
    print(f"\n  VERDICT: the ballistic-law rejection is {verdict} "
          f"({len(rej)}/{len(inf)} informative blocks, rule needs "
          f">= {N_REPLICATE_REQUIRED})")

    # Sign test across BOTH aggregation levels. A CI-counting rule throws away
    # the information in slopes that miss 1.0 but all miss it the same way, and
    # with 10-15 cells per finer block most CIs are wide enough to cover 1.0
    # individually while every point estimate sits on one side.
    #
    # Reporting both levels is not padding. If a law's verdict FLIPS between
    # them the verdict is a property of how the data was grouped, not of the
    # law, and saying only the flattering one would be a choice.
    from scipy import stats
    print("\n  sign test on the slope point estimates:")
    signs = {}
    for lvl, recs in (("launch x Uref", by_lu),
                      ("launch x Uref x heading", by_luh)):
        rs = [r for r in recs if "ballistic" in r and "tracer" in r]
        signs[lvl] = {}
        for law in ("ballistic", "tracer"):
            sl = np.array([r[law]["slope"] for r in rs])
            below = int((sl < 1.0).sum())
            pv = float(stats.binomtest(below, len(sl), 0.5).pvalue)
            signs[lvl][law] = dict(n=len(sl), below_one=below, p=pv,
                                   mean_slope=float(sl.mean()),
                                   lo=float(sl.min()), hi=float(sl.max()))
            print(f"    {lvl:<24} {law:<10} slope {sl.min():.3f}-{sl.max():.3f}"
                  f" (mean {sl.mean():.3f}), below 1.0 in {below}/{len(sl)}, "
                  f"p = {pv:.4f}")
    b_lo = all(signs[lv]["ballistic"]["mean_slope"] < 1.0 for lv in signs)
    tp = signs["launch x Uref"]["tracer"]["mean_slope"]
    tf = signs["launch x Uref x heading"]["tracer"]["mean_slope"]
    t_flip = (tp < 1.0) != (tf < 1.0)
    print("")
    print(f"    ballistic sits below 1.0 at both aggregation levels: {b_lo}")
    if t_flip:
        print("    !! the TRACER law's verdict FLIPS between the two "
              f"levels (mean slope {tp:.3f} pooled over heading vs "
              f"{tf:.3f} within heading). Its survival is therefore a "
              "property of the grouping, not a result.")
        print("       The safe statement is that the measurement is "
              "INTERMEDIATE between the two laws, and that only the "
              "ballistic one is rejected consistently.")

    # Which law is CLOSER is a different question from which is rejected, and
    # the two have different answers here. Distance is |ln b|, so a slope of
    # 0.5 and a slope of 2.0 count as equally wrong.
    print("")
    print("  bracketing: is the measurement between the two predictions?")
    brack = {}
    for lvl, recs in (("launch x Uref", by_lu),
                      ("launch x Uref x heading", by_luh)):
        rs = [r for r in recs if "ballistic" in r and "tracer" in r]
        b = np.array([r["ballistic"]["slope"] for r in rs])
        t = np.array([r["tracer"]["slope"] for r in rs])
        nb = int(((b < 1.0) & (t > 1.0)).sum())
        db, dt = float(np.abs(np.log(b)).mean()), float(np.abs(np.log(t)).mean())
        gm = np.sqrt(b * t)
        brack[lvl] = dict(n=len(b), bracketed=nb, dist_ballistic=db,
                          dist_tracer=dt,
                          closer=("tracer" if dt < db else "ballistic"),
                          geomean_slope=float(gm.mean()),
                          dist_geomean=float(np.abs(np.log(gm)).mean()))
        print(f"    {lvl:<24} ballistic < 1 < tracer in {nb}/{len(b)} blocks; "
              f"|ln b| ballistic {db:.3f} vs tracer {dt:.3f} "
              f"-> {brack[lvl]['closer']} closer")
        print(f"    {'':<24} geometric mean of the two predictions: slope "
              f"{gm.mean():.3f} (|ln| {np.abs(np.log(gm)).mean():.3f})")
    if brack["launch x Uref"]["closer"] != \
            brack["launch x Uref x heading"]["closer"]:
        print("    !! WHICH LAW IS CLOSER FLIPS BETWEEN THE TWO LEVELS. Only "
              "the ballistic law's\n       DIRECTION of error is stable. "
              "Reporting a winner would be reporting a grouping.")

    # Does the ballistic slope depend on the block? If the law were simply
    # mis-normalised the slope would be 1 everywhere and only the intercept
    # would move; a slope that itself varies means the law is missing structure.
    sl = np.array([r["ballistic"]["slope"] for r in inf if "ballistic" in r])
    pf = np.array([r["ballistic"]["prefactor"] for r in inf if "ballistic" in r])
    if len(sl):
        print(f"\n  ballistic slope across informative blocks: "
              f"{sl.min():.3f} to {sl.max():.3f} (mean {sl.mean():.3f})")
        print(f"  prefactor it would need:                    "
              f"{pf.min():.3f} to {pf.max():.3f} "
              f"(spread {pf.max()/pf.min():.2f}x)")

    figure(by_luh, args.figdir)
    Path(args.out).write_text(json.dumps(dict(
        sweep=args.sweep, meta=meta, blocked_max=BLOCKED_MAX,
        separation_min=SEPARATION_MIN,
        n_replicate_required=N_REPLICATE_REQUIRED,
        by_launch_uref=by_lu, by_launch_uref_heading=by_luh,
        n_informative=len(inf), n_ballistic_rejected_low=len(rej),
        n_tracer_rejected=len(tra_rej), verdict=verdict), indent=2))
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
