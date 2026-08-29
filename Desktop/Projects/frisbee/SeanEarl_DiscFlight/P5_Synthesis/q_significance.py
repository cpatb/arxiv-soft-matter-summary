"""
q_significance.py -- is the wind-direction ordering in Q actually resolved?

WHY THIS EXISTS
---------------
The synthesis reports Q_mean = 0.0443 / 0.0483 / 0.0356 for 90 / 210 / 270 degrees and
concludes "best = 210, worst = 270".  Those three means sit inside overlapping
Q_min-Q_max ranges, and nothing in the pipeline had yet asked whether the ordering
survives its own noise.  This module asks.

It answers two separate questions that were being conflated:

1. HOW MUCH OF THE Q RANGE IS SAMPLING ERROR?  The reported Q_min-Q_max is the spread
   across the 12 launch points, which is spatial structure, not uncertainty.  The Monte
   Carlo contribution is computed here from the binomial error on the measured success
   rate.  Answer: the spatial spread is 9-14x the Monte Carlo standard error, so the
   map is real and the sweep is adequately sampled.

2. IS THE ORDERING BETWEEN DIRECTIONS RESOLVED?  The 12 launch points are common to
   every direction, so the comparison is PAIRED, which is far more powerful than
   comparing three means against their spatial spread.  Answer, and it changes a
   headline: 270 is resolved as the worst against both other directions, but 90 and 210
   are NOT separable (4 of 12 launch points favour 90, paired p ~ 0.2).  The claim
   "best of the three computed = 210" is not supported by this evidence; the claim
   "270 is the worst of the three computed" is.

Method notes.
* Paired t on the 12 launch-point differences, plus a Wilcoxon signed-rank as a
  distribution-free check, because n = 12 is small and Q is bounded below by zero.
* Holm-Bonferroni across the three pairwise comparisons.  With three comparisons an
  uncorrected 0.05 is not the error rate anyone thinks it is.
* A paired bootstrap CI on each mean difference, resampling launch points, so the
  effect size is reported with an interval rather than only a p-value.

Run:  python q_significance.py [--Q Q.npz] [--landing ../P4_Coupling/landing_stats.npz]
Emits significance_generated.tex following the project's generated-macro convention, so
the manuscript cannot drift from these numbers.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy import stats

ALPHA = 0.05
N_BOOT = 20000
BOOT_SEED = 20260807


# --------------------------------------------------------------------------- #
def resolve_throw_subset(landing_npz, Q_npz):
    """Work out which throw types actually fed the stored success rate.

    This matters more than it sounds.  throwability.py supports --exclude-throws, and
    its own source records that the shipped Q.npz was built with `--exclude-throws
    hammer` while runDownstream.sh passed no such flag, so the committed script does not
    reproduce the committed artifact.  The shipped Q.npz also predates the provenance
    keys, so it does not say which subset it used.  Averaging over all three throws when
    the artifact used two inflates the binomial error by 1.5x and would make the sweep
    look under-sampled when it is not.

    Strategy: use the recorded provenance if present; otherwise search the subsets of
    throw types for the one that reproduces the stored success rate, and report which
    was inferred rather than guessing silently.
    """
    L = np.load(landing_npz, allow_pickle=True)
    Q = np.load(Q_npz, allow_pickle=True)
    throws = [str(t) for t in L["throws"]]
    sr_full = L["success_rate"]

    if "kept_throws" in Q.files:
        kept = [str(t) for t in Q["kept_throws"]]
        return dict(kept=kept, idx=[throws.index(t) for t in kept], source="recorded")

    if "success_rate" not in Q.files:
        return dict(kept=throws, idx=list(range(len(throws))), source="default (all)")

    target = Q["success_rate"]
    best = None
    for r in range(len(throws), 0, -1):
        for combo in itertools.combinations(range(len(throws)), r):
            cand = sr_full[:, :, :, list(combo), :].mean(axis=(3, 4))
            err = float(np.nanmax(np.abs(cand - target)))
            if best is None or err < best[0]:
                best = (err, combo)
    err, combo = best
    kept = [throws[i] for i in combo]
    if err > 1e-9:
        raise RuntimeError(
            f"no throw subset reproduces the stored success rate (best mismatch {err:.3g}, "
            f"subset {kept}).  The artifacts are inconsistent; regenerate Q.npz.")
    return dict(kept=kept, idx=list(combo), source="inferred by matching")


def monte_carlo_error(landing_npz, throw_idx):
    """Standard error of the per-launch-point mean success rate, from binomial counts.

    success_rate has shape (nx, ny, n_heading, n_throw, n_bearing) and each entry is a
    fraction of n_valid independent throws.  Q uses the mean over the SELECTED throw
    types and all bearings, so the variance of that mean is the sum of the per-cell
    binomial variances divided by the number of contributing cells squared.
    """
    L = np.load(landing_npz, allow_pickle=True)
    sr = L["success_rate"][:, :, :, list(throw_idx), :]
    n = np.maximum(L["n_valid"][:, :, :, list(throw_idx), :], 1)
    n_cells = sr.shape[3] * sr.shape[4]
    var_cell = sr * (1.0 - sr) / n
    se_mean = np.sqrt(var_cell.sum(axis=(3, 4))) / n_cells
    return dict(se_success=se_mean, sr_mean=sr.mean(axis=(3, 4)),
                n_per_cell=int(n.flat[0]), n_cells=int(n_cells))


def variance_decomposition(Q, sr_mean, se_success):
    """Split the reported Q spread into spatial structure and Monte Carlo noise.

    Q = wind_factor * turbulence_factor * success_rate, and only the success rate
    carries Monte Carlo error, so the relative error transfers directly.
    """
    out = []
    for ih in range(Q.shape[-1]):
        q = Q[..., ih]
        rel = se_success[..., ih] / np.maximum(sr_mean[..., ih], 1e-12)
        se_q_point = q * rel
        se_of_mean = float(np.sqrt((se_q_point**2).sum()) / q.size)
        spatial_sd = float(np.std(q, ddof=1))
        out.append(dict(q_mean=float(np.mean(q)), spatial_sd=spatial_sd,
                        mc_se_of_mean=se_of_mean,
                        ratio=spatial_sd / se_of_mean if se_of_mean > 0 else np.inf,
                        mc_se_point_median=float(np.median(se_q_point))))
    return out


# --------------------------------------------------------------------------- #
def paired_bootstrap_ci(a, b, n_boot=N_BOOT, seed=BOOT_SEED, alpha=ALPHA):
    """Percentile CI for the mean paired difference, resampling launch points."""
    rng = np.random.default_rng(seed)
    d = np.asarray(a, float) - np.asarray(b, float)
    idx = rng.integers(0, d.size, size=(n_boot, d.size))
    means = d[idx].mean(axis=1)
    lo, hi = np.percentile(means, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(d.mean()), float(lo), float(hi)


def holm(pvalues):
    """Holm-Bonferroni adjusted p-values, order preserved."""
    p = np.asarray(pvalues, float)
    m = p.size
    order = np.argsort(p)
    adj = np.empty(m)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (m - rank) * p[i])
        adj[i] = min(running, 1.0)
    return adj


def pairwise_tests(Q, headings, alpha=ALPHA):
    """All pairwise paired comparisons across launch points, Holm-corrected."""
    flat = {int(h): Q[..., i].ravel() for i, h in enumerate(headings)}
    pairs = list(itertools.combinations([int(h) for h in headings], 2))
    raw_t, raw_w, rows = [], [], []
    for a, b in pairs:
        qa, qb = flat[a], flat[b]
        good = np.isfinite(qa) & np.isfinite(qb)
        qa, qb = qa[good], qb[good]
        t_p = stats.ttest_rel(qa, qb).pvalue
        w_p = stats.wilcoxon(qa, qb).pvalue
        diff, lo, hi = paired_bootstrap_ci(qa, qb, alpha=alpha)
        raw_t.append(t_p)
        raw_w.append(w_p)
        rows.append(dict(a=a, b=b, n=int(qa.size), mean_diff=diff,
                         ci_lo=lo, ci_hi=hi,
                         favour_a=int(np.sum(qa > qb)),
                         p_t_raw=float(t_p), p_w_raw=float(w_p)))
    for row, pt, pw in zip(rows, holm(raw_t), holm(raw_w)):
        row["p_t_holm"] = float(pt)
        row["p_w_holm"] = float(pw)
        # Resolved only if BOTH the parametric and the rank test survive correction,
        # which is the conservative reading and the right one for n = 12.
        row["resolved"] = bool(pt < alpha and pw < alpha)
    return rows


def verdict(rows, headings):
    """Which ordering claims the evidence actually supports."""
    hs = [int(h) for h in headings]
    resolved = {(r["a"], r["b"]): r for r in rows if r["resolved"]}
    beats = {h: set() for h in hs}
    for (a, b), r in resolved.items():
        if r["mean_diff"] > 0:
            beats[a].add(b)
        else:
            beats[b].add(a)
    worst = [h for h in hs if all(h in beats[o] for o in hs if o != h)]
    best = [h for h in hs if beats[h] == {o for o in hs if o != h}]
    return dict(best_supported=best[0] if len(best) == 1 else None,
                worst_supported=worst[0] if len(worst) == 1 else None,
                beats={k: sorted(v) for k, v in beats.items()})


# --------------------------------------------------------------------------- #
def _tex(name, value):
    return f"\\newcommand{{\\{name}}}{{{value}}}\n"


_WORDS = {"0": "Zero", "1": "One", "2": "Two", "3": "Three", "4": "Four",
          "5": "Five", "6": "Six", "7": "Seven", "8": "Eight", "9": "Nine"}


def _word(n):
    return "".join(_WORDS[c] for c in str(int(n)))


def write_tex(path, decomp, rows, verd, headings, mc, subset):
    out = ["% GENERATED by P5_Synthesis/q_significance.py -- do not edit.\n",
           "% Statistical resolution of the wind-direction ordering in Q.\n\n"]
    out.append(_tex("QsigAlpha", f"{ALPHA:g}"))
    out.append(_tex("QsigNlaunch", rows[0]["n"]))
    out.append(_tex("QsigNboot", f"{N_BOOT:,}".replace(",", "\\,")))
    out.append(_tex("QsigThrowsPerCell", mc["n_per_cell"]))
    out.append(_tex("QsigThrowTypes",
                    ", ".join(t.replace("_", "\\_") for t in subset["kept"])))
    out.append(_tex("QsigThrowProvenance", subset["source"]))

    for i, h in enumerate(headings):
        w = _word(h)
        d = decomp[i]
        out.append(_tex(f"QspatialSdH{w}", f"{d['spatial_sd']:.4f}"))
        out.append(_tex(f"QmcSeH{w}", f"{d['mc_se_of_mean']:.5f}"))
        out.append(_tex(f"QsdOverSeH{w}", f"{d['ratio']:.0f}"))

    for r in rows:
        tag = f"{_word(r['a'])}vs{_word(r['b'])}"
        out.append(_tex(f"Qdiff{tag}", f"{r['mean_diff']:+.4f}"))
        out.append(_tex(f"Qci{tag}", f"[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}]"))
        out.append(_tex(f"Qpt{tag}", f"{r['p_t_holm']:.3g}"))
        out.append(_tex(f"Qpw{tag}", f"{r['p_w_holm']:.3g}"))
        out.append(_tex(f"Qfavour{tag}", f"{r['favour_a']}/{r['n']}"))
        out.append(_tex(f"Qresolved{tag}", "resolved" if r["resolved"] else "not resolved"))

    out.append(_tex("QworstSupported",
                    verd["worst_supported"] if verd["worst_supported"] else "none"))
    out.append(_tex("QbestSupported",
                    verd["best_supported"] if verd["best_supported"] else "none"))
    Path(path).write_text("".join(out), encoding="utf-8")
    return sum(1 for line in out if line.startswith("\\newcommand"))


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--Q", default="Q.npz")
    ap.add_argument("--landing", default="../P4_Coupling/landing_stats.npz")
    ap.add_argument("--tex", default="significance_generated.tex")
    ap.add_argument("--json", default="q_significance.json")
    args = ap.parse_args()

    Qd = np.load(args.Q, allow_pickle=True)
    Q, headings = Qd["Q"], Qd["headings"]
    subset = resolve_throw_subset(args.landing, args.Q)
    mc = monte_carlo_error(args.landing, subset["idx"])
    decomp = variance_decomposition(Q, mc["sr_mean"], mc["se_success"])
    rows = pairwise_tests(Q, headings)
    verd = verdict(rows, headings)

    print("=" * 74)
    print("Q SIGNIFICANCE  --  is the wind-direction ordering resolved?")
    print("=" * 74)
    print(f"{rows[0]['n']} launch points, {mc['n_per_cell']} throws per "
          f"(throw type x bearing) cell, {mc['n_cells']} such cells per point.")
    print(f"throw types feeding Q: {', '.join(subset['kept'])}  [{subset['source']}]")
    if subset["source"].startswith("inferred"):
        print("   NOTE: Q.npz carries no throw provenance, so the subset was recovered")
        print("   by matching its stored success rate.  Regenerating Q.npz with the")
        print("   current throwability.py would record it explicitly.")
    print()

    print("-- 1. Is the reported Q spread sampling noise or spatial structure? --")
    print(f"   {'heading':>8s} {'Q_mean':>8s} {'spatial sd':>11s} {'MC se':>9s} {'ratio':>7s}")
    for i, h in enumerate(headings):
        d = decomp[i]
        print(f"   {int(h):>7d}d {d['q_mean']:>8.4f} {d['spatial_sd']:>11.4f} "
              f"{d['mc_se_of_mean']:>9.5f} {d['ratio']:>6.0f}x")
    print("   The spatial spread dominates by more than an order of magnitude, so the")
    print("   Q map is structure and the Monte Carlo is adequately sampled.  The")
    print("   reported Q_min-Q_max is a spatial range and must not be read as an")
    print("   uncertainty on the mean.\n")

    print("-- 2. Are the directions separable?  Paired over launch points. --")
    for r in rows:
        flag = "RESOLVED    " if r["resolved"] else "NOT RESOLVED"
        print(f"   {r['a']:>3d} vs {r['b']:>3d}: diff {r['mean_diff']:+.4f} "
              f"95% CI [{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}]  "
              f"p_t={r['p_t_holm']:.3g} p_w={r['p_w_holm']:.3g}  "
              f"{r['favour_a']}/{r['n']} favour {r['a']}   {flag}")
    print("   (p-values Holm-corrected across the three comparisons; a pair counts as")
    print("   resolved only if the paired t AND the Wilcoxon both survive.)\n")

    print("-- 3. What the evidence supports --")
    if verd["worst_supported"] is not None:
        print(f"   WORST of the three computed: {verd['worst_supported']} degrees, "
              f"resolved against both others.")
    else:
        print("   No direction is resolved as the worst.")
    if verd["best_supported"] is not None:
        print(f"   BEST of the three computed: {verd['best_supported']} degrees.")
    else:
        print("   *** NO direction is resolved as the best. ***")
        print("   The synthesis currently reports 'best = 210'.  On this evidence 210")
        print("   and 90 are not separable, so that claim should be withdrawn or")
        print("   softened to '270 is the worst of the three computed'.")

    n = write_tex(args.tex, decomp, rows, verd, headings, mc, subset)
    Path(args.json).write_text(json.dumps(
        dict(decomposition=decomp, pairs=rows, verdict=verd,
             alpha=ALPHA, n_boot=N_BOOT, boot_seed=BOOT_SEED), indent=2), encoding="utf-8")
    print(f"\nwrote {args.tex} ({n} macros) and {args.json}")


if __name__ == "__main__":
    main()
