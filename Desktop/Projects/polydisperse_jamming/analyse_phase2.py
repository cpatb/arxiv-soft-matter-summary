"""analyse_phase2.py — turn the phi_J sweep into a verdict on mechanism M2.

Reads every out/phase2*.jsonl shard, fits phi_J per (cv, mu, seed) by extrapolating
phi(p_hat) -> 0, aggregates over seeds, and grades the PRE-REGISTERED M2 test from
SPEC_2026-08-31.md §5:

    M2 is SUFFICIENT if the measured shift in phi_J between the author's low-CV and
    high-CV samples changes (phi_J - phi)^-2 by more than ONE DECADE.

That test cannot be evaluated without naming a phi, because the amplification
(phi_J - phi)^-2 depends violently on how close the suspension sits to jamming. So the
script reports the amplification as a FUNCTION of phi, and additionally reports the
`phi_needed` at which the measured Delta phi_J would reach one decade. Comparing
`phi_needed` against the packing fractions the lab actually works at is the verdict.

Outputs
-------
  out/phase2_fits.csv        one row per (cv, mu, seed): phi_J and the exponent check
  out/phase2_summary.csv     one row per (cv, mu): phi_J mean +/- seed scatter
  out/phase2_verdict.txt     the graded M2 test, printed and saved
"""

from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

from phase2_phiJ import fit_phi_J

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "out")

# The reduced-pressure ladder is only in the critical regime at its low end. Points
# above this are dropped from the phi_J extrapolation and the drop is REPORTED.
P_HAT_FIT_MAX = 1.1e-3

# The exponent check needs a ladder that a monotonic power law could describe. TWO
# conditions, and the second is the one that actually bites:
#
#   range  > 2 x seed scatter   -- the ladder must move more than replicate noise
#   MONOTONE in p              -- phi = phi_J + c*p^b is monotonic by construction,
#                                 so fitting it to a non-monotonic series returns a
#                                 number with no meaning
#
# The frictionless ladders are the case in point: phi = 0.6434, 0.6319, 0.6352,
# 0.6313 over 1.5 decades. Range 0.012 against a seed scatter of 0.00111, so they
# clear the range gate ten times over -- but they go DOWN, UP, DOWN, and one of them
# returned b = 0.248 against a predicted 1. A range gate alone would have passed it.
# The frictional ladders are cleanly monotonic (0.5949, 0.5827, 0.5748).
EXPONENT_MIN_RANGE_IN_SD = 2.0

# Packing fractions the lab actually works at, for the amplification table.
PHI_OF_INTEREST = (0.40, 0.44, 0.48, 0.50, 0.52, 0.55, 0.58, 0.60)


# --- system size ---------------------------------------------------------------------
PRIMARY_N = 1500
"""The system size every claim in this project is quoted at.

`out/` also holds N = 3000 cells. Before 2026-09-08 nothing downstream distinguished them
and several analyses were silently mixing two sizes; see `N_GROUPING_2026-09-08.md`.
`analyse_phase2.py` now KEYS on N and grades per size. The consumers that cannot -- because
the N = 3000 coverage does not span the CV range their comparison needs -- restrict to this
size explicitly through `restrict_to_primary_N` rather than by accident.
"""


def restrict_to_primary_N(recs, why: str, n: int = PRIMARY_N):
    """Keep only cells at system size ``n``, and PRINT what that dropped.

    Never drop silently: the count and the sizes present are printed every run, so a
    later N sweep cannot quietly stop being included.
    """
    import collections as _c
    sizes = _c.Counter(r.get("n_requested") for r in recs)
    kept = [r for r in recs if r.get("n_requested") == n]
    dropped = len(recs) - len(kept)
    print(f"N RESTRICTION ({why}): keeping N = {n}; "
          f"{len(kept)} of {len(recs)} cells, {dropped} dropped. "
          f"sizes present: {dict(sorted(sizes.items(), key=lambda kv: (kv[0] is None, kv[0])))}")
    return kept


def load_all() -> list[dict]:
    """Load every shard, DEDUPLICATED by cell key.

    Shard files from an earlier, killed launch are still on disk (a 9-shard run was
    replaced by an 8-shard one), and the same cell can therefore appear twice. The
    cell key encodes cv, mu, seed, p_hat and N, and the protocol is deterministic
    given the seed, so duplicates are genuinely redundant -- but counting them twice
    would fake seed replication that does not exist.
    """
    by_key: dict[str, dict] = {}
    dupes = 0
    files = sorted(glob.glob(os.path.join(OUT, "phase2*.jsonl")))
    n_lines = 0
    for f in files:
        with open(f, encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                n_lines += 1
                k = r.get("key")
                if k is None:
                    continue
                if k in by_key:
                    dupes += 1
                    continue
                by_key[k] = r
    print(f"loaded {n_lines} lines from {len(files)} shard files -> "
          f"{len(by_key)} unique cells ({dupes} duplicates dropped)")
    return list(by_key.values())


def main() -> int:
    recs = load_all()
    if not recs:
        print("no records yet")
        return 1

    n_raw = len(recs)
    failed = [r for r in recs if r.get("failed")]
    recs = [r for r in recs if not r.get("failed")]
    print(f"N: {n_raw} raw -> {len(recs)} after dropping {len(failed)} failed cells")

    nonconv = [r for r in recs if not r.get("converged")]
    kept = [r for r in recs if r.get("converged")]
    print(f"N: -> {len(kept)} after dropping {len(nonconv)} non-converged cells")

    # The dropped cells are NOT a random sample -- non-convergence rises from 5% at
    # CV = 0 to 46% at CV = 0.5, and CV is the axis this whole project measures. That
    # is exactly the situation where a selection can manufacture a trend, so the bias
    # is MEASURED here rather than assumed away: for every (cv, mu, p_hat) cell,
    # compare mean phi with and without the non-converged members.
    bias = []
    cellmap: dict[tuple, list[dict]] = {}
    for r in recs:
        cellmap.setdefault((r.get("family", "lognormal"), r["cv"], r["mu_p"],
                            r["p_hat_target"], r["n_requested"]), []).append(r)
    for k, g in cellmap.items():
        c = [x for x in g if x.get("converged")]
        if c and len(c) < len(g):
            bias.append(float(np.mean([x["phi"] for x in c])
                              - np.mean([x["phi"] for x in g])))
    if bias:
        mx = max(abs(b) for b in bias)
        print(f"SELECTION CHECK: dropping the non-converged cells shifts phi by at "
              f"most {mx:.4f} across {len(bias)} affected cells "
              f"(median |shift| {np.median([abs(b) for b in bias]):.5f}). "
              f"{'NEGLIGIBLE' if mx < 0.005 else '*** NOT NEGLIGIBLE ***'}")
    recs = kept

    finite_size = [r for r in recs if r.get("d_max_over_L", 0) > 0.25]
    if finite_size:
        print(f"WARNING: {len(finite_size)} cells carry a finite-size flag "
              f"(d_max/L > 0.25); they are KEPT but listed in the csv")

    # ---- group and fit -------------------------------------------------------
    # N IS PART OF THE LADDER KEY. It was not until 2026-09-08, and out/ holds cells
    # at two system sizes, so 9 of 95 ladders were fitting N = 1500 and N = 3000 points
    # together at the same p_hat -- a system-size difference entering as ladder scatter.
    # diag_n_pooling.py measures what that cost: phi_J moved up to 1.30 seed scatters,
    # and the pre-registered exponent anchor FAILS either way (median b 0.567 pooled,
    # 0.488 split), so the pooling never explained the anchor's failure.
    groups: dict[tuple, list[dict]] = {}
    for r in recs:
        fam = r.get("family", "lognormal")
        groups.setdefault((fam, r["cv"], r["mu_p"], r["seed"], r["n_requested"]),
                          []).append(r)

    rows = []
    for (fam, cv, mu, seed, n_req), g in sorted(groups.items()):
        g_fit = [r for r in g if r["p_hat"] <= P_HAT_FIT_MAX]
        dropped = len(g) - len(g_fit)
        if len(g_fit) < 2:
            continue
        p = np.array([r["p_hat"] for r in g_fit])
        f = np.array([r["phi"] for r in g_fit])
        fit = fit_phi_J(p, f)
        cvr = float(np.mean([r["cv_realised"] for r in g_fit]))
        d32 = float(np.mean([r["d32_realised"] for r in g_fit]))
        rows.append({"family": fam, "cv_nom": cv, "cv_real": cvr,
                     "d32": d32, "mu_p": mu, "n_requested": n_req,
                     "seed": seed, "n_ladder": len(g_fit), "n_dropped_highp": dropped,
                     "phi_range": float(f.max() - f.min()),
                     "phi_monotone": bool(np.all(np.diff(f[np.argsort(p)]) >= 0)
                                          or np.all(np.diff(f[np.argsort(p)]) <= 0)),
                     **fit})

    if not rows:
        print("no complete ladders yet -- rerun when more cells have finished")
        return 1

    import csv
    fp = os.path.join(OUT, "phase2_fits.csv")
    with open(fp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {fp}  ({len(rows)} ladders)")

    # ---- the exponent anchor, gated on the ladder having any signal -----------
    lines = []
    lines.append("=" * 74)
    lines.append("PHASE 2 -- phi_J(CV) and the M2 verdict")
    lines.append("=" * 74)
    lines.append("")
    lines.append("ANCHOR: a harmonic contact law predicts p ~ (phi - phi_J)^1, so the")
    lines.append("free-exponent fit phi = phi_J + c*p_hat^b must return b ~ 1.")
    lines.append("")
    lines.append("GATE: a ladder counts only if its phi range exceeds "
                 f"{EXPONENT_MIN_RANGE_IN_SD}x the seed scatter AND it is")
    lines.append("MONOTONE in p. phi_J + c*p^b is monotonic by construction, so a")
    lines.append("non-monotonic ladder cannot be described by it at any b. The")
    lines.append("frictionless ladders go down/up/down and one returned b = 0.248.")
    # seed scatter of phi at fixed (family, cv, mu, p_hat)
    scat: dict[tuple, list[float]] = {}
    for r in recs:
        scat.setdefault((r.get("family", "lognormal"), r["cv"], r["mu_p"],
                         r["p_hat_target"], r["n_requested"]), []).append(r["phi"])
    sds = [float(np.std(v, ddof=1)) for v in scat.values() if len(v) > 1]
    seed_sd = float(np.median(sds)) if sds else float("nan")
    lines.append(f"  median seed scatter in phi = {seed_sd:.5f} "
                 f"(from {len(sds)} replicated cells)")
    if np.isfinite(seed_sd) and seed_sd > 0:
        def informative(r):
            return (r.get("phi_range", 0) > EXPONENT_MIN_RANGE_IN_SD * seed_sd
                    and r.get("phi_monotone", False))
        cand = [r for r in rows if np.isfinite(r.get("b_free", np.nan))]
        usable = [r for r in cand if informative(r)]
        excluded = [r for r in cand if not informative(r)]
    else:
        usable = [r for r in rows if np.isfinite(r.get("b_free", np.nan))]
        excluded = []
    n_flat = sum(1 for r in excluded
                 if r.get("phi_range", 0) <= EXPONENT_MIN_RANGE_IN_SD * seed_sd)
    n_nonmono = sum(1 for r in excluded if not r.get("phi_monotone", False))
    lines.append(f"  ladders with a free-exponent fit: {len(usable)+len(excluded)}; "
                 f"{len(excluded)} excluded ({n_flat} flat, {n_nonmono} "
                 f"non-monotonic)")
    b = np.array([r["b_free"] for r in usable], float)
    if b.size:
        lines.append(f"  b_free over {b.size} INFORMATIVE ladders: "
                     f"median {np.median(b):.3f}  mean {b.mean():.3f}  IQR "
                     f"[{np.percentile(b,25):.3f}, {np.percentile(b,75):.3f}]")
        verdict = "PASS" if 0.6 <= np.median(b) <= 1.6 else "FAIL"
        lines.append(f"  -> {verdict} (bar: median b in [0.6, 1.6])")
    else:
        lines.append("  -> UNINFORMATIVE: no ladder both had >= 4 points and moved")
        lines.append("     more than the seed scatter. The exponent is not measured,")
        lines.append("     and this is reported as such rather than graded.")
    if excluded:
        eb = np.array([r["b_free"] for r in excluded], float)
        lines.append(f"  (excluded ladders would have given median b = "
                     f"{np.median(eb):.3f} -- printed so the exclusion is visible)")
        lines.append("   NOTE: that median is not a measurement and must not be read")
        lines.append("   as evidence the gate is wrong. Fitting phi_J + c*p^b to one of")
        lines.append("   the actual excluded (non-monotonic) ladders was tested directly")
        lines.append("   in this session: it does not converge, and the returned b runs")
        lines.append("   to 9.7 / 10.1 / 13162 from starting guesses 0.3 / 1.0 / 2.0,")
        lines.append("   with the covariance unestimable. Individual values are")
        lines.append("   arbitrary, so a median over them is arbitrary too -- which is")
        lines.append("   the reason for excluding them, not an argument against it.")
    lines.append("")

    # ---- aggregate over seeds ------------------------------------------------
    # Also keyed on N: averaging phi_J over seeds must not average over system sizes,
    # and n_seeds would have counted the same seed twice at two N.
    agg: dict[tuple, list[dict]] = {}
    for r in rows:
        agg.setdefault((r["family"], r["cv_nom"], r["mu_p"], r["n_requested"]),
                       []).append(r)

    srows = []
    for (fam, cv, mu, n_req), g in sorted(agg.items()):
        pj = np.array([x["phi_J_lin"] for x in g], float)
        pj = pj[np.isfinite(pj)]
        if pj.size == 0:
            continue
        srows.append({
            "family": fam, "cv_nom": cv, "mu_p": mu, "n_requested": n_req,
            "cv_real": float(np.mean([x["cv_real"] for x in g])),
            "d32": float(np.mean([x["d32"] for x in g])),
            "n_seeds": int(pj.size),
            "phi_J": float(pj.mean()),
            "phi_J_sd": float(pj.std(ddof=1)) if pj.size > 1 else float("nan"),
        })
    fp = os.path.join(OUT, "phase2_summary.csv")
    with open(fp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=sorted(srows[0]))
        w.writeheader()
        w.writerows(srows)
    print(f"wrote {fp}")

    lines.append("phi_J by family, CV and friction (mean over seeds)")
    lines.append(f"{'family':>10} {'mu_p':>6} {'cv_nom':>7} {'cv_real':>8} "
                 f"{'d32':>7} {'n':>3} {'phi_J':>8} {'sd':>8}")
    for r in sorted(srows, key=lambda x: (x["family"], x["mu_p"], x["cv_nom"])):
        lines.append(f"{r['family']:>10} {r['mu_p']:6.2f} {r['cv_nom']:7.2f} "
                     f"{r['cv_real']:8.3f} {r['d32']:7.3f} {r['n_seeds']:3d} "
                     f"{r['phi_J']:8.4f} {r['phi_J_sd']:8.4f}")
    lines.append("")

    # ---- ANCHOR A1: frictionless monodisperse RCP ----------------------------
    a1 = [r for r in srows if r["cv_nom"] == 0.0 and r["mu_p"] == 0.0
          and r["family"] == "lognormal"]
    if a1:
        v = a1[0]["phi_J"]
        ok = abs(v - 0.64) <= 0.02
        lines.append(f"ANCHOR A1 (frictionless monodisperse RCP = 0.64 +/- 0.02): "
                     f"measured {v:.4f} -> {'PASS' if ok else 'FAIL'}")
    else:
        lines.append("ANCHOR A1: cv=0, mu=0 cell not present yet")

    # ---- ANCHOR A4: phi_J falls with friction --------------------------------
    # A4 IS GRADED PER SYSTEM SIZE. Before 2026-09-08 srows had one row per
    # (family, cv, mu_p); it now has one per N as well, and iterating it unfiltered put
    # two mu_p = 0.3 entries in a sequence that is supposed to be monotone IN mu_p --
    # which reported FAIL on a different question than the one A4 asks.
    for n_req in sorted({r["n_requested"] for r in srows}):
        mono = sorted([r for r in srows if r["cv_nom"] == 0.0
                       and r["family"] == "lognormal"
                       and r["n_requested"] == n_req],
                      key=lambda x: x["mu_p"])
        if len(mono) < 2:
            lines.append(f"ANCHOR A4 (N = {n_req}): only {len(mono)} mu_p value(s); "
                         f"not gradeable at this size")
            continue
        seq = [(r["mu_p"], r["phi_J"]) for r in mono]
        falling = all(seq[i + 1][1] <= seq[i][1] + 1e-9 for i in range(len(seq) - 1))
        lines.append(f"ANCHOR A4 (phi_J falls monotonically with mu_p), N = {n_req}: "
                     f"{[f'{m:.1f}:{p:.4f}' for m, p in seq]} -> "
                     f"{'PASS' if falling else 'FAIL'}")
    lines.append("")

    # ---- ANCHOR A5 + the M2 test --------------------------------------------
    # A5 sorts by cv_real and takes the FIRST and LAST row, so without N in the key it
    # could weld the low-CV row of one system size to the high-CV row of another. N is
    # part of the loop key for that reason, and every block below names its N.
    for fam, mu, n_req in sorted({(r["family"], r["mu_p"], r["n_requested"])
                                  for r in srows}):
        g = sorted([r for r in srows if r["mu_p"] == mu
                    and r["family"] == fam
                    and r["n_requested"] == n_req], key=lambda x: x["cv_real"])
        if len(g) < 2:
            continue
        lo, hi = g[0], g[-1]
        d_phiJ = hi["phi_J"] - lo["phi_J"]
        scatter = np.nanmax([lo.get("phi_J_sd", np.nan), hi.get("phi_J_sd", np.nan)])
        lines.append(f"--- family = {fam}, mu_p = {mu:.2f}, N = {n_req} ---")
        # Do NOT read two N blocks as a finite-size comparison unless they span the same
        # CV range. They need not: at 2026-09-08 the uniform family has N = 3000 cells
        # only at cv 0.4 and 0.5, while N = 1500 starts at 0.0, so the "same" anchor is
        # measuring a different interval at each size.
        _spans = {nn: (min(x["cv_nom"] for x in srows
                           if x["family"] == fam and x["mu_p"] == mu
                           and x["n_requested"] == nn),
                       max(x["cv_nom"] for x in srows
                           if x["family"] == fam and x["mu_p"] == mu
                           and x["n_requested"] == nn))
                  for nn in sorted({x["n_requested"] for x in srows
                                    if x["family"] == fam and x["mu_p"] == mu})}
        if len(_spans) > 1 and len(set(_spans.values())) > 1:
            lines.append("    !! CV SPAN DIFFERS BETWEEN SYSTEM SIZES: "
                         + ", ".join(f"N={nn}: cv {a:.2f}-{b:.2f}"
                                     for nn, (a, b) in sorted(_spans.items())))
            lines.append("    !! These blocks are NOT a matched finite-size comparison.")
        lines.append(f"ANCHOR A5 (phi_J RISES with polydispersity): "
                     f"cv {lo['cv_real']:.3f} -> {hi['cv_real']:.3f} gives phi_J "
                     f"{lo['phi_J']:.4f} -> {hi['phi_J']:.4f}, "
                     f"Delta = {d_phiJ:+.4f}")
        lines.append(f"    seed scatter (max sd) = {scatter:.4f}; "
                     f"Delta/scatter = {d_phiJ/scatter if scatter else float('nan'):.2f}")
        if abs(d_phiJ) < 2 * (scatter if np.isfinite(scatter) else 0):
            lines.append("    -> Delta phi_J is INSIDE 2x seed scatter: not resolved")

        lines.append("    M2 effect on eta ~ (phi_J - phi)^-2.")
        lines.append("    Ratio is eta(high CV) / eta(low CV). phi_J RISES with CV, so the")
        lines.append("    high-CV pack is FURTHER from jamming and the ratio is BELOW 1 --")
        lines.append("    i.e. polydispersity THINS the suspension, which is the direction")
        lines.append("    of the bench observation. The ratio runs to 0, never to infinity.")
        lines.append(f"      {'phi':>6} {'(phiJlo-phi)':>13} {'(phiJhi-phi)':>13} "
                     f"{'eta_hi/eta_lo':>14} {'decades':>8}")
        for phi in PHI_OF_INTEREST:
            a, bb = lo["phi_J"] - phi, hi["phi_J"] - phi
            if a <= 0 or bb <= 0:
                lines.append(f"      {phi:6.2f}   ABOVE phi_J -- already jammed")
                continue
            ratio = (a / bb) ** 2
            lines.append(f"      {phi:6.2f} {a:13.4f} {bb:13.4f} {ratio:14.3f} "
                         f"{np.log10(ratio):8.3f}")
        # The pre-registered bar is ONE DECADE of change. Because phi_J rises with CV,
        # that decade is a REDUCTION: solve (phiJlo-phi)^2/(phiJhi-phi)^2 = 1/10.
        #
        # The first version of this solved for a factor of +10 instead, and returned
        # phi = 0.6723 -- a value ABOVE BOTH phi_J's, which is unreachable by
        # construction. The sign of the pre-registered test was wrong, not the data.
        if d_phiJ > 0:
            k = np.sqrt(10.0)
            phi_needed = (k * lo["phi_J"] - hi["phi_J"]) / (k - 1.0)
            reachable = phi_needed < lo["phi_J"]
            lines.append(f"    phi at which the shift is worth ONE DECADE of thinning: "
                         f"{phi_needed:.4f}")
            lines.append(f"    (that is {'BELOW' if reachable else 'ABOVE'} "
                         f"phi_J(low CV) = {lo['phi_J']:.4f}, so it is "
                         f"{'reachable' if reachable else 'UNREACHABLE'})")
            at_lab = ((lo["phi_J"] - 0.44) / (hi["phi_J"] - 0.44)) ** 2
            lines.append(f"    at the lab's phi = 0.44 the effect is {at_lab:.3f}x "
                         f"({np.log10(at_lab):+.3f} decades)")
            lines.append("    -> M2 VERDICT: " + (
                "SUFFICIENT at the lab's phi" if abs(np.log10(at_lab)) >= 1.0 else
                f"INSUFFICIENT at the lab's phi = 0.44 -- the shift is real and "
                f"well-resolved, but worth only {at_lab:.2f}x there, not a decade. "
                f"A decade needs phi >= {phi_needed:.3f}."))
        else:
            lines.append("    -> M2 VERDICT: phi_J does not rise with CV here; "
                         "M2 cannot explain the observation in this direction")
        lines.append("")

    txt = "\n".join(lines)
    print(txt)
    fp = os.path.join(OUT, "phase2_verdict.txt")
    with open(fp, "w", encoding="utf-8") as fh:
        fh.write(txt + "\n")
    print(f"wrote {fp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
